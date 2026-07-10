"""Build identity-corrected, source-selected DEMI event-alignment evidence.

This script consumes the version-1 dual-source inventory rather than reading
EDF annotations directly. It is the bridge from the canonical identity and
event-source contract to the existing, historically grounded in-memory event
cleanup and timing-linkage audit.

The script deliberately reuses the pure sequence, trial-numbering, offset-join,
first-mismatch, and nearby-candidate functions in
02_audit_raw_event_sequences.py. It changes only the input evidence surface:

- EEG recording IDs are mapped to behavioural IDs through the tracked explicit
  identity contract;
- EEG 13 is joined to behavioural participant 11;
- EEG 11 remains explicitly unmapped and cannot join by numeric equality;
- coherent annotations remain preferred;
- EEG 94 uses its selected physical Trigger stream;
- raw source codes, comparison candidates, selected codes, and source reasons
  remain separate columns.

Inputs:

- _Data/eeg/event_source_inventory_v1/selected_event_evidence.csv;
- _Data/eeg/event_source_inventory_v1/event_source_file_inventory.csv;
- _Data/behavior/event_offsets/event_offsets_old_compatible.csv;
- analysis/eeg_mne/eeg_behavior_identity_contract.csv;
- existing annotation-only proposed join evidence for a local unaffected-case
  regression comparison when available.

Outputs below _Data/eeg/event_evidence_v1:

- selected_events_for_alignment.csv;
- proposed_selected_events.csv;
- proposed_trial_summary.csv;
- proposed_offset_join_audit.csv;
- first_mismatch_by_file.csv;
- nearby_trial_candidate_audit.csv;
- corrected_event_evidence_summary.md;
- corrected_event_evidence_run_manifest.json.

The old annotation-only directories are not overwritten. The old-compatible
behavioural offset CSV is read and hashed but never changed.

Safety boundaries:

- no raw EDF is opened or modified;
- no annotation or Trigger sample is mutated;
- no filtering, referencing, bad-channel work, ICA, CSD, or EEG derivative is
  created;
- proposed event labels remain separate from raw and selected labels;
- no final event-policy, inclusion, trial-rejection, or epoch decision is made.

Run from the repository root:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/06_build_corrected_event_evidence.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
from collections import Counter
from datetime import datetime
from pathlib import Path
from types import ModuleType
from typing import Any

import mne
import numpy as np
import pandas as pd

from event_source_contract import (
    IDENTITY_CONTRACT_VERSION,
    TASK_EVENT_CODES,
    load_identity_contract,
)


SOURCE_INVENTORY_DIR = (
    Path("_Data") / "eeg" / "event_source_inventory_v1"
)
SELECTED_SOURCE_EVENTS_PATH = (
    SOURCE_INVENTORY_DIR / "selected_event_evidence.csv"
)
SOURCE_FILE_INVENTORY_PATH = (
    SOURCE_INVENTORY_DIR / "event_source_file_inventory.csv"
)
OLD_COMPATIBLE_OFFSETS_PATH = (
    Path("_Data")
    / "behavior"
    / "event_offsets"
    / "event_offsets_old_compatible.csv"
)
OLD_ANNOTATION_ONLY_JOIN_PATH = (
    Path("_Data")
    / "eeg"
    / "event_sequence_audit"
    / "proposed_offset_join_audit.csv"
)
IDENTITY_CONTRACT_PATH = (
    Path("analysis") / "eeg_mne" / "eeg_behavior_identity_contract.csv"
)
SEQUENCE_SCRIPT_PATH = (
    Path("analysis") / "eeg_mne" / "02_audit_raw_event_sequences.py"
)
OUTPUT_DIR = Path("_Data") / "eeg" / "event_evidence_v1"

SELECTED_ALIGNMENT_EVENTS_FILENAME = "selected_events_for_alignment.csv"
PROPOSED_EVENTS_FILENAME = "proposed_selected_events.csv"
PROPOSED_TRIAL_SUMMARY_FILENAME = "proposed_trial_summary.csv"
PROPOSED_JOIN_FILENAME = "proposed_offset_join_audit.csv"
FIRST_MISMATCH_FILENAME = "first_mismatch_by_file.csv"
NEARBY_CANDIDATE_FILENAME = "nearby_trial_candidate_audit.csv"
SUMMARY_FILENAME = "corrected_event_evidence_summary.md"
RUN_MANIFEST_FILENAME = "corrected_event_evidence_run_manifest.json"

WATCHLIST_BEHAVIOURAL_IDS = {5, 11, 14, 54, 56, 65, 86, 89, 94, 100}
IDENTITY_OR_SOURCE_CHANGED_IDS = {11, 13, 94}

TASK_EVENT_NAMES = {
    "stim_on",
    "red_on",
    "trace_start",
    "trace_end",
    "accuracy_submit",
    "vividness_submit",
}

INTEGER_COLUMNS = {
    "eeg_recording_id",
    "behavioural_participant_id",
    "participant_id",
    "audit_participant_id",
    "offset_id",
    "event_row_order",
    "sample",
    "event_code",
    "raw_code",
    "normalized_candidate_code",
    "selected_normalized_code",
    "raw_trial_sequence",
    "old_update_trial_count",
    "join_trial_count",
    "proposed_raw_trial_sequence",
    "proposed_old_update_trial_count",
    "proposed_join_trial_count",
    "audit_trial_count",
    "offset_trial_count",
    "trial_count",
    "split_part",
}


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script path.

    Args:
        None.

    Returns:
        Absolute repository root.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def load_sequence_module(repo_root: Path) -> ModuleType:
    """Load script 02 as a module so its reviewed helpers can be reused.

    Args:
        repo_root: Repository root.

    Returns:
        Imported script-02 module.

    Side effects:
        Executes module-level definitions from script 02. Its main function is
        not called.
    """

    script_path = repo_root / SEQUENCE_SCRIPT_PATH
    spec = importlib.util.spec_from_file_location(
        "demi_raw_event_sequence_audit", script_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"unable to load event-sequence helpers from {script_path}"
        )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv_checked(
    path: Path,
    required_columns: set[str],
    label: str,
) -> pd.DataFrame:
    """Read a local CSV and validate required columns.

    Args:
        path: CSV path.
        required_columns: Columns that must be present.
        label: Human-readable source label.

    Returns:
        Loaded DataFrame.

    Side effects:
        Reads path. Raises RuntimeError for a missing file or schema.
    """

    if not path.exists():
        raise RuntimeError(f"required {label} is missing: {path}")
    frame = pd.read_csv(path, low_memory=False)
    missing = sorted(required_columns.difference(frame.columns))
    if missing:
        raise RuntimeError(
            f"{label} is missing required column(s): {', '.join(missing)}"
        )
    return coerce_integer_columns(frame)


def coerce_integer_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Coerce known integer-like columns to pandas nullable integers.

    Args:
        frame: Input DataFrame.

    Returns:
        Copy with known columns converted where present.

    Side effects:
        None.
    """

    out = frame.copy()
    for column in INTEGER_COLUMNS:
        if column in out.columns:
            out[column] = pd.to_numeric(
                out[column], errors="coerce"
            ).astype("Int64")
    return out


def bool_value(value: Any) -> bool:
    """Interpret bool-like CSV values conservatively.

    Args:
        value: Scalar CSV value.

    Returns:
        True only for explicit true-like values.

    Side effects:
        None.
    """

    if value is None or pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def participant_padded(value: Any) -> str:
    """Return a three-digit behavioural participant label.

    Args:
        value: Integer-like behavioural ID.

    Returns:
        Three-digit text or an empty string.

    Side effects:
        None.
    """

    if value is None or pd.isna(value):
        return ""
    return f"{int(value):03d}"


def selected_description(row: pd.Series) -> str:
    """Return a historical-compatible description for one selected event.

    Args:
        row: Selected event-source row.

    Returns:
        Original annotation text for annotation events, otherwise the selected
        physical raw code as text.

    Side effects:
        None.
    """

    if row["source_type"] == "annotation":
        return str(row["raw_value"])
    code = row.get("selected_normalized_code")
    return str(int(code)) if pd.notna(code) else str(row.get("raw_value", ""))


def build_selected_alignment_events(
    selected_events: pd.DataFrame,
    file_inventory: pd.DataFrame,
) -> pd.DataFrame:
    """Translate selected source rows into the script-02 input schema.

    Args:
        selected_events: Version-1 selected source-event rows.
        file_inventory: One-row-per-EDF source inventory.

    Returns:
        Alignment event table using behavioural participant IDs for all joins
        while retaining EEG/source provenance.

    Side effects:
        None. Raises RuntimeError if a selected event lacks an explicit
        behavioural mapping.
    """

    out = selected_events.copy()
    if out["behavioural_participant_id"].isna().any():
        bad = sorted(
            out.loc[
                out["behavioural_participant_id"].isna(),
                "eeg_recording_id",
            ]
            .dropna()
            .astype(int)
            .unique()
        )
        raise RuntimeError(
            "selected events cannot enter alignment without a behavioural "
            f"mapping; EEG ID(s): {', '.join(map(str, bad))}"
        )

    file_context = file_inventory[
        [
            "source_filename",
            "recording_duration_seconds",
            "source_discrepancy_class",
            "selected_source",
            "selection_status",
            "selection_reason",
            "selection_unresolved_issue",
            "annotation_event_count",
            "physical_event_count",
        ]
    ].drop_duplicates("source_filename")
    out = out.merge(
        file_context,
        on="source_filename",
        how="left",
        validate="many_to_one",
        suffixes=("", "_file"),
    )

    out["participant_id"] = out["behavioural_participant_id"].astype("Int64")
    out["participant_id_padded"] = out[
        "participant_id"
    ].map(participant_padded)
    out["watchlist_id"] = out["participant_id"].map(
        lambda value: (
            int(value) in WATCHLIST_BEHAVIOURAL_IDS
            if pd.notna(value)
            else False
        )
    )
    out["read_status"] = "ok"
    out["mne_orig_time"] = ""
    out["annotation_description"] = out.apply(
        selected_description, axis=1
    )
    out["event_code"] = out["selected_normalized_code"]
    out["event_name"] = out["selected_event_name"].fillna("")
    out["mapped_task_event"] = out["event_name"].isin(TASK_EVENT_NAMES)
    out["trial_numbering_candidate"] = out["mapped_task_event"]

    for column, default in [
        ("raw_trial_sequence", pd.NA),
        ("old_update_trial_count", pd.NA),
        ("join_trial_count", pd.NA),
        ("old_update_removed_reason", ""),
        ("fig_duration_seconds", np.nan),
        ("bad_start_like_event", False),
        ("practice_like_group", False),
        ("bad_start_like_group", False),
    ]:
        out[column] = default

    preferred = [
        "participant_id",
        "participant_id_padded",
        "eeg_recording_id",
        "eeg_recording_id_padded",
        "behavioural_participant_id",
        "identity_mapping_status",
        "identity_contract_version",
        "source_filename",
        "file_path",
        "file_role",
        "split_part",
        "watchlist_id",
        "read_status",
        "sampling_frequency_hz",
        "recording_duration_seconds",
        "session_date",
        "source_type",
        "selected_source",
        "selection_status",
        "selection_reason",
        "selection_unresolved_issue",
        "source_discrepancy_class",
        "event_row_order",
        "sample",
        "onset_seconds",
        "duration_seconds",
        "raw_value",
        "raw_code",
        "normalized_candidate_code",
        "normalization_candidate_rule",
        "normalization_candidate_applied",
        "selected_normalized_code",
        "normalization_applied",
        "annotation_description",
        "event_code",
        "event_name",
        "mapped_task_event",
        "trial_numbering_candidate",
        "raw_trial_sequence",
        "old_update_trial_count",
        "join_trial_count",
        "old_update_removed_reason",
        "fig_duration_seconds",
        "bad_start_like_event",
        "practice_like_group",
        "bad_start_like_group",
        "annotation_event_count",
        "physical_event_count",
    ]
    ordered = [column for column in preferred if column in out.columns] + [
        column for column in out.columns if column not in preferred
    ]
    return coerce_integer_columns(
        out[ordered].sort_values(
            [
                "behavioural_participant_id",
                "eeg_recording_id",
                "source_filename",
                "event_row_order",
            ],
            kind="mergesort",
        )
    ).reset_index(drop=True)


def selected_event_pattern(
    events: pd.DataFrame,
    source_filename: str,
) -> str:
    """Return compact selected semantic event counts for one EDF.

    Args:
        events: Selected alignment event table.
        source_filename: File to summarize.

    Returns:
        Semicolon-delimited task/special event counts.

    Side effects:
        None.
    """

    names = events.loc[
        events["source_filename"].eq(source_filename), "event_name"
    ]
    counts = Counter(names[names.fillna("").ne("")])
    ordered_names = [
        "stim_on",
        "red_on",
        "trace_start",
        "trace_end",
        "accuracy_submit",
        "vividness_submit",
        "file start",
    ]
    return ";".join(
        f"{name}={int(counts.get(name, 0))}"
        for name in ordered_names
        if counts.get(name, 0)
    )


def file_issue_codes(row: pd.Series) -> str:
    """Build factual selected-source file issue codes.

    Args:
        row: File inventory row with selected-source counts.

    Returns:
        Semicolon-delimited issue/context codes.

    Side effects:
        None.
    """

    codes: list[str] = []
    if row.get("read_status") != "ok":
        codes.append("mne_read_error")
    if row.get("identity_mapping_status") == "explicit_remap":
        codes.append("explicit_identity_remap")
    if row.get("identity_mapping_status") == "explicitly_unmapped":
        codes.append("explicitly_unmapped_eeg_recording")
    if not str(row.get("selected_source", "")).strip():
        codes.append("no_coherent_selected_event_source")
    if row.get("selected_source") == "physical_trigger":
        codes.append("physical_trigger_fallback_selected")
    if str(row.get("selection_unresolved_issue", "")).strip():
        codes.append(str(row["selection_unresolved_issue"]).strip())
    if row.get("file_role") == "split_part":
        codes.append("split_file_not_concatenated")
    if row.get("file_role") == "concatenated":
        codes.append("concatenated_file_not_repaired")
    return ";".join(dict.fromkeys(codes))


def build_file_summary(
    file_inventory: pd.DataFrame,
    selected_alignment_events: pd.DataFrame,
) -> pd.DataFrame:
    """Build the script-02 file-fact schema from source inventory rows.

    Args:
        file_inventory: One-row-per-EDF source inventory.
        selected_alignment_events: Selected events translated for alignment.

    Returns:
        One row per EDF. The compatibility field annotation_count contains the
        selected-source event count; explicit annotation and physical counts
        remain separate columns.

    Side effects:
        None.
    """

    out = file_inventory.copy()
    out["participant_id"] = out["behavioural_participant_id"].astype("Int64")
    out["participant_id_padded"] = out[
        "participant_id"
    ].map(participant_padded)
    out["watchlist_id"] = out["participant_id"].map(
        lambda value: (
            int(value) in WATCHLIST_BEHAVIOURAL_IDS
            if pd.notna(value)
            else False
        )
    )
    selected_counts = (
        selected_alignment_events.groupby("source_filename")
        .size()
        .rename("selected_event_count")
    )
    mapped_counts = (
        selected_alignment_events[
            selected_alignment_events["mapped_task_event"].map(bool_value)
        ]
        .groupby("source_filename")
        .size()
        .rename("selected_mapped_task_event_count")
    )
    out = out.merge(
        selected_counts,
        on="source_filename",
        how="left",
        validate="one_to_one",
    ).merge(
        mapped_counts,
        on="source_filename",
        how="left",
        validate="one_to_one",
    )
    out["selected_event_count"] = (
        out["selected_event_count"].fillna(0).astype(int)
    )
    out["selected_mapped_task_event_count"] = (
        out["selected_mapped_task_event_count"].fillna(0).astype(int)
    )

    # Script-02 helper names are retained for exact reuse, but their source is
    # explicit here to avoid implying that a physical fallback is an annotation.
    out["annotation_count"] = out["selected_event_count"]
    out["mapped_task_event_count"] = out[
        "selected_mapped_task_event_count"
    ]
    out["dominant_event_code_pattern"] = out["source_filename"].map(
        lambda name: selected_event_pattern(
            selected_alignment_events, str(name)
        )
    )
    out["file_issue_codes"] = out.apply(file_issue_codes, axis=1)

    preferred = [
        "participant_id",
        "participant_id_padded",
        "eeg_recording_id",
        "eeg_recording_id_padded",
        "behavioural_participant_id",
        "identity_mapping_status",
        "identity_contract_version",
        "source_filename",
        "file_path",
        "file_role",
        "split_part",
        "watchlist_id",
        "read_status",
        "sampling_frequency_hz",
        "recording_duration_seconds",
        "annotation_count",
        "mapped_task_event_count",
        "dominant_event_code_pattern",
        "file_issue_codes",
        "selected_event_count",
        "selected_mapped_task_event_count",
        "annotation_event_count",
        "physical_event_count",
        "source_discrepancy_class",
        "selected_source",
        "selection_status",
        "selection_reason",
        "selection_unresolved_issue",
        "normalization_applied",
    ]
    ordered = [column for column in preferred if column in out.columns] + [
        column for column in out.columns if column not in preferred
    ]
    return coerce_integer_columns(
        out[ordered].sort_values(
            ["eeg_recording_id", "source_filename"], kind="mergesort"
        )
    ).reset_index(drop=True)


def add_file_provenance(
    frame: pd.DataFrame,
    file_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Attach identity and selected-source provenance by source filename.

    Args:
        frame: Trial or join output with source_filename.
        file_summary: One-row-per-EDF file facts.

    Returns:
        Copy with contract provenance columns.

    Side effects:
        None.
    """

    context_columns = [
        "source_filename",
        "eeg_recording_id",
        "behavioural_participant_id",
        "identity_mapping_status",
        "identity_contract_version",
        "source_discrepancy_class",
        "selected_source",
        "selection_status",
        "selection_reason",
        "selection_unresolved_issue",
        "normalization_applied",
        "annotation_event_count",
        "physical_event_count",
        "selected_event_count",
    ]
    context = file_summary[
        [column for column in context_columns if column in file_summary.columns]
    ].drop_duplicates("source_filename")

    out = frame.copy()
    overlapping = [
        column
        for column in context.columns
        if column != "source_filename" and column in out.columns
    ]
    if overlapping:
        out = out.drop(columns=overlapping)
    return coerce_integer_columns(
        out.merge(
            context,
            on="source_filename",
            how="left",
            validate="many_to_one",
        )
    )


def add_join_identity_provenance(
    join: pd.DataFrame,
    file_summary: pd.DataFrame,
    identity_contract_path: Path,
) -> pd.DataFrame:
    """Attach file and inverse contract mapping to proposed join rows.

    Args:
        join: Proposed offset join audit.
        file_summary: One-row-per-EDF facts.
        identity_contract_path: Tracked identity CSV.

    Returns:
        Join audit with explicit EEG, behaviour, source, and mapping columns.

    Side effects:
        Reads the identity contract.
    """

    out = add_file_provenance(join, file_summary)
    mappings = load_identity_contract(identity_contract_path)
    inverse = {
        mapping.behavioural_participant_id: mapping
        for mapping in mappings.values()
        if mapping.behavioural_participant_id is not None
    }
    out["contract_eeg_recording_id"] = out["audit_participant_id"].map(
        lambda value: (
            inverse[int(value)].eeg_recording_id
            if pd.notna(value) and int(value) in inverse
            else pd.NA
        )
    )
    out["contract_mapping_status"] = out["audit_participant_id"].map(
        lambda value: (
            inverse[int(value)].mapping_status
            if pd.notna(value) and int(value) in inverse
            else ""
        )
    )
    out["identity_join_status"] = np.select(
        [
            out["source_filename"].notna()
            & out["eeg_recording_id"].notna()
            & out["audit_participant_id"].notna(),
            out["source_filename"].isna()
            & out["contract_eeg_recording_id"].notna(),
            out["source_filename"].isna()
            & out["contract_eeg_recording_id"].isna(),
        ],
        [
            "raw_event_row_uses_explicit_mapping",
            "offset_only_row_has_contract_eeg_mapping",
            "offset_only_row_without_raw_eeg_contract_mapping",
        ],
        default="unresolved_identity_join_state",
    )
    return coerce_integer_columns(out)


def sha256_file(path: Path) -> str:
    """Return a lowercase SHA-256 digest.

    Args:
        path: File to hash.

    Returns:
        Hexadecimal digest.

    Side effects:
        Reads path.
    """

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_regression_surface(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a stable unaffected-case join surface for comparison.

    Args:
        frame: Proposed join audit.

    Returns:
        Sorted string-normalized subset excluding the identity/source changed
        behavioural IDs 11, 13, and 94.

    Side effects:
        None.
    """

    columns = [
        "audit_participant_id",
        "audit_trial_count",
        "source_filename",
        "join_status",
        "old_trace_epoch_duration_mismatch",
        "old_stimulus_duration_mismatch",
        "missing_expected_event_names",
        "extra_expected_event_names",
        "alignment_problem_codes",
        "clean_timing_row",
    ]
    available = [column for column in columns if column in frame.columns]
    out = frame[
        ~frame["audit_participant_id"].isin(IDENTITY_OR_SOURCE_CHANGED_IDS)
    ][available].copy()
    for column in available:
        out[column] = out[column].fillna("").astype(str)
    return out.sort_values(available, kind="mergesort").reset_index(drop=True)


def unaffected_regression_result(
    old_join_path: Path,
    corrected_join: pd.DataFrame,
) -> dict[str, Any]:
    """Compare corrected and annotation-only joins outside changed IDs.

    Args:
        old_join_path: Existing annotation-only proposed join CSV.
        corrected_join: Corrected proposed join audit.

    Returns:
        Regression dictionary. If the old local evidence is absent, status is
        unavailable rather than guessed.

    Side effects:
        Reads old_join_path when present.
    """

    if not old_join_path.exists():
        return {
            "status": "unavailable",
            "old_rows": None,
            "corrected_rows": len(
                normalized_regression_surface(corrected_join)
            ),
            "surfaces_equal": None,
        }

    old_join = read_csv_checked(
        old_join_path,
        {"audit_participant_id", "audit_trial_count", "join_status"},
        "annotation-only proposed join audit",
    )
    old_surface = normalized_regression_surface(old_join)
    corrected_surface = normalized_regression_surface(corrected_join)
    return {
        "status": "compared",
        "old_rows": len(old_surface),
        "corrected_rows": len(corrected_surface),
        "surfaces_equal": old_surface.equals(corrected_surface),
    }


def join_metrics(frame: pd.DataFrame) -> dict[str, Any]:
    """Return compact aggregate join and strict-clean counts.

    Args:
        frame: Proposed join audit.

    Returns:
        JSON-serializable counts.

    Side effects:
        None.
    """

    status_counts = {
        str(label): int(count)
        for label, count in frame["join_status"].value_counts().items()
    }
    clean_column = (
        "strict_clean_timing"
        if "strict_clean_timing" in frame.columns
        else "clean_timing_row"
    )
    return {
        "rows": int(len(frame)),
        "join_status_counts": status_counts,
        "strict_clean_timing": int(
            frame.get(
                clean_column,
                pd.Series(False, index=frame.index),
            )
            .map(bool_value)
            .sum()
        ),
        "direct_alignment_problem_rows": int(
            frame.get(
                "alignment_problem_codes",
                pd.Series("", index=frame.index),
            )
            .fillna("")
            .ne("")
            .sum()
        ),
    }


def old_join_metrics(old_join_path: Path) -> dict[str, Any] | None:
    """Read aggregate metrics from the superseded annotation-only join.

    Args:
        old_join_path: Annotation-only join CSV.

    Returns:
        Metrics dictionary or None when unavailable.

    Side effects:
        Reads old_join_path when present.
    """

    if not old_join_path.exists():
        return None
    old_join = read_csv_checked(
        old_join_path,
        {"join_status"},
        "annotation-only proposed join audit",
    )
    return join_metrics(old_join)


def format_counts(series: pd.Series) -> list[str]:
    """Format sorted value counts as Markdown bullets.

    Args:
        series: Values to count.

    Returns:
        Markdown bullet lines.

    Side effects:
        None.
    """

    counts = series.fillna("").replace("", "none").value_counts().sort_index()
    return [f"- {label}: {int(count)}" for label, count in counts.items()]


def format_id_list(values: pd.Series) -> str:
    """Format unique integer-like IDs as sorted text.

    Args:
        values: ID values.

    Returns:
        Comma-delimited IDs or none.

    Side effects:
        None.
    """

    ids = sorted(
        {
            int(value)
            for value in values
            if value is not None and not pd.isna(value)
        }
    )
    return ", ".join(str(value) for value in ids) if ids else "none"


def id_join_summary(
    join: pd.DataFrame,
    behavioural_id: int,
) -> dict[str, int]:
    """Return join-status counts for one behavioural participant.

    Args:
        join: Corrected proposed join audit.
        behavioural_id: Behavioural participant ID.

    Returns:
        Count dictionary keyed by join status.

    Side effects:
        None.
    """

    rows = join[join["audit_participant_id"].eq(behavioural_id)]
    return {
        str(label): int(count)
        for label, count in rows["join_status"].value_counts().items()
    }


def build_summary(
    *,
    file_inventory: pd.DataFrame,
    selected_alignment_events: pd.DataFrame,
    proposed_events: pd.DataFrame,
    proposed_trials: pd.DataFrame,
    corrected_join: pd.DataFrame,
    old_metrics: dict[str, Any] | None,
    regression: dict[str, Any],
    started_at: datetime,
) -> str:
    """Build the corrected local event-evidence summary.

    Args:
        file_inventory: Source inventory file rows.
        selected_alignment_events: Selected source rows translated for join.
        proposed_events: Events after historical in-memory cleanup.
        proposed_trials: Proposed trial summary.
        corrected_join: Corrected offset join audit.
        old_metrics: Superseded annotation-only aggregate metrics.
        regression: Unaffected-case regression result.
        started_at: Run start time.

    Returns:
        Markdown summary.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    corrected_metrics = join_metrics(corrected_join)
    fallback = file_inventory[
        file_inventory["selected_source"].eq("physical_trigger")
    ]
    unresolved = file_inventory[
        file_inventory["selection_unresolved_issue"].fillna("").ne("")
    ]

    lines = [
        "# Identity-Corrected, Source-Selected Event Evidence",
        "",
        f"Generated: {finished_at.isoformat(timespec='seconds')}",
        f"Duration seconds: {(finished_at - started_at).total_seconds():.1f}",
        "",
        "## Safety boundary",
        "",
        "- No raw EDF was opened or modified by this script.",
        "- No source event was mutated.",
        "- Proposed cleanup labels remain separate from raw and selected labels.",
        "- No preprocessing, event-policy, inclusion, or epoch decision was made.",
        "",
        "## Identity contract",
        "",
        f"- Contract version: {IDENTITY_CONTRACT_VERSION}",
        "- EEG 13 maps explicitly to behavioural participant 11.",
        "- EEG 11 is explicitly unmapped and cannot join to behavioural 11.",
        "- All joins use behavioural_participant_id from the contract; there is no numeric fallback.",
        "",
        "## Source selection",
        "",
        *format_counts(file_inventory["selected_source"]),
        "",
        f"- Physical fallback files: {len(fallback)} "
        f"({format_id_list(fallback['eeg_recording_id'])})",
        f"- Unresolved source files: {len(unresolved)} "
        f"({format_id_list(unresolved['eeg_recording_id'])})",
        f"- Selected alignment events: {len(selected_alignment_events)}",
        f"- Proposed event rows: {len(proposed_events)}",
        f"- Proposed trial-summary rows: {len(proposed_trials)}",
        "",
        "## Corrected offset join",
        "",
        "Join status counts:",
        "",
    ]
    for label, count in sorted(
        corrected_metrics["join_status_counts"].items()
    ):
        lines.append(f"- {label}: {count}")
    lines.extend(
        [
            "",
            f"- Strict clean-timing rows: {corrected_metrics['strict_clean_timing']}",
            f"- Direct alignment-problem rows: {corrected_metrics['direct_alignment_problem_rows']}",
            "",
            "Behavioural ID 11 join counts:",
            "",
        ]
    )
    for label, count in sorted(id_join_summary(corrected_join, 11).items()):
        lines.append(f"- {label}: {count}")
    lines.extend(
        [
            "",
            "Behavioural ID 94 join counts:",
            "",
        ]
    )
    for label, count in sorted(id_join_summary(corrected_join, 94).items()):
        lines.append(f"- {label}: {count}")

    lines.extend(
        [
            "",
            "## Change from annotation-only evidence",
            "",
        ]
    )
    if old_metrics is None:
        lines.append("- Annotation-only aggregate metrics unavailable.")
    else:
        old_counts = old_metrics["join_status_counts"]
        all_statuses = sorted(
            set(old_counts).union(corrected_metrics["join_status_counts"])
        )
        for status in all_statuses:
            old_value = int(old_counts.get(status, 0))
            new_value = int(
                corrected_metrics["join_status_counts"].get(status, 0)
            )
            lines.append(
                f"- {status}: {old_value} -> {new_value} "
                f"(change {new_value - old_value:+d})"
            )
        lines.append(
            "- strict_clean_timing: "
            f"{old_metrics['strict_clean_timing']} -> "
            f"{corrected_metrics['strict_clean_timing']} "
            f"(change {corrected_metrics['strict_clean_timing'] - old_metrics['strict_clean_timing']:+d})"
        )

    lines.extend(
        [
            "",
            "## Unaffected-case regression",
            "",
            f"- Status: {regression['status']}",
            f"- Old rows: {regression['old_rows']}",
            f"- Corrected rows: {regression['corrected_rows']}",
            f"- Surfaces equal outside IDs 11, 13, and 94: {regression['surfaces_equal']}",
            "",
            "## Remaining boundary",
            "",
            "- These outputs are ready to inform a private event-policy table only after Tony reviews this remediation.",
            "- Split-file continuity, ID 5 file-start policy, partial recordings, and residual raw-only rows remain policy questions.",
            "",
        ]
    )
    return "\n".join(lines)


def build_manifest(
    *,
    repo_root: Path,
    selected_alignment_events: pd.DataFrame,
    proposed_events: pd.DataFrame,
    proposed_trials: pd.DataFrame,
    corrected_join: pd.DataFrame,
    old_metrics: dict[str, Any] | None,
    regression: dict[str, Any],
    started_at: datetime,
) -> dict[str, Any]:
    """Build run-level corrected-evidence provenance.

    Args:
        repo_root: Repository root.
        selected_alignment_events: Selected input events.
        proposed_events: Proposed cleanup events.
        proposed_trials: Trial summary.
        corrected_join: Corrected join audit.
        old_metrics: Annotation-only aggregate metrics.
        regression: Unaffected-case comparison.
        started_at: Run start time.

    Returns:
        JSON-serializable manifest.

    Side effects:
        Hashes the identity contract, selected-event CSV, and frozen offset CSV.
    """

    return {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "started": started_at.isoformat(timespec="seconds"),
        "script": "analysis/eeg_mne/06_build_corrected_event_evidence.py",
        "source_inventory_directory": SOURCE_INVENTORY_DIR.as_posix(),
        "selected_source_events": SELECTED_SOURCE_EVENTS_PATH.as_posix(),
        "selected_source_events_sha256": sha256_file(
            repo_root / SELECTED_SOURCE_EVENTS_PATH
        ),
        "identity_contract": IDENTITY_CONTRACT_PATH.as_posix(),
        "identity_contract_sha256": sha256_file(
            repo_root / IDENTITY_CONTRACT_PATH
        ),
        "identity_contract_version": IDENTITY_CONTRACT_VERSION,
        "old_compatible_offsets": OLD_COMPATIBLE_OFFSETS_PATH.as_posix(),
        "old_compatible_offsets_sha256": sha256_file(
            repo_root / OLD_COMPATIBLE_OFFSETS_PATH
        ),
        "output_directory": OUTPUT_DIR.as_posix(),
        "selected_alignment_event_rows": int(
            len(selected_alignment_events)
        ),
        "proposed_event_rows": int(len(proposed_events)),
        "proposed_trial_summary_rows": int(len(proposed_trials)),
        "corrected_join_metrics": join_metrics(corrected_join),
        "annotation_only_join_metrics": old_metrics,
        "unaffected_regression": regression,
        "python_version": platform.python_version(),
        "mne_version": mne.__version__,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "raw_eeg_opened": False,
        "raw_event_mutation_allowed": False,
        "preprocessing_allowed": False,
        "epoching_allowed": False,
    }


def validate_output_directory(output_dir: Path, repo_root: Path) -> None:
    """Require output below the versioned local EEG evidence path.

    Args:
        output_dir: Proposed output directory.
        repo_root: Repository root.

    Returns:
        None.

    Side effects:
        None. Raises RuntimeError for unsafe paths.
    """

    resolved = output_dir.resolve()
    required_parent = (repo_root / "_Data" / "eeg").resolve()
    try:
        resolved.relative_to(required_parent)
    except ValueError as error:
        raise RuntimeError(
            f"corrected event evidence must stay below {required_parent}"
        ) from error
    if resolved == (repo_root / "_Data" / "eeg" / "raw").resolve():
        raise RuntimeError("corrected event evidence cannot overwrite raw EEG")


def write_outputs(
    *,
    repo_root: Path,
    selected_alignment_events: pd.DataFrame,
    proposed_events: pd.DataFrame,
    proposed_trials: pd.DataFrame,
    proposed_join: pd.DataFrame,
    first_mismatch: pd.DataFrame,
    nearby_candidates: pd.DataFrame,
    summary: str,
    manifest: dict[str, Any],
) -> None:
    """Write versioned local corrected event evidence.

    Args:
        repo_root: Repository root.
        selected_alignment_events: Selected events in alignment schema.
        proposed_events: Events with in-memory cleanup columns.
        proposed_trials: Proposed trial summary.
        proposed_join: Corrected offset join audit.
        first_mismatch: First direct problem by file.
        nearby_candidates: Nearby timing candidates.
        summary: Markdown summary.
        manifest: Run manifest.

    Returns:
        None.

    Side effects:
        Writes ignored CSV, Markdown, and JSON files below OUTPUT_DIR.
    """

    output_dir = repo_root / OUTPUT_DIR
    validate_output_directory(output_dir, repo_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_alignment_events.to_csv(
        output_dir / SELECTED_ALIGNMENT_EVENTS_FILENAME, index=False
    )
    proposed_events.to_csv(
        output_dir / PROPOSED_EVENTS_FILENAME, index=False
    )
    proposed_trials.to_csv(
        output_dir / PROPOSED_TRIAL_SUMMARY_FILENAME, index=False
    )
    proposed_join.to_csv(
        output_dir / PROPOSED_JOIN_FILENAME, index=False
    )
    first_mismatch.to_csv(
        output_dir / FIRST_MISMATCH_FILENAME, index=False
    )
    nearby_candidates.to_csv(
        output_dir / NEARBY_CANDIDATE_FILENAME, index=False
    )
    (output_dir / SUMMARY_FILENAME).write_text(summary, encoding="utf-8")
    (output_dir / RUN_MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_corrected_evidence(
    repo_root: Path,
) -> dict[str, Any]:
    """Build all corrected evidence in memory.

    Args:
        repo_root: Repository root.

    Returns:
        Dictionary containing all output tables, summary, manifest, and
        regression metrics.

    Side effects:
        Reads local selected-source inventory, frozen offsets, tracked identity
        contract, and optional annotation-only join evidence. Does not write.
    """

    started_at = datetime.now()
    selected_source_events = read_csv_checked(
        repo_root / SELECTED_SOURCE_EVENTS_PATH,
        {
            "eeg_recording_id",
            "behavioural_participant_id",
            "identity_mapping_status",
            "source_filename",
            "source_type",
            "event_row_order",
            "sample",
            "onset_seconds",
            "raw_value",
            "raw_code",
            "normalized_candidate_code",
            "selected_normalized_code",
            "selected_event_name",
            "selected_source",
            "selection_reason",
            "normalization_applied",
        },
        "selected event-source evidence",
    )
    file_inventory = read_csv_checked(
        repo_root / SOURCE_FILE_INVENTORY_PATH,
        {
            "eeg_recording_id",
            "behavioural_participant_id",
            "identity_mapping_status",
            "source_filename",
            "file_role",
            "read_status",
            "sampling_frequency_hz",
            "recording_duration_seconds",
            "annotation_event_count",
            "physical_event_count",
            "source_discrepancy_class",
            "selected_source",
            "selection_status",
            "selection_reason",
            "selection_unresolved_issue",
            "normalization_applied",
        },
        "event-source file inventory",
    )

    selected_alignment_events = build_selected_alignment_events(
        selected_source_events, file_inventory
    )
    file_summary = build_file_summary(
        file_inventory, selected_alignment_events
    )
    sequence = load_sequence_module(repo_root)
    offsets = sequence.read_old_compatible_offsets(
        repo_root / OLD_COMPATIBLE_OFFSETS_PATH
    )

    proposed_events = sequence.propose_annotation_event_sequence(
        selected_alignment_events
    )
    proposed_events = sequence.assign_proposed_trial_numbers(
        proposed_events
    )
    proposed_trials = sequence.build_proposed_trial_summary(
        proposed_events, file_summary
    )
    proposed_trials = add_file_provenance(
        proposed_trials, file_summary
    )
    proposed_join = sequence.build_proposed_offset_join_audit(
        proposed_trials, offsets, file_summary
    )
    proposed_join = add_join_identity_provenance(
        proposed_join,
        file_summary,
        repo_root / IDENTITY_CONTRACT_PATH,
    )
    first_mismatch = sequence.build_first_mismatch_by_file(
        proposed_join, file_summary
    )
    first_mismatch = add_file_provenance(
        first_mismatch, file_summary
    )
    nearby_candidates = sequence.build_nearby_trial_candidate_audit(
        proposed_join, offsets
    )

    old_metrics = old_join_metrics(
        repo_root / OLD_ANNOTATION_ONLY_JOIN_PATH
    )
    regression = unaffected_regression_result(
        repo_root / OLD_ANNOTATION_ONLY_JOIN_PATH,
        proposed_join,
    )
    if (
        regression["status"] == "compared"
        and regression["surfaces_equal"] is not True
    ):
        raise RuntimeError(
            "unaffected proposed-join regression differs outside IDs "
            "11, 13, and 94"
        )

    summary = build_summary(
        file_inventory=file_inventory,
        selected_alignment_events=selected_alignment_events,
        proposed_events=proposed_events,
        proposed_trials=proposed_trials,
        corrected_join=proposed_join,
        old_metrics=old_metrics,
        regression=regression,
        started_at=started_at,
    )
    manifest = build_manifest(
        repo_root=repo_root,
        selected_alignment_events=selected_alignment_events,
        proposed_events=proposed_events,
        proposed_trials=proposed_trials,
        corrected_join=proposed_join,
        old_metrics=old_metrics,
        regression=regression,
        started_at=started_at,
    )
    return {
        "selected_alignment_events": selected_alignment_events,
        "proposed_events": proposed_events,
        "proposed_trials": proposed_trials,
        "proposed_join": proposed_join,
        "first_mismatch": first_mismatch,
        "nearby_candidates": nearby_candidates,
        "summary": summary,
        "manifest": manifest,
        "file_inventory": file_inventory,
        "regression": regression,
    }


def main() -> None:
    """Build and write the corrected version-1 event evidence.

    Args:
        None.

    Returns:
        None.

    Side effects:
        Reads ignored local event evidence and writes versioned ignored outputs
        below _Data/eeg/event_evidence_v1. Raw EEG is never opened.
    """

    repo_root = repo_root_from_script()
    try:
        result = run_corrected_evidence(repo_root)
        write_outputs(
            repo_root=repo_root,
            selected_alignment_events=result["selected_alignment_events"],
            proposed_events=result["proposed_events"],
            proposed_trials=result["proposed_trials"],
            proposed_join=result["proposed_join"],
            first_mismatch=result["first_mismatch"],
            nearby_candidates=result["nearby_candidates"],
            summary=result["summary"],
            manifest=result["manifest"],
        )
    except Exception as error:
        raise SystemExit(
            f"FAIL corrected event evidence: {error}"
        ) from error

    metrics = result["manifest"]["corrected_join_metrics"]
    print(
        f"Wrote corrected event evidence to {repo_root / OUTPUT_DIR}; "
        f"join rows={metrics['rows']}, "
        f"strict clean={metrics['strict_clean_timing']}"
    )


if __name__ == "__main__":
    main()
