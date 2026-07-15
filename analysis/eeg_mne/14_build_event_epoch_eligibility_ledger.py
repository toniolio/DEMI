"""Build the accepted DEMI event/epoch eligibility ledger, without epochs.

This is the deterministic policy-implementation stage between corrected event
evidence, completed continuous preprocessing, and any future authorization to
construct MNE epochs.  It consumes the accepted private event-policy v1.0
authority and the existing identity/source/alignment evidence.  It does not
reinterpret or rediscover event policy.

The canonical ledger has one row per row of the corrected proposed
trial/offset join surface.  Raw-present rows are linked back to their proposed
trial-event group and use the proposed ``trace_start`` event as the stable
ledger anchor when available.  A visibly labelled first-available-event anchor
is used only for malformed groups that lack ``trace_start``; it does not make
those rows eligible.  Offset-only rows remain explicit and have no invented
raw event.

Expected inputs, all resolved from the repository root:

- ``_Private/inventory/eeg_event_policy_v1.0.csv``;
- ``_Private/inventory/demi_source_id_crosswalk.csv``;
- ``analysis/eeg_mne/eeg_behavior_identity_contract.csv``;
- corrected event evidence below ``_Data/eeg/event_evidence_v1/``;
- selected-source inventory below ``_Data/eeg/event_source_inventory_v1/``;
- the reviewed raw-only classification below
  ``_Data/eeg/event_special_case_audit/``;
- the frozen old-compatible offset table below
  ``_Data/behavior/event_offsets/``;
- the finalized complete-surface run and per-recording manifests below
  ``_Data/eeg/mne_preprocessing/continuous_v2/``.

Generated local outputs below ``_Data/eeg/event_epoch_eligibility_v1/``:

- ``event_epoch_eligibility_ledger.parquet``: canonical durable table;
- ``event_epoch_eligibility_ledger.csv``: review surface;
- ``source_status.csv``: all continuous source files, including files with no
  candidate-event rows;
- count/reason tables in CSV;
- ``event_epoch_eligibility_summary.json`` and ``.md``;
- ``event_epoch_eligibility_run_manifest.json``.

The primary event-policy surface is exactly 8,905 rows.  The strict-clean-only
sensitivity surface is exactly 8,896 rows.  These remain separate from
continuous derivative and post-ICA availability.  In particular, files 49 and
54_1 retain non-excluding continuous QC warnings, while ID 86 retains its
accepted pre-ICA-only historical stop and is not silently analytically included
or globally excluded.

Safety boundaries:

- no EDF or FIF signal is opened;
- no raw annotation, Trigger sample, trial number, or offset is changed;
- no participant or analytic inclusion decision is made;
- no MNE Epochs, epoch FIF, AutoReject, CSD, time-frequency, feature, repaired
  EDF, or other signal derivative is created;
- all generated files are local/ignored tables and summaries only;
- outputs are written through same-directory temporary files and atomically
  replaced.

Run from the repository root:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/14_build_event_epoch_eligibility_ledger.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from event_epoch_eligibility import (
    EXPECTED_POSSIBLE_DISCREPANCY_COUNT,
    EXPECTED_PRIMARY_COUNT,
    EXPECTED_STRICT_CLEAN_COUNT,
    EXPECTED_WARNING_ONLY_COUNT,
    LEDGER_VERSION,
    POLICY_VERSION,
    SPLIT_EEG_IDS,
    assert_no_signal_outputs,
    bool_value,
    canonical_event_key,
    nullable_int,
    parse_continuous_recording_manifest,
    stable_frame_hash,
    status_fields,
    text_value,
    validate_ledger_contract,
)


POLICY_PATH = Path("_Private/inventory/eeg_event_policy_v1.0.csv")
PRIVATE_CROSSWALK_PATH = Path("_Private/inventory/demi_source_id_crosswalk.csv")
IDENTITY_CONTRACT_PATH = Path("analysis/eeg_mne/eeg_behavior_identity_contract.csv")
EVENT_EVIDENCE_DIR = Path("_Data/eeg/event_evidence_v1")
PROPOSED_JOIN_PATH = EVENT_EVIDENCE_DIR / "proposed_offset_join_audit.csv"
PROPOSED_EVENTS_PATH = EVENT_EVIDENCE_DIR / "proposed_selected_events.csv"
SELECTED_ALIGNMENT_EVENTS_PATH = EVENT_EVIDENCE_DIR / "selected_events_for_alignment.csv"
EVENT_EVIDENCE_MANIFEST_PATH = EVENT_EVIDENCE_DIR / "corrected_event_evidence_run_manifest.json"
SOURCE_FILE_INVENTORY_PATH = Path(
    "_Data/eeg/event_source_inventory_v1/event_source_file_inventory.csv"
)
RAW_ONLY_CLASSIFICATION_PATH = Path(
    "_Data/eeg/event_special_case_audit/raw_only_row_classification.csv"
)
OLD_COMPATIBLE_OFFSETS_PATH = Path(
    "_Data/behavior/event_offsets/event_offsets_old_compatible.csv"
)
CONTINUOUS_ROOT = Path("_Data/eeg/mne_preprocessing/continuous_v2")
OUTPUT_DIR = Path("_Data/eeg/event_epoch_eligibility_v1")

LEDGER_PARQUET_FILENAME = "event_epoch_eligibility_ledger.parquet"
LEDGER_CSV_FILENAME = "event_epoch_eligibility_ledger.csv"
SOURCE_STATUS_FILENAME = "source_status.csv"
COUNT_SUMMARY_FILENAME = "eligibility_count_summary.csv"
REASON_COUNTS_FILENAME = "eligibility_reason_counts.csv"
EVENT_POLICY_REASON_COUNTS_FILENAME = "event_policy_reason_counts.csv"
CONTINUOUS_STATUS_COUNTS_FILENAME = "continuous_status_counts.csv"
SUMMARY_JSON_FILENAME = "event_epoch_eligibility_summary.json"
SUMMARY_MD_FILENAME = "event_epoch_eligibility_summary.md"
RUN_MANIFEST_FILENAME = "event_epoch_eligibility_run_manifest.json"

MATCHED_STATUS = "raw_annotation_and_offset"
RAW_ONLY_STATUS = "raw_annotation_without_offset"
OFFSET_ONLY_STATUS = "offset_without_raw_annotation_trial"

ANCHOR_EVENT_COLUMNS = (
    "event_row_order",
    "sample",
    "onset_seconds",
    "duration_seconds",
    "source_type",
    "raw_value",
    "raw_code",
    "normalized_candidate_code",
    "normalization_candidate_rule",
    "normalization_candidate_applied",
    "selected_normalized_code",
    "normalization_applied",
    "selected_event_name",
    "raw_event_name",
    "raw_annotation_description",
    "proposed_event_name",
    "proposed_sequence_action",
    "proposed_sequence_reason",
    "proposed_label_differs_from_raw",
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script path.

    Returns:
        Absolute repository root.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for one input or output file.

    Args:
        path: File to hash.

    Returns:
        Lower-case SHA-256 hexadecimal digest.

    Side effects:
        Reads the file without modifying it.
    """

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv_checked(path: Path, required_columns: Iterable[str], label: str) -> pd.DataFrame:
    """Read a required CSV and fail clearly when its schema is incomplete.

    Args:
        path: CSV path.
        required_columns: Columns required by the ledger contract.
        label: Human-readable input label used in errors.

    Returns:
        Parsed data frame.

    Raises:
        FileNotFoundError: If the required input is absent.
        RuntimeError: If required columns are absent.

    Side effects:
        Reads one CSV.
    """

    if not path.is_file():
        raise FileNotFoundError(f"required {label} is missing: {path}")
    frame = pd.read_csv(path, low_memory=False)
    missing = sorted(set(required_columns) - set(frame.columns))
    if missing:
        raise RuntimeError(f"{label} is missing required columns: {missing}")
    return frame


def load_policy(policy_path: Path) -> dict[str, dict[str, Any]]:
    """Load and validate the accepted private event-policy table.

    Args:
        policy_path: Accepted event-policy v1.0 CSV.

    Returns:
        Dictionary keyed by policy ID.

    Raises:
        RuntimeError: If the version, accepted state, policy IDs, or expected
            row count differs from accepted v1.0.

    Side effects:
        Reads the private policy CSV only.
    """

    frame = read_csv_checked(
        policy_path,
        {
            "policy_id",
            "policy_version",
            "policy_state",
            "accepted_date",
            "proposed_first_pass_policy",
            "current_direct_evidence",
            "implementation_requirement",
        },
        "accepted event policy",
    )
    if len(frame) != 23:
        raise RuntimeError(f"accepted v1.0 policy must have 23 rows, found {len(frame)}")
    if set(frame["policy_version"].astype(str)) != {POLICY_VERSION}:
        raise RuntimeError("event-policy version differs from accepted v1.0")
    if set(frame["policy_state"].astype(str)) != {"accepted"}:
        raise RuntimeError("every event-policy row must be accepted")
    if frame["policy_id"].duplicated().any():
        raise RuntimeError("accepted event-policy IDs are not unique")
    expected_ids = {f"EP-{value:03d}" for value in range(1, 24)}
    observed_ids = set(frame["policy_id"].astype(str))
    if observed_ids != expected_ids:
        raise RuntimeError(
            f"accepted event-policy ID surface differs: {sorted(observed_ids ^ expected_ids)}"
        )
    return {
        str(row["policy_id"]): row.to_dict() for _, row in frame.iterrows()
    }


def load_identity_contract(path: Path) -> pd.DataFrame:
    """Load and validate the tracked explicit EEG/behaviour identity contract.

    Args:
        path: Tracked identity-contract CSV.

    Returns:
        Validated contract data frame.

    Raises:
        RuntimeError: If keys are duplicated or the 11/13 mapping contract is
            not exact.

    Side effects:
        Reads the tracked CSV.
    """

    frame = read_csv_checked(
        path,
        {
            "eeg_recording_id",
            "behavioural_participant_id",
            "mapping_status",
            "contract_version",
        },
        "tracked identity contract",
    )
    if frame["eeg_recording_id"].duplicated().any():
        raise RuntimeError("identity contract contains duplicate EEG IDs")
    eeg11 = frame[frame["eeg_recording_id"].eq(11)]
    eeg13 = frame[frame["eeg_recording_id"].eq(13)]
    if len(eeg11) != 1 or pd.notna(eeg11.iloc[0]["behavioural_participant_id"]):
        raise RuntimeError("EEG 11 must remain explicitly unavailable for behavioural joining")
    if (
        len(eeg13) != 1
        or nullable_int(eeg13.iloc[0]["behavioural_participant_id"]) != 11
        or str(eeg13.iloc[0]["mapping_status"]) != "explicit_remap"
    ):
        raise RuntimeError("behavioural participant 11 must map explicitly to EEG 13")
    return frame


def exploded_private_crosswalk(path: Path) -> pd.DataFrame:
    """Expand the reviewed private source crosswalk to one row per EEG file.

    Args:
        path: Private crosswalk CSV.

    Returns:
        File-level crosswalk retaining task, behaviour, figure, EEG, and
        mapping provenance fields.  Split parts remain separate filenames.

    Side effects:
        Reads the private crosswalk.
    """

    frame = read_csv_checked(
        path,
        {
            "task_database_id",
            "task_filename",
            "behavioural_participant_id",
            "figure_participant_id",
            "eeg_recording_id",
            "raw_eeg_file",
            "mapping_status",
            "pipeline_mapping_allowed",
        },
        "private source crosswalk",
    )
    rows: list[dict[str, Any]] = []
    for _, row in frame.iterrows():
        filenames = [
            item.strip()
            for item in str(row.get("raw_eeg_file", "")).split(";")
            if item.strip() and item.strip().lower() != "nan"
        ]
        if not filenames:
            continue
        for filename in filenames:
            record = row.to_dict()
            record["source_filename"] = filename
            rows.append(record)
    out = pd.DataFrame(rows)
    if out.empty or out["source_filename"].duplicated().any():
        raise RuntimeError("private crosswalk must resolve to unique EEG source filenames")
    eeg11 = out[out["eeg_recording_id"].eq(11)]
    eeg13 = out[out["eeg_recording_id"].eq(13)]
    if (
        len(eeg11) != 1
        or bool_value(eeg11.iloc[0]["pipeline_mapping_allowed"])
        or len(eeg13) != 1
        or nullable_int(eeg13.iloc[0]["behavioural_participant_id"]) != 11
    ):
        raise RuntimeError("private crosswalk disagrees with accepted EEG 11/13 identity authority")
    return out


def discover_finalized_continuous_run(continuous_root: Path, run_id: str | None) -> Path:
    """Select the finalized complete-surface continuous-v2 run manifest.

    Args:
        continuous_root: Continuous-v2 output root.
        run_id: Optional exact run directory name.  When absent, the
            lexicographically latest finalized 95-file verified run is used.

    Returns:
        Selected run-manifest path.

    Raises:
        RuntimeError: If no requested/qualifying run exists.

    Side effects:
        Reads candidate JSON manifests.
    """

    if run_id:
        candidates = [continuous_root / "runs" / run_id / "run_manifest.json"]
    else:
        candidates = sorted(
            (continuous_root / "runs").glob("*/run_manifest.json"), reverse=True
        )
    for path in candidates:
        if not path.is_file():
            continue
        manifest = json.loads(path.read_text(encoding="utf-8"))
        terminal = manifest.get("current_surface_terminal_counts", {})
        validation = manifest.get("validation", {})
        if (
            manifest.get("run_state") == "finalized"
            and len(manifest.get("recordings", [])) == 95
            and terminal == {
                "complete": 94,
                "failed": 0,
                "missing_or_incomplete": 0,
                "stopped": 1,
            }
            and bool(validation.get("all_saved_derivatives_reopened_and_valid"))
            and bool(validation.get("all_terminal_raw_sources_unchanged"))
            and validation.get("epochs_written") is False
            and validation.get("csd_written") is False
        ):
            return path
    requested = f" run {run_id}" if run_id else ""
    raise RuntimeError(f"no finalized verified 95-file continuous-v2{requested} manifest found")


def load_continuous_provenance(
    run_manifest_path: Path, repo_root: Path
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load flat file-level provenance from the selected continuous run.

    Args:
        run_manifest_path: Finalized continuous-v2 run manifest.
        repo_root: Repository root.

    Returns:
        Pair of a 95-row flat recording table and parsed run manifest.

    Raises:
        RuntimeError: If referenced manifests, source filenames, terminal
            counts, warning routes, or ID-86 route differ from authority.

    Side effects:
        Reads JSON manifests and checks canonical derivative file presence.  It
        never opens FIF content.
    """

    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for recording in run_manifest.get("recordings", []):
        manifest_path = Path(str(recording["manifest_path"]))
        if not manifest_path.is_absolute():
            manifest_path = repo_root / manifest_path
        if not manifest_path.is_file():
            raise RuntimeError(f"continuous recording manifest missing: {manifest_path}")
        parsed = json.loads(manifest_path.read_text(encoding="utf-8"))
        flat = parse_continuous_recording_manifest(parsed, manifest_path, repo_root)
        derivative_path = str(flat["continuous_derivative_path"])
        if derivative_path and not (repo_root / derivative_path).is_file():
            raise RuntimeError(f"canonical continuous derivative missing: {derivative_path}")
        rows.append(flat)
    frame = pd.DataFrame(rows).sort_values("source_filename", kind="mergesort").reset_index(drop=True)
    if len(frame) != 95 or frame["source_filename"].duplicated().any():
        raise RuntimeError("continuous-v2 provenance must contain 95 unique source files")
    for filename in ("demi_49 Data.edf", "demi_54_1 Data.edf"):
        row = frame[frame["source_filename"].eq(filename)]
        if (
            len(row) != 1
            or str(row.iloc[0]["continuous_v2_completion_class"]) != "complete_with_qc_warning"
            or "global_bad_proportion_above_25_percent"
            not in str(row.iloc[0]["continuous_qc_warning_codes"])
            or nullable_int(row.iloc[0]["interpolated_channel_count"]) != 8
        ):
            raise RuntimeError(f"continuous warning authority differs for {filename}")
    id86 = frame[frame["source_filename"].eq("demi_86 Data.edf")]
    if (
        len(id86) != 1
        or id86.iloc[0]["continuous_derivative_kind"] != "retained_pre_ica"
        or id86.iloc[0]["ica_terminal_route"]
        != "accepted_historical_zero_component_stop_pre_ica_retained"
        or nullable_int(id86.iloc[0]["ica_exclusion_count"]) != 0
    ):
        raise RuntimeError("ID 86 continuous authority differs from the accepted historical stop")
    return frame, run_manifest


def json_scalar(value: Any) -> Any:
    """Convert pandas/numpy scalars into deterministic JSON-safe values.

    Args:
        value: Scalar from an event row.

    Returns:
        Plain Python scalar or ``None`` for missing values.

    Side effects:
        None.
    """

    if value is None or value is pd.NA:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, np.generic):
        return value.item()
    return value


def build_trial_event_anchors(events: pd.DataFrame) -> pd.DataFrame:
    """Build one provenance anchor and event-list link per raw trial group.

    Args:
        events: Corrected proposed selected-event evidence.

    Returns:
        One row per behavioural ID/source/proposed raw-trial group.  The anchor
        is ``trace_start`` when present and otherwise the first available
        proposed task event.  ``trial_event_provenance_json`` retains every
        selected row in the group without changing raw/proposed fields.

    Side effects:
        None.
    """

    candidates = events[
        events["proposed_raw_trial_sequence"].notna()
        & events["proposed_trial_numbering_candidate"].map(bool_value)
    ].copy()
    group_columns = [
        "behavioural_participant_id",
        "source_filename",
        "proposed_raw_trial_sequence",
    ]
    rows: list[dict[str, Any]] = []
    for key, group in candidates.groupby(group_columns, sort=True, dropna=False):
        group = group.sort_values("event_row_order", kind="mergesort")
        trace_start = group[group["proposed_event_name"].eq("trace_start")]
        if not trace_start.empty:
            anchor = trace_start.iloc[0]
            anchor_role = "trace_start"
        else:
            usable = group[group["proposed_event_name"].fillna("").astype(str).ne("")]
            anchor = (usable if not usable.empty else group).iloc[0]
            anchor_role = "first_available_trial_event_fallback"

        provenance_rows = []
        for _, event_row in group.iterrows():
            provenance_rows.append(
                {
                    column: json_scalar(event_row.get(column))
                    for column in ANCHOR_EVENT_COLUMNS
                }
            )
        record: dict[str, Any] = {
            "behavioural_participant_id": nullable_int(key[0]),
            "source_filename": str(key[1]),
            "raw_trial_sequence": nullable_int(key[2]),
            "ledger_anchor_role": anchor_role,
            "trial_event_provenance_json": json.dumps(
                provenance_rows, sort_keys=True, separators=(",", ":")
            ),
        }
        for column in ANCHOR_EVENT_COLUMNS:
            record[f"anchor_{column}"] = json_scalar(anchor.get(column))
        rows.append(record)
    out = pd.DataFrame(rows)
    if out.duplicated(
        ["behavioural_participant_id", "source_filename", "raw_trial_sequence"]
    ).any():
        raise RuntimeError("trial-event anchor keys are not unique")
    return out


def corrected_raw_only_classification(
    join: pd.DataFrame, reviewed_classification: pd.DataFrame
) -> pd.DataFrame:
    """Port reviewed raw-only categories onto the corrected evidence surface.

    Args:
        join: Corrected proposed join audit.
        reviewed_classification: Existing special-case raw-only classification.

    Returns:
        Unique key/category table for corrected raw-only rows.  The reviewed 72
        possible-discrepancy keys are reused verbatim; newly corrected EEG-13
        and EEG-94 raw-only rows receive the already accepted frozen-target
        category rather than a new repair interpretation.

    Raises:
        RuntimeError: If the exact reviewed 72-row key set is not found in the
            corrected raw-only surface.

    Side effects:
        None.
    """

    key_columns = ["audit_participant_id", "source_filename", "raw_trial_sequence"]
    reviewed = reviewed_classification[
        key_columns
        + [
            "raw_only_category",
            "raw_only_category_reason",
            "best_candidate_offset_trial_count",
            "best_candidate_offset_delta",
            "best_candidate_timing_score_seconds",
            "best_candidate_clean_by_timing",
            "best_candidate_better_than_nominal",
        ]
    ].copy()
    possible = reviewed[
        reviewed["raw_only_category"].eq("possible_event_sequence_problem")
    ]
    if len(possible) != EXPECTED_POSSIBLE_DISCREPANCY_COUNT:
        raise RuntimeError("reviewed raw-only classification no longer contains exactly 72 discrepancy rows")

    raw_only = join[join["join_status"].eq(RAW_ONLY_STATUS)][key_columns].copy()
    merged = raw_only.merge(reviewed, on=key_columns, how="left", validate="one_to_one")
    missing_review = merged["raw_only_category"].isna()
    merged.loc[missing_review, "raw_only_category"] = "likely_old_compatible_row_surface_difference"
    merged.loc[missing_review, "raw_only_category_reason"] = (
        "identity/source-corrected raw event is absent from the frozen old-compatible offset target"
    )
    found_possible = merged["raw_only_category"].eq("possible_event_sequence_problem")
    if int(found_possible.sum()) != EXPECTED_POSSIBLE_DISCREPANCY_COUNT:
        raise RuntimeError("the reviewed 72 discrepancy keys do not survive on the corrected raw-only surface")
    return merged


def id5_segment(value: Any, file_start_onset: float) -> str:
    """Classify a row around ID 5's accepted ``file start`` boundary.

    Args:
        value: First available trial-event onset.
        file_start_onset: Accepted marker onset from selected event evidence.

    Returns:
        ``before_file_start``, ``after_file_start``, or ``unavailable``.

    Side effects:
        None.
    """

    try:
        onset = float(value)
    except (TypeError, ValueError):
        return "unavailable"
    if pd.isna(onset):
        return "unavailable"
    return "before_file_start" if onset < file_start_onset else "after_file_start"


def policy_ids_for_row(row: Mapping[str, Any]) -> tuple[list[str], str]:
    """Map one corrected join row to accepted v1.0 policy IDs.

    Args:
        row: Ledger row with join, case, warning, and classification fields.

    Returns:
        Ordered policy-ID list plus the most specific policy ID used for the
        exact accepted action/reason text.

    Side effects:
        None.
    """

    status = str(row.get("join_status", ""))
    behavioural_id = nullable_int(row.get("behavioural_participant_id"))
    eeg_id = nullable_int(row.get("eeg_recording_id"))
    file_role = str(row.get("file_role", ""))
    category = str(row.get("raw_only_category", ""))
    policy_ids = ["EP-014"]
    specific = "EP-015"

    if status == MATCHED_STATUS:
        policy_ids.extend(["EP-001", "EP-015"])
        if file_role == "split_part" or eeg_id in SPLIT_EEG_IDS:
            policy_ids.append("EP-010")
            specific = "EP-010"
        elif behavioural_id == 5:
            policy_ids.append("EP-009")
            specific = "EP-009"
        elif eeg_id == 94:
            policy_ids.append("EP-008")
            specific = "EP-008"
        elif behavioural_id == 11:
            policy_ids.append("EP-002")
            specific = "EP-002"
        if bool_value(row.get("duration_warning_only")):
            policy_ids.append("EP-022")
            specific = "EP-022"
    elif status == RAW_ONLY_STATUS:
        policy_ids.append("EP-016")
        specific = "EP-016"
        if behavioural_id in {5, 54, 56, 65}:
            policy_ids.append("EP-018")
            specific = "EP-018"
            if behavioural_id == 5:
                policy_ids.append("EP-009")
            else:
                policy_ids.append("EP-010")
        elif behavioural_id in {24, 89, 100}:
            policy_ids.append("EP-017")
            case_policy = {24: "EP-013", 89: "EP-011", 100: "EP-012"}[behavioural_id]
            policy_ids.append(case_policy)
            specific = case_policy
        elif category == "possible_event_sequence_problem":
            policy_ids.append("EP-019")
            specific = "EP-019"
        if behavioural_id == 89 and "missing" in str(row.get("trial_issue_codes", "")):
            policy_ids.append("EP-020")
            specific = "EP-020"
        if eeg_id == 94:
            policy_ids.append("EP-008")
            specific = "EP-008"
        if behavioural_id == 11:
            policy_ids.append("EP-002")
    elif status == OFFSET_ONLY_STATUS:
        policy_ids = ["EP-021"]
        specific = "EP-021"
        if behavioural_id in {6, 7}:
            policy_ids.append("EP-023")
            specific = "EP-023"
        elif behavioural_id in SPLIT_EEG_IDS:
            policy_ids.append("EP-010")
            specific = "EP-010"
        elif behavioural_id == 71:
            policy_ids.append("EP-019")
            specific = "EP-019"
    return list(dict.fromkeys(policy_ids)), specific


def event_policy_decision(row: Mapping[str, Any]) -> tuple[bool, str, str]:
    """Apply accepted v1.0 candidacy gates without continuous-status inference.

    Args:
        row: Corrected join row with explicit mapping/source/case fields.

    Returns:
        Boolean candidacy, machine-readable reason code, and readable reason.

    Side effects:
        None.
    """

    status = str(row.get("join_status", ""))
    behavioural_id = nullable_int(row.get("behavioural_participant_id"))
    eeg_id = nullable_int(row.get("eeg_recording_id"))
    file_role = str(row.get("file_role", ""))
    category = str(row.get("raw_only_category", ""))

    if status == MATCHED_STATUS:
        if file_role == "split_part" or eeg_id in SPLIT_EEG_IDS:
            return (
                False,
                "split_recording_excluded_first_pass",
                "Accepted v1.0 excludes split IDs 54, 56, and 65 from the first-pass candidate surface without merging or renumbering.",
            )
        if behavioural_id == 5 and str(row.get("id5_segment")) != "before_file_start":
            return (
                False,
                "id5_outside_accepted_pre_file_start_segment",
                "Accepted v1.0 permits only ID 5 matched target rows before the file-start marker.",
            )
        if eeg_id is None or not str(row.get("source_filename", "")).strip():
            return (
                False,
                "selected_raw_source_unavailable",
                "A matched candidate requires an explicit EEG source file.",
            )
        if not str(row.get("identity_mapping_status", "")).strip():
            return (
                False,
                "explicit_identity_mapping_unavailable",
                "A matched candidate requires an explicit reviewed identity mapping.",
            )
        return (
            True,
            "accepted_v1_primary_candidate",
            "Matched non-split row satisfies the accepted v1.0 first-pass event-policy gates.",
        )

    if status == RAW_ONLY_STATUS:
        if category == "possible_event_sequence_problem":
            return (
                False,
                "possible_event_sequence_discrepancy_unavailable",
                "Accepted v1.0 keeps the reviewed possible sequence-discrepancy row unavailable and defers repair or reassignment.",
            )
        if behavioural_id in {24, 89, 100}:
            return (
                False,
                "no_reviewed_old_compatible_target",
                "The accepted first-pass policy has no reviewed old-compatible target for this raw-only row.",
            )
        if behavioural_id in {5, 54, 56, 65}:
            return (
                False,
                "concatenated_or_split_raw_only_unavailable",
                "Accepted v1.0 keeps concatenated/split raw-only context outside the first-pass candidate surface.",
            )
        return (
            False,
            "no_frozen_offset_match",
            "Raw event presence alone cannot enter the frozen behavioural/tracing target without an accepted offset match.",
        )

    if status == OFFSET_ONLY_STATUS:
        if behavioural_id in {6, 7}:
            return (
                False,
                "no_raw_eeg_recording",
                "The frozen target row has no raw EEG recording under accepted v1.0.",
            )
        if behavioural_id in SPLIT_EEG_IDS:
            return (
                False,
                "split_offset_only_unavailable",
                "Split-file reconstruction is deferred and no raw event is assigned to this offset-only row.",
            )
        return (
            False,
            "no_selected_raw_event_match",
            "An offset-only row cannot become an epoch candidate without an accepted selected raw event.",
        )
    return False, "unrecognized_corrected_join_status", "The corrected join status is not recognized by the accepted ledger contract."


def add_continuous_statuses(row: Mapping[str, Any]) -> dict[str, Any]:
    """Create distinct derivative, post-ICA, and future-readiness statuses.

    Args:
        row: Ledger row after continuous provenance has been joined.

    Returns:
        Explicit status/reason families for continuous derivative availability,
        post-ICA availability, and future epoch-construction readiness.

    Side effects:
        None.
    """

    derivative_kind = str(row.get("continuous_derivative_kind", "") or "none")
    primary = bool_value(row.get("primary_ledger_eligibility"))
    available = derivative_kind in {"post_ica", "retained_pre_ica"}
    post_available = derivative_kind == "post_ica"
    eeg_id = nullable_int(row.get("eeg_recording_id"))
    status = str(row.get("join_status", ""))

    if available:
        continuous_true_code = (
            "post_ica_continuous_derivative_available"
            if post_available
            else "retained_pre_ica_continuous_derivative_available"
        )
        continuous_true_reason = (
            "Canonical post-ICA continuous FIF is available."
            if post_available
            else "Accepted ID-86 route retains a valid pre-ICA continuous FIF and ICA evidence."
        )
        continuous_false_code = ""
        continuous_false_reason = ""
    elif eeg_id is None:
        continuous_true_code = ""
        continuous_true_reason = ""
        continuous_false_code = "no_raw_eeg_recording_mapping"
        continuous_false_reason = "No EEG recording is mapped to this frozen offset row."
    elif status == OFFSET_ONLY_STATUS and eeg_id in SPLIT_EEG_IDS:
        continuous_true_code = ""
        continuous_true_reason = ""
        continuous_false_code = "split_offset_has_no_unique_continuous_file"
        continuous_false_reason = "The offset-only split row cannot be assigned to one continuous part without reconstruction."
    else:
        continuous_true_code = ""
        continuous_true_reason = ""
        continuous_false_code = "continuous_manifest_or_derivative_unavailable"
        continuous_false_reason = "No canonical continuous derivative is linked to this ledger row."

    out = status_fields(
        "continuous_derivative_availability",
        available,
        continuous_true_code,
        continuous_true_reason,
        continuous_false_code,
        continuous_false_reason,
    )

    if not post_available and derivative_kind == "retained_pre_ica":
        post_false_code = "accepted_id86_stop_no_post_ica_derivative"
        post_false_reason = "ID 86 retains pre-ICA signal and ICA evidence with zero exclusions; no post-ICA derivative is authorized."
    elif not post_available and not available:
        post_false_code = "continuous_derivative_unavailable"
        post_false_reason = "Post-ICA availability cannot be established without a linked continuous derivative."
    else:
        post_false_code = ""
        post_false_reason = ""
    out.update(
        status_fields(
            "post_ica_derivative_availability",
            post_available,
            "post_ica_continuous_derivative_available",
            "Canonical post-ICA continuous FIF is available.",
            post_false_code,
            post_false_reason,
        )
    )

    ready = primary and post_available
    if not primary:
        readiness_false_code = f"event_policy_{row.get('event_policy_reason_code', 'unavailable')}"
        readiness_false_reason = "The row is not in the accepted primary event-policy surface; continuous status does not change that result."
    elif not post_available and derivative_kind == "retained_pre_ica":
        readiness_false_code = "accepted_id86_stop_requires_later_epoch_route_decision"
        readiness_false_reason = "The event is an accepted policy candidate, but ID 86 has no post-ICA derivative; later epoch/analytic handling remains a Tony decision."
    elif not post_available:
        readiness_false_code = "post_ica_derivative_unavailable"
        readiness_false_reason = "The accepted event candidate lacks a post-ICA continuous derivative required by the ordinary future epoch route."
    else:
        readiness_false_code = ""
        readiness_false_reason = ""
    out.update(
        status_fields(
            "future_epoch_construction_readiness",
            ready,
            "accepted_event_candidate_with_post_ica_derivative",
            "Accepted primary event candidate and ordinary post-ICA continuous derivative are both available; epoch construction is still not authorized by this ledger.",
            readiness_false_code,
            readiness_false_reason,
        )
    )
    return out


def build_ledger(repo_root: Path, continuous_run_id: str | None = None) -> dict[str, Any]:
    """Build and validate all ledger tables in memory.

    Args:
        repo_root: Repository root.
        continuous_run_id: Optional exact continuous-v2 run ID.

    Returns:
        Dictionary containing canonical ledger, source status, summaries,
        reason tables, selected continuous manifest, policy, and input paths.

    Side effects:
        Reads accepted private/local/tracked tabular and JSON evidence.  It does
        not write outputs or open signal data.
    """

    paths = {
        "policy": repo_root / POLICY_PATH,
        "private_crosswalk": repo_root / PRIVATE_CROSSWALK_PATH,
        "identity_contract": repo_root / IDENTITY_CONTRACT_PATH,
        "proposed_join": repo_root / PROPOSED_JOIN_PATH,
        "proposed_events": repo_root / PROPOSED_EVENTS_PATH,
        "selected_alignment_events": repo_root / SELECTED_ALIGNMENT_EVENTS_PATH,
        "source_file_inventory": repo_root / SOURCE_FILE_INVENTORY_PATH,
        "raw_only_classification": repo_root / RAW_ONLY_CLASSIFICATION_PATH,
        "old_compatible_offsets": repo_root / OLD_COMPATIBLE_OFFSETS_PATH,
        "event_evidence_manifest": repo_root / EVENT_EVIDENCE_MANIFEST_PATH,
    }
    for label, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"required {label} input is missing: {path}")

    policy = load_policy(paths["policy"])
    identity = load_identity_contract(paths["identity_contract"])
    crosswalk = exploded_private_crosswalk(paths["private_crosswalk"])
    join = read_csv_checked(
        paths["proposed_join"],
        {
            "audit_participant_id",
            "audit_trial_count",
            "join_status",
            "alignment_problem_codes",
            "clean_timing_row",
            "source_filename",
            "file_role",
            "split_part",
            "raw_trial_sequence",
            "offset_trial_count",
            "trace_start_onset_seconds",
            "eeg_recording_id",
            "identity_mapping_status",
            "selected_source",
            "selection_reason",
            "identity_join_status",
        },
        "corrected proposed join audit",
    )
    events = read_csv_checked(
        paths["proposed_events"],
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
            "raw_event_name",
            "proposed_event_name",
            "proposed_sequence_action",
            "proposed_sequence_reason",
            "proposed_raw_trial_sequence",
            "proposed_trial_numbering_candidate",
            "selected_source",
        },
        "corrected proposed selected events",
    )
    selected_alignment = read_csv_checked(
        paths["selected_alignment_events"],
        {"eeg_recording_id", "source_filename", "source_type", "selected_source", "sample"},
        "selected alignment events",
    )
    source_inventory = read_csv_checked(
        paths["source_file_inventory"],
        {
            "eeg_recording_id",
            "behavioural_participant_id",
            "identity_mapping_status",
            "source_filename",
            "file_role",
            "split_part",
            "selected_source",
            "selection_status",
            "selection_reason",
        },
        "event-source file inventory",
    )
    reviewed_raw_only = read_csv_checked(
        paths["raw_only_classification"],
        {
            "audit_participant_id",
            "source_filename",
            "raw_trial_sequence",
            "raw_only_category",
            "raw_only_category_reason",
        },
        "reviewed raw-only classification",
    )

    continuous_run_path = discover_finalized_continuous_run(
        repo_root / CONTINUOUS_ROOT, continuous_run_id
    )
    continuous, continuous_run = load_continuous_provenance(
        continuous_run_path, repo_root
    )
    verification_paths = sorted((repo_root / CONTINUOUS_ROOT / "verifications").glob("*.json"))
    if not verification_paths:
        raise RuntimeError("continuous-v2 verification record is missing")
    continuous_verification_path = verification_paths[-1]
    verification = json.loads(continuous_verification_path.read_text(encoding="utf-8"))
    if not bool(verification.get("all_valid", verification.get("valid", True))):
        raise RuntimeError("latest continuous-v2 verification is not valid")

    expected_join_counts = {
        MATCHED_STATUS: 9_197,
        RAW_ONLY_STATUS: 1_313,
        OFFSET_ONLY_STATUS: 288,
    }
    if join["join_status"].value_counts().to_dict() != expected_join_counts:
        raise RuntimeError("corrected join surface counts differ from accepted event evidence")
    if len(join) != sum(expected_join_counts.values()):
        raise RuntimeError("corrected join row total is inconsistent")

    anchors = build_trial_event_anchors(events)
    classifications = corrected_raw_only_classification(join, reviewed_raw_only)
    ledger = join.copy()
    ledger["behavioural_participant_id"] = pd.to_numeric(
        ledger["audit_participant_id"], errors="coerce"
    ).astype("Int64")
    inverse_identity = {
        nullable_int(row["behavioural_participant_id"]): nullable_int(row["eeg_recording_id"])
        for _, row in identity.iterrows()
        if nullable_int(row["behavioural_participant_id"]) is not None
    }
    ledger["eeg_recording_id"] = pd.to_numeric(
        ledger["eeg_recording_id"], errors="coerce"
    ).astype("Int64")
    ledger["eeg_recording_id"] = ledger["eeg_recording_id"].fillna(
        ledger["behavioural_participant_id"].map(inverse_identity).astype("Int64")
    )
    ledger["identity_mapping_status"] = ledger["identity_mapping_status"].fillna("").astype(str)
    ledger.loc[
        ledger["identity_mapping_status"].str.strip().eq("")
        & ledger["behavioural_participant_id"].isin(inverse_identity),
        "identity_mapping_status",
    ] = ledger.loc[
        ledger["identity_mapping_status"].str.strip().eq("")
        & ledger["behavioural_participant_id"].isin(inverse_identity),
        "behavioural_participant_id",
    ].map(
        {
            nullable_int(row["behavioural_participant_id"]): str(row["mapping_status"])
            for _, row in identity.iterrows()
            if nullable_int(row["behavioural_participant_id"]) is not None
        }
    )
    ledger.loc[
        ledger["identity_mapping_status"].str.strip().eq("")
        & ledger["behavioural_participant_id"].isin([6, 7]),
        "identity_mapping_status",
    ] = "no_raw_eeg_recording"

    ledger = ledger.merge(
        anchors,
        on=["behavioural_participant_id", "source_filename", "raw_trial_sequence"],
        how="left",
        validate="many_to_one",
    )
    raw_present = ledger["source_filename"].notna()
    if ledger.loc[raw_present, "ledger_anchor_role"].isna().any():
        missing = ledger.loc[raw_present & ledger["ledger_anchor_role"].isna()]
        raise RuntimeError(f"{len(missing)} raw-present join rows lack trial-event provenance anchors")
    ledger = ledger.merge(
        classifications,
        on=["audit_participant_id", "source_filename", "raw_trial_sequence"],
        how="left",
        validate="many_to_one",
    )

    file_start_rows = events[
        events["eeg_recording_id"].eq(5)
        & events["raw_annotation_description"].fillna("").astype(str).str.lower().eq("file start")
    ]
    if len(file_start_rows) != 1:
        raise RuntimeError("ID 5 must have exactly one selected file-start marker")
    file_start_onset = float(file_start_rows.iloc[0]["onset_seconds"])
    ledger["id5_file_start_onset_seconds"] = np.where(
        ledger["behavioural_participant_id"].eq(5), file_start_onset, np.nan
    )
    first_onset = ledger["stim_on_onset_seconds"].fillna(ledger["anchor_onset_seconds"])
    ledger["id5_segment"] = "not_id5"
    id5_rows = ledger["behavioural_participant_id"].eq(5)
    ledger.loc[id5_rows, "id5_segment"] = first_onset[id5_rows].map(
        lambda value: id5_segment(value, file_start_onset)
    )

    alignment_codes = ledger["alignment_problem_codes"].fillna("").astype(str)
    ledger["duration_warning_only"] = (
        ledger["join_status"].eq(MATCHED_STATUS)
        & ~ledger["clean_timing_row"].map(bool_value)
        & alignment_codes.eq("old_trace_epoch_duration_mismatch")
        & ~ledger["file_role"].eq("split_part")
    )

    policy_ids_values: list[str] = []
    policy_specific_values: list[str] = []
    policy_action_values: list[str] = []
    policy_reason_values: list[str] = []
    candidate_status_rows: list[dict[str, Any]] = []
    event_reason_codes: list[str] = []
    event_reasons: list[str] = []
    strict_status_rows: list[dict[str, Any]] = []
    for _, row in ledger.iterrows():
        policy_ids, specific = policy_ids_for_row(row)
        candidate, event_code, event_reason = event_policy_decision(row)
        policy_ids_values.append(";".join(policy_ids))
        policy_specific_values.append(specific)
        policy_action_values.append(str(policy[specific]["proposed_first_pass_policy"]))
        policy_reason_values.append(str(policy[specific]["current_direct_evidence"]))
        event_reason_codes.append(event_code)
        event_reasons.append(event_reason)
        candidate_status_rows.append(
            status_fields(
                "accepted_event_policy_candidate",
                candidate,
                event_code,
                event_reason,
                event_code,
                event_reason,
            )
        )
        strict = candidate and bool_value(row.get("clean_timing_row"))
        if strict:
            strict_code = "accepted_v1_strict_clean_candidate"
            strict_reason = "Accepted primary candidate also satisfies the strict-clean timing definition."
        elif candidate and bool_value(row.get("duration_warning_only")):
            strict_code = "duration_warning_excluded_strict_clean_sensitivity"
            strict_reason = "Accepted duration-warning row remains primary-eligible but is omitted from the required strict-clean-only sensitivity surface."
        else:
            strict_code = f"event_policy_{event_code}"
            strict_reason = "The row is not a primary event-policy candidate and therefore cannot enter the strict-clean sensitivity surface."
        strict_status_rows.append(
            status_fields(
                "strict_clean_only_sensitivity_eligibility",
                strict,
                strict_code,
                strict_reason,
                strict_code,
                strict_reason,
            )
        )

    ledger["accepted_policy_ids"] = policy_ids_values
    ledger["accepted_specific_policy_id"] = policy_specific_values
    ledger["accepted_policy_action"] = policy_action_values
    ledger["accepted_policy_reason"] = policy_reason_values
    ledger["event_policy_reason_code"] = event_reason_codes
    ledger["event_policy_reason"] = event_reasons
    ledger = pd.concat(
        [ledger.reset_index(drop=True), pd.DataFrame(candidate_status_rows), pd.DataFrame(strict_status_rows)],
        axis=1,
    )
    ledger["primary_ledger_eligibility"] = ledger["accepted_event_policy_candidate"]
    ledger["primary_ledger_eligibility_status"] = ledger["accepted_event_policy_candidate_status"]
    ledger["primary_ledger_eligibility_reason_code"] = ledger["accepted_event_policy_candidate_reason_code"]
    ledger["primary_ledger_eligibility_reason"] = ledger["accepted_event_policy_candidate_reason"]

    source_to_eeg = source_inventory.set_index("source_filename")["eeg_recording_id"].map(nullable_int).to_dict()
    eeg_to_sources: dict[int, list[str]] = {}
    for source_filename, eeg_id in source_to_eeg.items():
        if eeg_id is not None:
            eeg_to_sources.setdefault(eeg_id, []).append(str(source_filename))
    ledger["continuous_join_source_filename"] = ledger["source_filename"].fillna("").astype(str)
    offset_only_missing_source = ledger["continuous_join_source_filename"].str.strip().eq("")
    for index in ledger.index[offset_only_missing_source]:
        eeg_id = nullable_int(ledger.at[index, "eeg_recording_id"])
        sources = eeg_to_sources.get(eeg_id, []) if eeg_id is not None else []
        if len(sources) == 1:
            ledger.at[index, "continuous_join_source_filename"] = sources[0]
    ledger = ledger.merge(
        continuous,
        left_on="continuous_join_source_filename",
        right_on="source_filename",
        how="left",
        validate="many_to_one",
        suffixes=("", "_continuous"),
    )
    ledger["continuous_derivative_kind"] = ledger["continuous_derivative_kind"].fillna("none")
    continuous_status_rows = [add_continuous_statuses(row) for _, row in ledger.iterrows()]
    ledger = pd.concat(
        [ledger.reset_index(drop=True), pd.DataFrame(continuous_status_rows)], axis=1
    )

    ledger["accepted_raw_event_label"] = ledger["anchor_raw_event_name"]
    ledger["accepted_proposed_event_label"] = ledger["anchor_proposed_event_name"]
    ledger["source_type"] = ledger["anchor_source_type"]
    ledger["proposed_cleanup_action"] = ledger["anchor_proposed_sequence_action"]
    ledger["proposed_cleanup_reason"] = ledger["anchor_proposed_sequence_reason"]
    ledger["canonical_event_key_basis"] = ledger.apply(
        lambda row: (
            f"raw|eeg={nullable_int(row.get('eeg_recording_id'))}|beh={nullable_int(row.get('behavioural_participant_id'))}|"
            f"file={text_value(row.get('source_filename'))}|raw_trial={nullable_int(row.get('raw_trial_sequence'))}|"
            f"sample={nullable_int(row.get('anchor_sample'))}"
            if text_value(row.get("source_filename"))
            else f"offset_only|eeg={nullable_int(row.get('eeg_recording_id'))}|beh={nullable_int(row.get('behavioural_participant_id'))}|offset_trial={nullable_int(row.get('offset_trial_count'))}"
        ),
        axis=1,
    )
    ledger["canonical_event_key"] = ledger.apply(canonical_event_key, axis=1)
    ledger["ledger_version"] = LEDGER_VERSION
    ledger["event_policy_version"] = POLICY_VERSION
    ledger["continuous_v2_run_id"] = str(continuous_run["run_id"])

    crosswalk_columns = [
        "source_filename",
        "task_database_id",
        "task_filename",
        "figure_participant_id",
        "mapping_confidence",
        "pipeline_mapping_allowed",
    ]
    ledger = ledger.merge(
        crosswalk[crosswalk_columns],
        on="source_filename",
        how="left",
        validate="many_to_one",
    )

    preferred_columns = [
        "ledger_version",
        "event_policy_version",
        "canonical_event_key",
        "canonical_event_key_basis",
        "behavioural_participant_id",
        "task_database_id",
        "task_filename",
        "figure_participant_id",
        "eeg_recording_id",
        "identity_mapping_status",
        "mapping_confidence",
        "pipeline_mapping_allowed",
        "identity_join_status",
        "source_filename",
        "file_path",
        "file_role",
        "split_part",
        "source_type",
        "selected_source",
        "selection_status",
        "selection_reason",
        "ledger_anchor_role",
        "anchor_event_row_order",
        "anchor_sample",
        "anchor_onset_seconds",
        "anchor_duration_seconds",
        "anchor_raw_value",
        "anchor_raw_code",
        "accepted_raw_event_label",
        "accepted_proposed_event_label",
        "proposed_cleanup_action",
        "proposed_cleanup_reason",
        "trial_event_provenance_json",
        "audit_trial_count",
        "raw_trial_sequence",
        "old_update_trial_count",
        "join_trial_count",
        "offset_trial_count",
        "join_status",
        "join_audit_issue_codes",
        "alignment_problem_codes",
        "clean_timing_row",
        "duration_warning_only",
        "stim_on_onset_seconds",
        "red_on_onset_seconds",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
        "trace_epoch_duration_seconds",
        "trace_epoch_end_difference_seconds",
        "stimulus_duration_difference_seconds",
        "old_trace_epoch_duration_mismatch",
        "old_stimulus_duration_mismatch",
        "raw_only_category",
        "raw_only_category_reason",
        "id5_segment",
        "id5_file_start_onset_seconds",
        "accepted_policy_ids",
        "accepted_specific_policy_id",
        "accepted_policy_action",
        "accepted_policy_reason",
        "event_policy_reason_code",
        "event_policy_reason",
        "accepted_event_policy_candidate",
        "accepted_event_policy_candidate_status",
        "accepted_event_policy_candidate_reason_code",
        "accepted_event_policy_candidate_reason",
        "primary_ledger_eligibility",
        "primary_ledger_eligibility_status",
        "primary_ledger_eligibility_reason_code",
        "primary_ledger_eligibility_reason",
        "strict_clean_only_sensitivity_eligibility",
        "strict_clean_only_sensitivity_eligibility_status",
        "strict_clean_only_sensitivity_eligibility_reason_code",
        "strict_clean_only_sensitivity_eligibility_reason",
        "continuous_v2_run_id",
        "continuous_v2_terminal_status",
        "continuous_v2_completion_class",
        "continuous_v2_manifest_path",
        "continuous_v2_manifest_id",
        "continuous_derivative_kind",
        "continuous_derivative_path",
        "continuous_qc_warning_codes",
        "continuous_qc_warnings_json",
        "interpolated_channel_count",
        "interpolated_channel_proportion",
        "interpolation_proportion_denominator",
        "ica_terminal_route",
        "ica_proposal_count",
        "ica_exclusion_count",
        "ica_application_performed",
        "continuous_stop_reason_code",
        "continuous_stop_reason",
        "continuous_derivative_availability",
        "continuous_derivative_availability_status",
        "continuous_derivative_availability_reason_code",
        "continuous_derivative_availability_reason",
        "post_ica_derivative_availability",
        "post_ica_derivative_availability_status",
        "post_ica_derivative_availability_reason_code",
        "post_ica_derivative_availability_reason",
        "future_epoch_construction_readiness",
        "future_epoch_construction_readiness_status",
        "future_epoch_construction_readiness_reason_code",
        "future_epoch_construction_readiness_reason",
    ]
    ordered = [column for column in preferred_columns if column in ledger.columns]
    ordered.extend(column for column in ledger.columns if column not in ordered and not column.endswith("_continuous"))
    ledger = ledger[ordered].sort_values(
        [
            "behavioural_participant_id",
            "eeg_recording_id",
            "source_filename",
            "raw_trial_sequence",
            "offset_trial_count",
        ],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)

    counts = validate_ledger_contract(ledger)
    validate_project_specific_contracts(ledger, selected_alignment, source_inventory)

    source_status = build_source_status(source_inventory, crosswalk, continuous, ledger)
    reason_counts = build_reason_counts(ledger)
    policy_reason_counts = build_event_policy_reason_counts(ledger)
    continuous_status_counts = build_continuous_status_counts(ledger)
    count_summary = build_count_summary(ledger)

    input_paths = dict(paths)
    input_paths["continuous_run_manifest"] = continuous_run_path
    input_paths["continuous_verification"] = continuous_verification_path
    input_hashes = {
        label: {
            "path": path.resolve().relative_to(repo_root.resolve()).as_posix(),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for label, path in sorted(input_paths.items())
    }
    summary = build_summary_payload(
        ledger, source_status, counts, continuous_run, input_hashes
    )
    return {
        "ledger": ledger,
        "source_status": source_status,
        "count_summary": count_summary,
        "reason_counts": reason_counts,
        "event_policy_reason_counts": policy_reason_counts,
        "continuous_status_counts": continuous_status_counts,
        "summary": summary,
        "policy": policy,
        "continuous_run_manifest": continuous_run,
        "continuous_run_manifest_path": continuous_run_path,
        "continuous_verification_path": continuous_verification_path,
        "input_hashes": input_hashes,
    }


def validate_project_specific_contracts(
    ledger: pd.DataFrame,
    selected_alignment: pd.DataFrame,
    source_inventory: pd.DataFrame,
) -> None:
    """Validate accepted DEMI exceptional-case behavior on the built ledger.

    Args:
        ledger: Canonical ledger after all statuses are assigned.
        selected_alignment: Selected event evidence used for source checks.
        source_inventory: File-level selected-source inventory.

    Raises:
        RuntimeError: If an accepted high-risk special-case contract differs.

    Side effects:
        None.
    """

    primary = ledger["primary_ledger_eligibility"].map(bool_value)
    strict = ledger["strict_clean_only_sensitivity_eligibility"].map(bool_value)
    if ledger.loc[primary, "eeg_recording_id"].isin(SPLIT_EEG_IDS).any():
        raise RuntimeError("split IDs entered the first-pass primary surface")
    id5 = ledger[primary & ledger["behavioural_participant_id"].eq(5)]
    if len(id5) != 114 or set(id5["id5_segment"]) != {"before_file_start"}:
        raise RuntimeError("ID 5 primary surface differs from the accepted 114-row pre-file-start segment")
    mapped11 = ledger[primary & ledger["behavioural_participant_id"].eq(11)]
    if len(mapped11) != 99 or set(mapped11["eeg_recording_id"].dropna().astype(int)) != {13}:
        raise RuntimeError("behavioural 11 primary rows do not map exactly to EEG 13")
    if ledger["eeg_recording_id"].eq(11).any():
        raise RuntimeError("EEG 11 must not enter the event-row ledger through numeric equality")
    id94 = ledger[primary & ledger["eeg_recording_id"].eq(94)]
    if len(id94) != 104 or set(id94["source_type"]) != {"physical_trigger"}:
        raise RuntimeError("EEG 94 primary surface must contain 104 physical-Trigger rows")
    selected94 = selected_alignment[selected_alignment["eeg_recording_id"].eq(94)]
    if len(selected94) != 617 or set(selected94["source_type"]) != {"physical_trigger"}:
        raise RuntimeError("EEG 94 selected evidence no longer contains the validated 617 physical events")
    inventory94 = source_inventory[source_inventory["eeg_recording_id"].eq(94)]
    if len(inventory94) != 1 or inventory94.iloc[0]["selected_source"] != "physical_trigger":
        raise RuntimeError("EEG 94 source inventory no longer selects physical Trigger")
    warnings = ledger[primary & ~strict]
    if len(warnings) != 9 or not warnings["duration_warning_only"].map(bool_value).all():
        raise RuntimeError("primary/strict difference is not exactly the nine accepted duration warnings")
    for filename in ("demi_49 Data.edf", "demi_54_1 Data.edf"):
        rows = ledger[ledger["source_filename"].eq(filename)]
        if rows.empty or not rows["continuous_qc_warning_codes"].str.contains(
            "global_bad_proportion_above_25_percent", na=False
        ).all():
            raise RuntimeError(f"continuous warning was not propagated for {filename}")
        if filename == "demi_49 Data.edf" and not rows.loc[
            rows["primary_ledger_eligibility"].map(bool_value),
            "future_epoch_construction_readiness",
        ].map(bool_value).all():
            raise RuntimeError("file 49 warning incorrectly blocks future readiness")
    id86 = ledger[ledger["eeg_recording_id"].eq(86)]
    if (
        id86.empty
        or not id86["continuous_derivative_availability"].map(bool_value).all()
        or id86["post_ica_derivative_availability"].map(bool_value).any()
        or id86["future_epoch_construction_readiness"].map(bool_value).any()
        or set(id86["continuous_derivative_kind"]) != {"retained_pre_ica"}
    ):
        raise RuntimeError("ID 86 continuous/post-ICA/readiness distinction is not explicit")


def build_source_status(
    source_inventory: pd.DataFrame,
    crosswalk: pd.DataFrame,
    continuous: pd.DataFrame,
    ledger: pd.DataFrame,
) -> pd.DataFrame:
    """Build a 95-file source-status surface including zero-ledger-row files.

    Args:
        source_inventory: Current selected-source file inventory.
        crosswalk: Exploded private source crosswalk.
        continuous: Parsed continuous-v2 recording provenance.
        ledger: Canonical row-level ledger.

    Returns:
        One row per source recording with identity/source/continuous status and
        ledger-row counts.  EEG 11 remains explicitly unmapped here even though
        it correctly has no canonical event row.

    Side effects:
        None.
    """

    counts = (
        ledger.groupby("source_filename", dropna=False)
        .agg(
            ledger_row_count=("canonical_event_key", "size"),
            accepted_event_policy_candidate_count=("accepted_event_policy_candidate", lambda x: int(x.map(bool_value).sum())),
            primary_ledger_eligible_count=("primary_ledger_eligibility", lambda x: int(x.map(bool_value).sum())),
            strict_clean_only_count=("strict_clean_only_sensitivity_eligibility", lambda x: int(x.map(bool_value).sum())),
            future_epoch_ready_count=("future_epoch_construction_readiness", lambda x: int(x.map(bool_value).sum())),
        )
        .reset_index()
    )
    out = source_inventory.merge(
        crosswalk[
            [
                "source_filename",
                "task_database_id",
                "task_filename",
                "figure_participant_id",
                "mapping_confidence",
                "pipeline_mapping_allowed",
            ]
        ],
        on="source_filename",
        how="left",
        validate="one_to_one",
    ).merge(
        continuous,
        on="source_filename",
        how="left",
        validate="one_to_one",
    ).merge(
        counts,
        on="source_filename",
        how="left",
        validate="one_to_one",
    )
    count_columns = [
        "ledger_row_count",
        "accepted_event_policy_candidate_count",
        "primary_ledger_eligible_count",
        "strict_clean_only_count",
        "future_epoch_ready_count",
    ]
    out[count_columns] = out[count_columns].fillna(0).astype(int)
    out["behavioural_join_available"] = out["behavioural_participant_id"].notna() & out[
        "pipeline_mapping_allowed"
    ].map(bool_value)
    out["behavioural_join_reason_code"] = np.where(
        out["behavioural_join_available"],
        "explicit_reviewed_mapping_available",
        np.where(
            out["eeg_recording_id"].eq(11),
            "eeg11_explicitly_unmapped_from_behavioural_join",
            "behavioural_mapping_unavailable",
        ),
    )
    out["behavioural_join_reason"] = np.where(
        out["behavioural_join_available"],
        "Tracked identity contract and reviewed crosswalk permit the behavioural join.",
        np.where(
            out["eeg_recording_id"].eq(11),
            "EEG 11 is explicitly unavailable for behavioural joining and cannot use numeric equality.",
            "No permitted behavioural mapping is available for this source file.",
        ),
    )
    eeg11 = out[out["eeg_recording_id"].eq(11)]
    if (
        len(out) != 95
        or out["source_filename"].duplicated().any()
        or len(eeg11) != 1
        or bool_value(eeg11.iloc[0]["behavioural_join_available"])
        or int(eeg11.iloc[0]["ledger_row_count"]) != 0
    ):
        raise RuntimeError("source-status surface does not preserve the 95-file/EEG-11 contract")
    return out.sort_values(["eeg_recording_id", "source_filename"], kind="mergesort").reset_index(drop=True)


def build_reason_counts(ledger: pd.DataFrame) -> pd.DataFrame:
    """Count every eligibility status/reason family.

    Args:
        ledger: Canonical eligibility ledger.

    Returns:
        Long-form status/reason count table.

    Side effects:
        None.
    """

    dimensions = [
        "accepted_event_policy_candidate",
        "primary_ledger_eligibility",
        "strict_clean_only_sensitivity_eligibility",
        "continuous_derivative_availability",
        "post_ica_derivative_availability",
        "future_epoch_construction_readiness",
    ]
    rows: list[dict[str, Any]] = []
    for dimension in dimensions:
        grouped = ledger.groupby(
            [
                dimension,
                f"{dimension}_status",
                f"{dimension}_reason_code",
                f"{dimension}_reason",
            ],
            dropna=False,
            sort=True,
        ).size()
        for key, count in grouped.items():
            rows.append(
                {
                    "eligibility_dimension": dimension,
                    "value": bool_value(key[0]),
                    "status": key[1],
                    "reason_code": key[2],
                    "reason": key[3],
                    "row_count": int(count),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["eligibility_dimension", "value", "reason_code"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def build_event_policy_reason_counts(ledger: pd.DataFrame) -> pd.DataFrame:
    """Count accepted-policy outcomes with primary/strict subdivisions.

    Args:
        ledger: Canonical eligibility ledger.

    Returns:
        One row per event-policy reason code and readable reason.

    Side effects:
        None.
    """

    return (
        ledger.groupby(["event_policy_reason_code", "event_policy_reason"], sort=True, dropna=False)
        .agg(
            row_count=("canonical_event_key", "size"),
            primary_row_count=("primary_ledger_eligibility", lambda x: int(x.map(bool_value).sum())),
            strict_clean_row_count=("strict_clean_only_sensitivity_eligibility", lambda x: int(x.map(bool_value).sum())),
        )
        .reset_index()
        .sort_values("event_policy_reason_code", kind="mergesort")
        .reset_index(drop=True)
    )


def build_continuous_status_counts(ledger: pd.DataFrame) -> pd.DataFrame:
    """Summarize row counts across continuous terminal/derivative routes.

    Args:
        ledger: Canonical eligibility ledger.

    Returns:
        Continuous terminal/completion/derivative/ICA route breakdown with
        primary, strict, and future-ready counts.

    Side effects:
        None.
    """

    group_columns = [
        "continuous_v2_terminal_status",
        "continuous_v2_completion_class",
        "continuous_derivative_kind",
        "ica_terminal_route",
    ]
    filled = ledger.copy()
    filled[group_columns] = filled[group_columns].fillna("unavailable")
    return (
        filled.groupby(group_columns, sort=True, dropna=False)
        .agg(
            ledger_row_count=("canonical_event_key", "size"),
            primary_row_count=("primary_ledger_eligibility", lambda x: int(x.map(bool_value).sum())),
            strict_clean_row_count=("strict_clean_only_sensitivity_eligibility", lambda x: int(x.map(bool_value).sum())),
            future_epoch_ready_count=("future_epoch_construction_readiness", lambda x: int(x.map(bool_value).sum())),
        )
        .reset_index()
    )


def build_count_summary(ledger: pd.DataFrame) -> pd.DataFrame:
    """Build the compact canonical surface count table.

    Args:
        ledger: Canonical eligibility ledger.

    Returns:
        Named surface/count rows.

    Side effects:
        None.
    """

    primary = ledger["primary_ledger_eligibility"].map(bool_value)
    strict = ledger["strict_clean_only_sensitivity_eligibility"].map(bool_value)
    return pd.DataFrame(
        [
            {"surface": "canonical_ledger_rows", "row_count": len(ledger)},
            {"surface": "accepted_event_policy_candidates", "row_count": int(ledger["accepted_event_policy_candidate"].map(bool_value).sum())},
            {"surface": "primary_ledger_eligibility", "row_count": int(primary.sum())},
            {"surface": "strict_clean_only_sensitivity_eligibility", "row_count": int(strict.sum())},
            {"surface": "primary_warning_only_difference", "row_count": int((primary & ~strict).sum())},
            {"surface": "continuous_derivative_available", "row_count": int(ledger["continuous_derivative_availability"].map(bool_value).sum())},
            {"surface": "post_ica_derivative_available", "row_count": int(ledger["post_ica_derivative_availability"].map(bool_value).sum())},
            {"surface": "future_epoch_construction_ready", "row_count": int(ledger["future_epoch_construction_readiness"].map(bool_value).sum())},
            {"surface": "possible_event_sequence_discrepancies_unavailable", "row_count": int(ledger["event_policy_reason_code"].eq("possible_event_sequence_discrepancy_unavailable").sum())},
        ]
    )


def build_summary_payload(
    ledger: pd.DataFrame,
    source_status: pd.DataFrame,
    counts: Mapping[str, int],
    continuous_run: Mapping[str, Any],
    input_hashes: Mapping[str, Any],
) -> dict[str, Any]:
    """Build deterministic machine-readable ledger summary content.

    Args:
        ledger: Canonical eligibility ledger.
        source_status: File-level source surface.
        counts: Validated accepted counts.
        continuous_run: Selected continuous-v2 run manifest.
        input_hashes: Stable input path/hash records.

    Returns:
        JSON-safe deterministic summary dictionary.

    Side effects:
        None.
    """

    primary = ledger["primary_ledger_eligibility"].map(bool_value)
    strict = ledger["strict_clean_only_sensitivity_eligibility"].map(bool_value)
    continuous_breakdown = (
        ledger.assign(
            continuous_v2_terminal_status=ledger["continuous_v2_terminal_status"].fillna("unavailable"),
            continuous_derivative_kind=ledger["continuous_derivative_kind"].fillna("none"),
        )
        .groupby(["continuous_v2_terminal_status", "continuous_derivative_kind"], sort=True)
        .agg(
            ledger_rows=("canonical_event_key", "size"),
            primary_rows=("primary_ledger_eligibility", lambda x: int(x.map(bool_value).sum())),
            strict_clean_rows=("strict_clean_only_sensitivity_eligibility", lambda x: int(x.map(bool_value).sum())),
            future_ready_rows=("future_epoch_construction_readiness", lambda x: int(x.map(bool_value).sum())),
        )
        .reset_index()
        .to_dict("records")
    )
    exceptional: dict[str, Any] = {}
    for filename in ("demi_49 Data.edf", "demi_54_1 Data.edf", "demi_86 Data.edf"):
        rows = ledger[ledger["source_filename"].eq(filename)]
        exceptional[filename] = {
            "ledger_rows": len(rows),
            "primary_rows": int(rows["primary_ledger_eligibility"].map(bool_value).sum()),
            "strict_clean_rows": int(rows["strict_clean_only_sensitivity_eligibility"].map(bool_value).sum()),
            "future_ready_rows": int(rows["future_epoch_construction_readiness"].map(bool_value).sum()),
            "continuous_terminal_status": str(rows["continuous_v2_terminal_status"].dropna().iloc[0]) if not rows.empty else "",
            "continuous_completion_class": str(rows["continuous_v2_completion_class"].dropna().iloc[0]) if not rows.empty else "",
            "continuous_derivative_kind": str(rows["continuous_derivative_kind"].dropna().iloc[0]) if not rows.empty else "",
            "continuous_qc_warning_codes": str(rows["continuous_qc_warning_codes"].dropna().iloc[0]) if not rows.empty else "",
            "interpolated_channel_count": nullable_int(rows["interpolated_channel_count"].dropna().iloc[0]) if not rows.empty else None,
            "interpolated_channel_proportion": float(rows["interpolated_channel_proportion"].dropna().iloc[0]) if not rows.empty else None,
            "ica_terminal_route": str(rows["ica_terminal_route"].dropna().iloc[0]) if not rows.empty else "",
        }
    return {
        "schema_version": 1,
        "ledger_version": LEDGER_VERSION,
        "accepted_event_policy_version": POLICY_VERSION,
        "accepted_event_policy_date": "2026-07-11",
        "continuous_v2_run_id": str(continuous_run["run_id"]),
        "canonical_ledger_row_count": len(ledger),
        "source_status_row_count": len(source_status),
        "primary_candidate_row_count": int(primary.sum()),
        "strict_clean_only_row_count": int(strict.sum()),
        "warning_only_difference_row_count": int((primary & ~strict).sum()),
        "possible_event_sequence_discrepancy_unavailable_row_count": int(
            ledger["event_policy_reason_code"].eq("possible_event_sequence_discrepancy_unavailable").sum()
        ),
        "continuous_derivative_available_row_count": int(ledger["continuous_derivative_availability"].map(bool_value).sum()),
        "post_ica_derivative_available_row_count": int(ledger["post_ica_derivative_availability"].map(bool_value).sum()),
        "future_epoch_construction_ready_row_count": int(ledger["future_epoch_construction_readiness"].map(bool_value).sum()),
        "continuous_status_breakdown": continuous_breakdown,
        "exceptional_continuous_routes": exceptional,
        "input_authorities": input_hashes,
        "canonical_ledger_content_sha256": stable_frame_hash(ledger),
        "epochs_written": False,
        "signal_derivatives_written": False,
        "participant_inclusion_decided": False,
        "analytic_inclusion_decided": False,
    }


def build_markdown_summary(summary: Mapping[str, Any]) -> str:
    """Render a concise deterministic Markdown review summary.

    Args:
        summary: Machine-readable summary payload.

    Returns:
        Markdown text.

    Side effects:
        None.
    """

    exceptional = summary["exceptional_continuous_routes"]
    lines = [
        "# Event/epoch eligibility ledger v1",
        "",
        "The accepted event-policy ledger is complete. It records event-policy, continuous-derivative, post-ICA, and future epoch-readiness statuses separately. It does not construct epochs or decide participant/analytic inclusion.",
        "",
        "## Accepted surfaces",
        "",
        f"- Canonical ledger rows: {summary['canonical_ledger_row_count']:,}",
        f"- Primary accepted event-policy candidates: {summary['primary_candidate_row_count']:,}",
        f"- Strict-clean-only sensitivity candidates: {summary['strict_clean_only_row_count']:,}",
        f"- Primary-only duration-warning difference: {summary['warning_only_difference_row_count']:,}",
        f"- Reviewed possible sequence discrepancies kept unavailable: {summary['possible_event_sequence_discrepancy_unavailable_row_count']:,}",
        f"- Future epoch-construction ready under the ordinary post-ICA route: {summary['future_epoch_construction_ready_row_count']:,}",
        "",
        "## Continuous status",
        "",
        f"- Continuous-v2 run: `{summary['continuous_v2_run_id']}`",
        f"- Rows with a canonical continuous derivative: {summary['continuous_derivative_available_row_count']:,}",
        f"- Rows with a post-ICA derivative: {summary['post_ica_derivative_available_row_count']:,}",
    ]
    for filename in ("demi_49 Data.edf", "demi_54_1 Data.edf"):
        item = exceptional[filename]
        lines.append(
            f"- `{filename}` remains `{item['continuous_completion_class']}` with `{item['continuous_qc_warning_codes']}` ({item['interpolated_channel_count']}/30 interpolated); the warning does not remove accepted event rows."
        )
    id86 = exceptional["demi_86 Data.edf"]
    lines.append(
        f"- `demi_86 Data.edf` remains `{id86['continuous_terminal_status']}` on `{id86['ica_terminal_route']}` with a retained pre-ICA derivative and no post-ICA derivative. Its {id86['primary_rows']:,} accepted event candidates are not marked future-ready; later analytic handling remains undecided."
    )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "No EDF event was mutated. No split file was merged or renumbered. No epoch, AutoReject, CSD, time-frequency, feature, repaired EDF, participant-inclusion, or analytic-inclusion output was created. Actual epoch construction remains a separate authorization boundary.",
            "",
        ]
    )
    return "\n".join(lines)


def atomic_write_text(path: Path, text: str) -> None:
    """Atomically replace one UTF-8 text file in its final directory.

    Args:
        path: Final output path.
        text: Complete text payload.

    Side effects:
        Writes a temporary sibling and replaces ``path`` atomically.
    """

    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Atomically replace one deterministic CSV table.

    Args:
        path: Final CSV path.
        frame: Data frame in final row/column order.

    Side effects:
        Writes a temporary sibling and replaces ``path`` atomically.
    """

    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n", na_rep="")
    os.replace(temporary, path)


def atomic_write_parquet(path: Path, frame: pd.DataFrame) -> None:
    """Atomically replace the canonical Parquet ledger.

    Args:
        path: Final Parquet path.
        frame: Canonical ledger.

    Side effects:
        Writes a temporary sibling and replaces ``path`` atomically.
    """

    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    os.replace(temporary, path)


def write_outputs(result: Mapping[str, Any], output_dir: Path, repo_root: Path) -> dict[str, Any]:
    """Write every versioned local ledger artifact atomically.

    Args:
        result: In-memory output dictionary returned by ``build_ledger``.
        output_dir: Final local versioned namespace.
        repo_root: Repository root for stable path rendering.

    Returns:
        Deterministic run-manifest payload including output table hashes.

    Side effects:
        Creates the ignored output directory and writes Parquet, CSV, JSON,
        and Markdown files only.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    ledger = result["ledger"]
    tables = {
        LEDGER_CSV_FILENAME: ledger,
        SOURCE_STATUS_FILENAME: result["source_status"],
        COUNT_SUMMARY_FILENAME: result["count_summary"],
        REASON_COUNTS_FILENAME: result["reason_counts"],
        EVENT_POLICY_REASON_COUNTS_FILENAME: result["event_policy_reason_counts"],
        CONTINUOUS_STATUS_COUNTS_FILENAME: result["continuous_status_counts"],
    }
    atomic_write_parquet(output_dir / LEDGER_PARQUET_FILENAME, ledger)
    for filename, frame in tables.items():
        atomic_write_csv(output_dir / filename, frame)
    atomic_write_text(
        output_dir / SUMMARY_JSON_FILENAME,
        json.dumps(result["summary"], indent=2, sort_keys=True) + "\n",
    )
    atomic_write_text(
        output_dir / SUMMARY_MD_FILENAME, build_markdown_summary(result["summary"])
    )

    output_records = {}
    for path in sorted(output_dir.iterdir()):
        if not path.is_file() or path.name == RUN_MANIFEST_FILENAME or path.name.startswith("."):
            continue
        output_records[path.name] = {
            "path": path.resolve().relative_to(repo_root.resolve()).as_posix(),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
    manifest = {
        "schema_version": 1,
        "ledger_version": LEDGER_VERSION,
        "accepted_event_policy_version": POLICY_VERSION,
        "continuous_v2_run_id": result["summary"]["continuous_v2_run_id"],
        "input_authorities": result["input_hashes"],
        "outputs": output_records,
        "canonical_ledger_content_sha256": stable_frame_hash(ledger),
        "primary_candidate_row_count": EXPECTED_PRIMARY_COUNT,
        "strict_clean_only_row_count": EXPECTED_STRICT_CLEAN_COUNT,
        "warning_only_difference_row_count": EXPECTED_WARNING_ONLY_COUNT,
        "possible_event_sequence_discrepancy_unavailable_row_count": EXPECTED_POSSIBLE_DISCREPANCY_COUNT,
        "epochs_written": False,
        "signal_derivatives_written": False,
    }
    atomic_write_text(
        output_dir / RUN_MANIFEST_FILENAME,
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    assert_no_signal_outputs(output_dir)
    return manifest


def parse_args() -> argparse.Namespace:
    """Parse optional exact continuous-run selection.

    Returns:
        Parsed command-line namespace.

    Side effects:
        Reads process arguments.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--continuous-run-id",
        help="Exact finalized continuous_v2 run directory name; default selects the latest qualifying run.",
    )
    return parser.parse_args()


def main() -> int:
    """Build, validate, and atomically write the accepted eligibility ledger.

    Returns:
        Process exit status zero on success.

    Side effects:
        Reads accepted local/tracked evidence and writes only the versioned
        ignored ledger namespace.
    """

    args = parse_args()
    repo_root = repo_root_from_script()
    result = build_ledger(repo_root, args.continuous_run_id)
    output_dir = repo_root / OUTPUT_DIR
    write_outputs(result, output_dir, repo_root)
    summary = result["summary"]
    print(f"Wrote canonical ledger: {output_dir / LEDGER_PARQUET_FILENAME}")
    print(f"Canonical rows: {summary['canonical_ledger_row_count']:,}")
    print(f"Primary candidates: {summary['primary_candidate_row_count']:,}")
    print(f"Strict-clean-only candidates: {summary['strict_clean_only_row_count']:,}")
    print(f"Warning-only difference: {summary['warning_only_difference_row_count']:,}")
    print("Epochs written: false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
