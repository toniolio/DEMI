"""Inventory DEMI EDF annotations and physical Trigger events side by side.

This script is the read-only event-source discovery layer for the restarted
DEMI EEG analysis. Historical conversion code used EDF annotations when they
were sufficiently populated and otherwise called MNE find_events on the
physical Trigger channel. The first new MNE audits inspected annotations only,
which made an annotation-export failure look like an event-acquisition
failure. This script restores the historically grounded dual-source view.

Inputs:

- raw EDF files below _Data/eeg/raw;
- analysis/eeg_mne/eeg_behavior_identity_contract.csv;
- MNE and pandas from the EEG Python environment.

Local outputs below _Data/eeg/event_source_inventory_v1:

- dual_source_events.csv: every annotation and physical event with raw values,
  raw codes, comparison-only candidate codes, samples, onsets, extraction
  parameters, file identity, and mapping provenance;
- selected_event_evidence.csv: copied rows from the conservative selected
  source, with selection reason and normalization status;
- event_source_file_inventory.csv: one row per EDF with source assessments,
  discrepancy class, selected source, and unresolved issue;
- event_source_inventory_summary.md: compact aggregate results;
- event_source_run_manifest.json: code/environment/extraction provenance.

The physical extraction uses the same core MNE settings as Austin's historical
edf2bids.py: shortest_event=1, mask=65280, and mask_type=not_and. The Trigger
channel is named explicitly, and all effective parameters are recorded.

Safety boundaries:

- EDF files are opened read-only.
- The script never sets or rewrites Raw annotations.
- It never modifies Trigger samples.
- It performs no filtering, notch filtering, referencing, interpolation, ICA,
  CSD, event repair, epoch construction, or participant/trial decision.
- Candidate source-code normalization is recorded but never silently applied.
- Output is versioned and local-only; annotation-only audit directories are
  not overwritten.

Run from the repository root:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/05_inventory_event_sources.py
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import mne
import numpy as np
import pandas as pd

from event_source_contract import (
    EVENT_CODE_TO_NAME,
    IDENTITY_CONTRACT_VERSION,
    TASK_EVENT_CODES,
    IdentityMapping,
    StreamAssessment,
    StreamComparison,
    SourceSelection,
    compare_event_streams,
    event_name_for_code,
    integer_or_none,
    load_identity_contract,
    physical_code_candidate,
    require_identity_mapping,
    selected_event_rows,
    select_event_source,
    validate_identity_coverage,
)


RAW_EEG_DIR = Path("_Data") / "eeg" / "raw"
OUTPUT_DIR = Path("_Data") / "eeg" / "event_source_inventory_v1"
IDENTITY_CONTRACT_PATH = (
    Path("analysis") / "eeg_mne" / "eeg_behavior_identity_contract.csv"
)

DUAL_SOURCE_EVENTS_FILENAME = "dual_source_events.csv"
SELECTED_EVENTS_FILENAME = "selected_event_evidence.csv"
FILE_INVENTORY_FILENAME = "event_source_file_inventory.csv"
SUMMARY_FILENAME = "event_source_inventory_summary.md"
RUN_MANIFEST_FILENAME = "event_source_run_manifest.json"

TRIGGER_CHANNEL_NAME = "Trigger"
PHYSICAL_EXTRACTION_PARAMETERS = {
    "function": "mne.find_events",
    "stim_channel": TRIGGER_CHANNEL_NAME,
    "output": "onset",
    "consecutive": "increasing",
    "min_duration": 0.0,
    "shortest_event": 1,
    "mask": 65280,
    "uint_cast": False,
    "mask_type": "not_and",
    "initial_event": False,
}
ANNOTATION_EXTRACTION_PARAMETERS = {
    "function": "raw.annotations iteration",
    "sample_conversion": "raw.time_as_index(use_rounding=True)+raw.first_samp",
}

EDF_FILENAME_RE = re.compile(
    r"^demi_(?P<eeg_recording_id>\d{1,3})"
    r"(?:_(?P<split_part>\d+))?"
    r"(?P<label>.*)\.edf$",
    flags=re.IGNORECASE,
)

EVENT_OUTPUT_COLUMNS = [
    "eeg_recording_id",
    "eeg_recording_id_padded",
    "behavioural_participant_id",
    "identity_mapping_status",
    "identity_contract_version",
    "source_filename",
    "file_path",
    "file_role",
    "split_part",
    "session_date",
    "sampling_frequency_hz",
    "source_type",
    "event_row_order",
    "sample",
    "onset_seconds",
    "duration_seconds",
    "raw_value",
    "previous_raw_value",
    "raw_code",
    "previous_code",
    "normalized_candidate_code",
    "normalization_candidate_rule",
    "normalization_candidate_applied",
    "candidate_event_name",
    "task_event_candidate",
    "source_extraction_parameters_json",
    "selected_for_event_evidence",
    "selected_source",
    "selection_status",
    "selection_reason",
    "selection_unresolved_issue",
    "normalization_applied",
    "selected_normalized_code",
    "selected_event_name",
]


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


def relative_to_repo(path: Path, repo_root: Path) -> str:
    """Return a stable repository-relative POSIX path.

    Args:
        path: Path to describe.
        repo_root: Repository root.

    Returns:
        Relative POSIX path when path is inside the repository, otherwise the
        absolute POSIX path.

    Side effects:
        None.
    """

    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def discover_edf_files(raw_dir: Path) -> list[Path]:
    """Return sorted raw EDF files directly below raw_dir.

    Args:
        raw_dir: Expected raw EDF directory.

    Returns:
        Sorted EDF paths.

    Side effects:
        Reads the directory. Raises RuntimeError when it is missing or empty.
    """

    if not raw_dir.exists():
        raise RuntimeError(f"raw EEG directory is missing: {raw_dir}")
    files = sorted(
        path
        for path in raw_dir.iterdir()
        if path.is_file() and path.suffix.lower() == ".edf"
    )
    if not files:
        raise RuntimeError(f"no raw EDF files found in {raw_dir}")
    return files


def parse_edf_filename(path: Path) -> dict[str, Any]:
    """Parse EEG ID, split part, and file role from a DEMI EDF filename.

    Args:
        path: EDF path; only the filename is parsed.

    Returns:
        Parsed identity dictionary.

    Side effects:
        None. Raises RuntimeError when the filename is not recognized.
    """

    match = EDF_FILENAME_RE.match(path.name)
    if match is None:
        raise RuntimeError(
            f"raw EDF filename does not match the DEMI contract: {path.name}"
        )

    eeg_id = int(match.group("eeg_recording_id"))
    split_text = match.group("split_part")
    split_part = int(split_text) if split_text else None
    label = match.group("label").lower()
    if split_part is not None:
        file_role = "split_part"
    elif "concatenated" in label:
        file_role = "concatenated"
    else:
        file_role = "single"
    return {
        "eeg_recording_id": eeg_id,
        "eeg_recording_id_padded": f"{eeg_id:03d}",
        "file_role": file_role,
        "split_part": split_part,
    }


def numeric_annotation_code(description: Any) -> int | None:
    """Return a numeric annotation description as an integer.

    Args:
        description: MNE annotation description.

    Returns:
        Integer code or None for labels such as impedance checks and file
        start.

    Side effects:
        None.
    """

    text = str(description).strip()
    if not text:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if not value.is_integer():
        return None
    return int(value)


def annotation_event_name(description: str, raw_code: int | None) -> str:
    """Return the semantic name for one annotation row.

    Args:
        description: Raw annotation description.
        raw_code: Parsed numeric code, if any.

    Returns:
        DEMI task event name, file-start marker name, or empty string.

    Side effects:
        None.
    """

    if description.strip().lower() == "file start":
        return "file start"
    return event_name_for_code(raw_code)


def session_date_from_raw(raw: mne.io.BaseRaw) -> str:
    """Return the acquisition date from MNE metadata when available.

    Args:
        raw: Open MNE Raw object.

    Returns:
        ISO date text or an empty string.

    Side effects:
        None.
    """

    meas_date = raw.info.get("meas_date")
    if meas_date is None:
        return ""
    return meas_date.date().isoformat()


def event_base(
    *,
    parsed: dict[str, Any],
    mapping: IdentityMapping,
    edf_path: Path,
    repo_root: Path,
    session_date: str,
    sampling_frequency_hz: float,
) -> dict[str, Any]:
    """Build file and identity provenance shared by every event row.

    Args:
        parsed: Parsed EDF filename fields.
        mapping: Explicit identity-contract row.
        edf_path: Source EDF path.
        repo_root: Repository root.
        session_date: Acquisition date text.
        sampling_frequency_hz: Raw sampling frequency.

    Returns:
        Shared event provenance dictionary.

    Side effects:
        None.
    """

    return {
        **parsed,
        "behavioural_participant_id": mapping.behavioural_participant_id,
        "identity_mapping_status": mapping.mapping_status,
        "identity_contract_version": mapping.contract_version,
        "source_filename": edf_path.name,
        "file_path": relative_to_repo(edf_path, repo_root),
        "session_date": session_date,
        "sampling_frequency_hz": sampling_frequency_hz,
    }


def extract_annotation_events(
    raw: mne.io.BaseRaw,
    base: dict[str, Any],
) -> list[dict[str, Any]]:
    """Extract annotation rows without changing the Raw object.

    Args:
        raw: Open MNE Raw object.
        base: Shared file and identity provenance.

    Returns:
        One dictionary per EDF annotation.

    Side effects:
        Reads Raw annotation metadata only.
    """

    extraction_json = json.dumps(
        ANNOTATION_EXTRACTION_PARAMETERS, sort_keys=True
    )
    rows: list[dict[str, Any]] = []
    for event_row_order, annotation in enumerate(raw.annotations, start=1):
        description = str(annotation["description"])
        raw_code = numeric_annotation_code(description)
        candidate_code = raw_code
        candidate_name = annotation_event_name(description, raw_code)
        onset = float(annotation["onset"])
        sample = int(
            raw.time_as_index([onset], use_rounding=True)[0] + raw.first_samp
        )
        rows.append(
            {
                **base,
                "source_type": "annotation",
                "event_row_order": event_row_order,
                "sample": sample,
                "onset_seconds": onset,
                "duration_seconds": float(annotation["duration"]),
                "raw_value": description,
                "previous_raw_value": "",
                "raw_code": raw_code,
                "previous_code": None,
                "normalized_candidate_code": candidate_code,
                "normalization_candidate_rule": "annotation_label_preserved",
                "normalization_candidate_applied": False,
                "candidate_event_name": candidate_name,
                "task_event_candidate": raw_code in TASK_EVENT_CODES,
                "source_extraction_parameters_json": extraction_json,
                "selected_for_event_evidence": False,
                "selected_source": "",
                "selection_status": "",
                "selection_reason": "",
                "selection_unresolved_issue": "",
                "normalization_applied": False,
                "selected_normalized_code": None,
                "selected_event_name": "",
            }
        )
    return rows


def extract_physical_trigger_events(
    raw: mne.io.BaseRaw,
    base: dict[str, Any],
) -> list[dict[str, Any]]:
    """Extract physical Trigger transitions using historical MNE settings.

    Args:
        raw: Open MNE Raw object.
        base: Shared file and identity provenance.

    Returns:
        One dictionary per MNE physical event. Stored Trigger values, masked
        codes, comparison-only candidates, samples, and onsets are retained.

    Side effects:
        Reads the Trigger channel into memory. Does not preload or modify other
        channels and does not alter the Raw object.
    """

    if TRIGGER_CHANNEL_NAME not in raw.ch_names:
        return []

    events = mne.find_events(
        raw,
        stim_channel=TRIGGER_CHANNEL_NAME,
        output="onset",
        consecutive="increasing",
        min_duration=0.0,
        shortest_event=1,
        mask=65280,
        uint_cast=False,
        mask_type="not_and",
        initial_event=False,
        verbose="ERROR",
    )
    trigger = raw.get_data(
        picks=[TRIGGER_CHANNEL_NAME],
        reject_by_annotation=None,
        verbose="ERROR",
    )[0]
    extraction_json = json.dumps(
        PHYSICAL_EXTRACTION_PARAMETERS, sort_keys=True
    )
    sampling_frequency_hz = float(raw.info["sfreq"])

    rows: list[dict[str, Any]] = []
    for event_row_order, event in enumerate(events, start=1):
        sample = int(event[0])
        data_index = sample - int(raw.first_samp)
        raw_value = int(round(float(trigger[data_index])))
        previous_raw_value = (
            int(round(float(trigger[data_index - 1])))
            if data_index > 0
            else None
        )
        raw_code = int(event[2])
        previous_code = int(event[1])
        candidate_code, candidate_rule = physical_code_candidate(raw_code)
        rows.append(
            {
                **base,
                "source_type": "physical_trigger",
                "event_row_order": event_row_order,
                "sample": sample,
                "onset_seconds": (
                    sample - int(raw.first_samp)
                )
                / sampling_frequency_hz,
                "duration_seconds": 0.0,
                "raw_value": raw_value,
                "previous_raw_value": previous_raw_value,
                "raw_code": raw_code,
                "previous_code": previous_code,
                "normalized_candidate_code": candidate_code,
                "normalization_candidate_rule": candidate_rule,
                "normalization_candidate_applied": False,
                "candidate_event_name": event_name_for_code(candidate_code),
                "task_event_candidate": candidate_code in TASK_EVENT_CODES,
                "source_extraction_parameters_json": extraction_json,
                "selected_for_event_evidence": False,
                "selected_source": "",
                "selection_status": "",
                "selection_reason": "",
                "selection_unresolved_issue": "",
                "normalization_applied": False,
                "selected_normalized_code": None,
                "selected_event_name": "",
            }
        )
    return rows


def assessment_fields(
    prefix: str,
    assessment: StreamAssessment,
) -> dict[str, Any]:
    """Flatten one stream assessment into stable output columns.

    Args:
        prefix: Column prefix such as annotation or physical_raw.
        assessment: Stream assessment to flatten.

    Returns:
        Dictionary of scalar/JSON output values.

    Side effects:
        None.
    """

    return {
        f"{prefix}_stream_status": assessment.status,
        f"{prefix}_stream_coherent": assessment.coherent,
        f"{prefix}_task_event_count": assessment.task_event_count,
        f"{prefix}_core_code_counts_json": json.dumps(
            assessment.core_code_counts, sort_keys=True
        ),
        f"{prefix}_optional_code_counts_json": json.dumps(
            assessment.optional_code_counts, sort_keys=True
        ),
        f"{prefix}_minimum_core_count": assessment.minimum_core_count,
        f"{prefix}_maximum_core_count": assessment.maximum_core_count,
    }


def counts_json(rows: list[dict[str, Any]], field: str) -> str:
    """Return sorted JSON counts for one event-row field.

    Args:
        rows: Event dictionaries.
        field: Field to count.

    Returns:
        JSON object text.

    Side effects:
        None.
    """

    values = [
        str(row[field])
        for row in rows
        if row.get(field) is not None and str(row.get(field)).strip()
    ]
    return json.dumps(dict(sorted(Counter(values).items())), sort_keys=True)


def apply_selection_flags(
    all_rows: list[dict[str, Any]],
    selected_rows: list[dict[str, Any]],
) -> None:
    """Mark selected rows in the combined dual-source table.

    Args:
        all_rows: Mutable annotation plus physical event rows.
        selected_rows: Copied rows returned by selected_event_rows.

    Returns:
        None.

    Side effects:
        Mutates only in-memory dictionaries in all_rows. Raw MNE objects and
        EDF files are untouched.
    """

    selected_keys = {
        (
            row["source_type"],
            int(row["event_row_order"]),
            int(row["sample"]),
        )
        for row in selected_rows
    }
    if not selected_rows:
        return

    provenance = selected_rows[0]
    for row in all_rows:
        key = (
            row["source_type"],
            int(row["event_row_order"]),
            int(row["sample"]),
        )
        row["selected_for_event_evidence"] = key in selected_keys
        row["selected_source"] = provenance["selected_source"]
        row["selection_status"] = provenance["selection_status"]
        row["selection_reason"] = provenance["selection_reason"]
        row["selection_unresolved_issue"] = provenance[
            "selection_unresolved_issue"
        ]
        if key in selected_keys:
            row["normalization_applied"] = provenance[
                "normalization_applied"
            ]
            row["selected_normalized_code"] = row["raw_code"]
            row["selected_event_name"] = (
                annotation_event_name(str(row["raw_value"]), row["raw_code"])
                if row["source_type"] == "annotation"
                else event_name_for_code(row["raw_code"])
            )


def file_inventory_row(
    *,
    base: dict[str, Any],
    edf_path: Path,
    raw: mne.io.BaseRaw,
    annotation_rows: list[dict[str, Any]],
    physical_rows: list[dict[str, Any]],
    comparison: StreamComparison,
    selection: SourceSelection,
) -> dict[str, Any]:
    """Build one file-level source comparison row.

    Args:
        base: Shared file/identity provenance.
        edf_path: Source EDF path.
        raw: Open MNE Raw object.
        annotation_rows: Annotation-source event rows.
        physical_rows: Physical-source event rows.
        comparison: Stream comparison.
        selection: Conservative source selection.

    Returns:
        File-level inventory dictionary.

    Side effects:
        Reads filesystem stat metadata for edf_path.
    """

    stat = edf_path.stat()
    candidate_difference_count = sum(
        integer_or_none(row.get("raw_code"))
        != integer_or_none(row.get("normalized_candidate_code"))
        for row in physical_rows
    )
    return {
        **base,
        "read_status": "ok",
        "read_error_type": "",
        "read_error_message": "",
        "file_size_bytes": stat.st_size,
        "file_mtime_ns": stat.st_mtime_ns,
        "recording_duration_seconds": float(
            raw.n_times / float(raw.info["sfreq"])
        ),
        "n_channels": len(raw.ch_names),
        "trigger_channel_present": TRIGGER_CHANNEL_NAME in raw.ch_names,
        "annotation_event_count": len(annotation_rows),
        "annotation_raw_code_counts_json": counts_json(
            annotation_rows, "raw_code"
        ),
        "physical_event_count": len(physical_rows),
        "physical_raw_code_counts_json": counts_json(
            physical_rows, "raw_code"
        ),
        "physical_candidate_code_counts_json": counts_json(
            physical_rows, "normalized_candidate_code"
        ),
        "physical_candidate_difference_count": candidate_difference_count,
        **assessment_fields("annotation", comparison.annotation),
        **assessment_fields("physical_raw", comparison.physical_raw),
        **assessment_fields(
            "physical_candidate", comparison.physical_candidate
        ),
        "source_discrepancy_class": comparison.discrepancy_class,
        "compared_event_count": comparison.compared_event_count,
        "same_sample_and_candidate_code_count": (
            comparison.same_sample_and_candidate_code_count
        ),
        "raw_code_difference_count": comparison.raw_code_difference_count,
        "maximum_absolute_sample_difference": (
            comparison.maximum_absolute_sample_difference
        ),
        "ordered_candidate_streams_match": (
            comparison.ordered_candidate_streams_match
        ),
        "selected_source": selection.selected_source or "",
        "selection_status": selection.selection_status,
        "selection_reason": selection.selection_reason,
        "selection_unresolved_issue": selection.unresolved_issue,
        "normalization_applied": selection.normalization_applied,
    }


def error_inventory_row(
    *,
    parsed: dict[str, Any],
    mapping: IdentityMapping,
    edf_path: Path,
    repo_root: Path,
    error: Exception,
) -> dict[str, Any]:
    """Build a file-level row when one EDF cannot be inspected.

    Args:
        parsed: Parsed EDF filename fields.
        mapping: Explicit identity mapping.
        edf_path: Source EDF path.
        repo_root: Repository root.
        error: Caught inspection error.

    Returns:
        Error row retaining file and identity provenance.

    Side effects:
        Reads filesystem stat metadata when available.
    """

    stat = edf_path.stat()
    return {
        **parsed,
        "behavioural_participant_id": mapping.behavioural_participant_id,
        "identity_mapping_status": mapping.mapping_status,
        "identity_contract_version": mapping.contract_version,
        "source_filename": edf_path.name,
        "file_path": relative_to_repo(edf_path, repo_root),
        "session_date": "",
        "sampling_frequency_hz": np.nan,
        "read_status": "mne_read_error",
        "read_error_type": type(error).__name__,
        "read_error_message": str(error),
        "file_size_bytes": stat.st_size,
        "file_mtime_ns": stat.st_mtime_ns,
        "selection_status": "read_error_no_source_selected",
        "selection_reason": "EDF source inventory failed",
        "selection_unresolved_issue": "mne_read_error",
        "selected_source": "",
        "normalization_applied": False,
    }


def inspect_one_edf(
    *,
    edf_path: Path,
    repo_root: Path,
    mapping: IdentityMapping,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Inspect one EDF and return dual-source, selected, and file rows.

    Args:
        edf_path: Raw EDF path.
        repo_root: Repository root.
        mapping: Explicit EEG-to-behaviour mapping.

    Returns:
        Tuple of combined source rows, selected rows, and one file summary row.

    Side effects:
        Opens the EDF read-only, reads annotations and the Trigger channel, and
        closes the Raw object. No file is modified.
    """

    parsed = parse_edf_filename(edf_path)
    raw: mne.io.BaseRaw | None = None
    try:
        raw = mne.io.read_raw_edf(
            edf_path, preload=False, verbose="ERROR"
        )
        sampling_frequency_hz = float(raw.info["sfreq"])
        base = event_base(
            parsed=parsed,
            mapping=mapping,
            edf_path=edf_path,
            repo_root=repo_root,
            session_date=session_date_from_raw(raw),
            sampling_frequency_hz=sampling_frequency_hz,
        )
        annotation_rows = extract_annotation_events(raw, base)
        physical_rows = extract_physical_trigger_events(raw, base)
        comparison = compare_event_streams(
            annotation_rows, physical_rows, sample_tolerance=1
        )
        selection = select_event_source(comparison)
        selected_rows = selected_event_rows(
            annotation_rows, physical_rows, selection
        )
        all_rows = annotation_rows + physical_rows
        apply_selection_flags(all_rows, selected_rows)
        inventory = file_inventory_row(
            base=base,
            edf_path=edf_path,
            raw=raw,
            annotation_rows=annotation_rows,
            physical_rows=physical_rows,
            comparison=comparison,
            selection=selection,
        )
        return all_rows, selected_rows, inventory
    finally:
        if raw is not None and hasattr(raw, "close"):
            raw.close()


def dataframe_with_columns(
    rows: list[dict[str, Any]],
    preferred_columns: list[str],
) -> pd.DataFrame:
    """Build a DataFrame with stable preferred columns.

    Args:
        rows: Row dictionaries.
        preferred_columns: Columns to place first and retain when rows are
            empty.

    Returns:
        DataFrame with preferred columns followed by any extras.

    Side effects:
        None.
    """

    if not rows:
        return pd.DataFrame(columns=preferred_columns)
    frame = pd.DataFrame(rows)
    ordered = preferred_columns + [
        column for column in frame.columns if column not in preferred_columns
    ]
    return frame.reindex(columns=ordered)


def run_inventory(
    repo_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the complete dual-source inventory in memory.

    Args:
        repo_root: Repository root.

    Returns:
        Dual-source event table, selected event table, and file inventory.

    Side effects:
        Reads all raw EDFs and the tracked identity contract. Does not write
        output files.
    """

    raw_dir = repo_root / RAW_EEG_DIR
    identity_path = repo_root / IDENTITY_CONTRACT_PATH
    edf_files = discover_edf_files(raw_dir)
    parsed_files = [(path, parse_edf_filename(path)) for path in edf_files]
    mappings = load_identity_contract(identity_path)
    validate_identity_coverage(
        mappings,
        [parsed["eeg_recording_id"] for _, parsed in parsed_files],
    )

    all_event_rows: list[dict[str, Any]] = []
    selected_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    for edf_path, parsed in parsed_files:
        mapping = require_identity_mapping(
            mappings, parsed["eeg_recording_id"]
        )
        try:
            events, selected, file_row = inspect_one_edf(
                edf_path=edf_path,
                repo_root=repo_root,
                mapping=mapping,
            )
        except Exception as error:
            events = []
            selected = []
            file_row = error_inventory_row(
                parsed=parsed,
                mapping=mapping,
                edf_path=edf_path,
                repo_root=repo_root,
                error=error,
            )
        all_event_rows.extend(events)
        selected_rows.extend(selected)
        file_rows.append(file_row)

    event_frame = dataframe_with_columns(
        all_event_rows, EVENT_OUTPUT_COLUMNS
    )
    selected_frame = dataframe_with_columns(
        selected_rows, EVENT_OUTPUT_COLUMNS
    )
    file_frame = pd.DataFrame(file_rows).sort_values(
        ["eeg_recording_id", "source_filename"],
        kind="mergesort",
    )
    return (
        event_frame.reset_index(drop=True),
        selected_frame.reset_index(drop=True),
        file_frame.reset_index(drop=True),
    )


def format_counts(series: pd.Series) -> list[str]:
    """Format value counts as Markdown bullets.

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
        values: Series of IDs.

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


def build_summary(
    file_inventory: pd.DataFrame,
    dual_events: pd.DataFrame,
    selected_events: pd.DataFrame,
    started_at: datetime,
) -> str:
    """Build the local Markdown inventory summary.

    Args:
        file_inventory: One-row-per-EDF source inventory.
        dual_events: Combined annotation and physical events.
        selected_events: Events copied from selected sources.
        started_at: Run start time.

    Returns:
        Markdown text.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    fallback = file_inventory[
        file_inventory["selected_source"].eq("physical_trigger")
    ]
    unresolved = file_inventory[
        file_inventory["selection_unresolved_issue"].fillna("").ne("")
    ]
    remapped = file_inventory[
        file_inventory["identity_mapping_status"].eq("explicit_remap")
    ]
    unmapped = file_inventory[
        file_inventory["identity_mapping_status"].eq("explicitly_unmapped")
    ]

    lines = [
        "# DEMI Dual-Source Event Inventory",
        "",
        f"Generated: {finished_at.isoformat(timespec='seconds')}",
        f"Duration seconds: {(finished_at - started_at).total_seconds():.1f}",
        "",
        "## Safety boundary",
        "",
        "- Read-only EDF inspection.",
        "- No raw annotation or Trigger mutation.",
        "- No preprocessing, event repair, epoching, or inclusion decision.",
        "- Candidate physical-code normalization was not applied.",
        "",
        "## Identity contract",
        "",
        f"- Contract version: {IDENTITY_CONTRACT_VERSION}",
        f"- EDF files inspected: {len(file_inventory)}",
        f"- EEG recording IDs: {file_inventory['eeg_recording_id'].nunique()}",
        f"- Explicitly remapped EDF rows: {len(remapped)} "
        f"({format_id_list(remapped['eeg_recording_id'])})",
        f"- Explicitly unmapped EDF rows: {len(unmapped)} "
        f"({format_id_list(unmapped['eeg_recording_id'])})",
        "",
        "## Source events",
        "",
        f"- Combined source-event rows: {len(dual_events)}",
        f"- Annotation rows: {int(dual_events['source_type'].eq('annotation').sum())}",
        f"- Physical Trigger rows: {int(dual_events['source_type'].eq('physical_trigger').sum())}",
        f"- Selected event rows: {len(selected_events)}",
        "",
        "Source discrepancy classes:",
        "",
        *format_counts(file_inventory["source_discrepancy_class"]),
        "",
        "Selected sources:",
        "",
        *format_counts(file_inventory["selected_source"]),
        "",
        f"- Physical-fallback EDFs: {len(fallback)} "
        f"({format_id_list(fallback['eeg_recording_id'])})",
        f"- EDFs with unresolved source issues: {len(unresolved)} "
        f"({format_id_list(unresolved['eeg_recording_id'])})",
        "",
        "## Extraction contract",
        "",
        f"- Physical parameters: {json.dumps(PHYSICAL_EXTRACTION_PARAMETERS, sort_keys=True)}",
        f"- Annotation parameters: {json.dumps(ANNOTATION_EXTRACTION_PARAMETERS, sort_keys=True)}",
        "- Physical candidate-code differences remain separate from raw codes.",
        "- Source selection is evidence selection, not final event policy.",
        "",
    ]
    return "\n".join(lines)


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a small contract file.

    Args:
        path: File to hash.

    Returns:
        Lowercase hexadecimal digest.

    Side effects:
        Reads path.
    """

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_run_manifest(
    *,
    repo_root: Path,
    file_inventory: pd.DataFrame,
    dual_events: pd.DataFrame,
    selected_events: pd.DataFrame,
    started_at: datetime,
) -> dict[str, Any]:
    """Build run-level source and environment provenance.

    Args:
        repo_root: Repository root.
        file_inventory: One-row-per-EDF inventory.
        dual_events: Combined event rows.
        selected_events: Selected event rows.
        started_at: Run start time.

    Returns:
        JSON-serializable manifest dictionary.

    Side effects:
        Hashes the small tracked identity contract.
    """

    identity_path = repo_root / IDENTITY_CONTRACT_PATH
    return {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "started": started_at.isoformat(timespec="seconds"),
        "script": "analysis/eeg_mne/05_inventory_event_sources.py",
        "identity_contract": IDENTITY_CONTRACT_PATH.as_posix(),
        "identity_contract_version": IDENTITY_CONTRACT_VERSION,
        "identity_contract_sha256": sha256_file(identity_path),
        "raw_eeg_directory": RAW_EEG_DIR.as_posix(),
        "output_directory": OUTPUT_DIR.as_posix(),
        "edf_files": int(len(file_inventory)),
        "eeg_recording_ids": int(
            file_inventory["eeg_recording_id"].nunique()
        ),
        "dual_source_event_rows": int(len(dual_events)),
        "selected_event_rows": int(len(selected_events)),
        "mne_version": mne.__version__,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "python_version": platform.python_version(),
        "physical_extraction_parameters": PHYSICAL_EXTRACTION_PARAMETERS,
        "annotation_extraction_parameters": (
            ANNOTATION_EXTRACTION_PARAMETERS
        ),
        "raw_mutation_allowed": False,
        "event_repair_allowed": False,
        "epoching_allowed": False,
        "normalization_applied": False,
    }


def validate_output_directory(output_dir: Path, repo_root: Path) -> None:
    """Require the versioned output directory below local _Data/eeg.

    Args:
        output_dir: Proposed output directory.
        repo_root: Repository root.

    Returns:
        None.

    Side effects:
        None. Raises RuntimeError for unsafe output paths.
    """

    resolved = output_dir.resolve()
    required_parent = (repo_root / "_Data" / "eeg").resolve()
    raw_dir = (repo_root / RAW_EEG_DIR).resolve()
    try:
        resolved.relative_to(required_parent)
    except ValueError as error:
        raise RuntimeError(
            f"event inventory output must stay below {required_parent}"
        ) from error
    if resolved == raw_dir or raw_dir in resolved.parents:
        raise RuntimeError("event inventory output cannot be inside raw EDF data")


def write_outputs(
    *,
    repo_root: Path,
    dual_events: pd.DataFrame,
    selected_events: pd.DataFrame,
    file_inventory: pd.DataFrame,
    summary: str,
    manifest: dict[str, Any],
) -> None:
    """Write versioned local event-source evidence.

    Args:
        repo_root: Repository root.
        dual_events: Combined source events.
        selected_events: Selected source rows.
        file_inventory: File-level comparisons.
        summary: Markdown summary.
        manifest: Run manifest.

    Returns:
        None.

    Side effects:
        Creates OUTPUT_DIR and writes CSV, Markdown, and JSON files below
        ignored _Data. No raw file is touched.
    """

    output_dir = repo_root / OUTPUT_DIR
    validate_output_directory(output_dir, repo_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    dual_events.to_csv(
        output_dir / DUAL_SOURCE_EVENTS_FILENAME, index=False
    )
    selected_events.to_csv(
        output_dir / SELECTED_EVENTS_FILENAME, index=False
    )
    file_inventory.to_csv(
        output_dir / FILE_INVENTORY_FILENAME, index=False
    )
    (output_dir / SUMMARY_FILENAME).write_text(summary, encoding="utf-8")
    (output_dir / RUN_MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    """Run the dual-source inventory and write versioned local evidence.

    Args:
        None.

    Returns:
        None.

    Side effects:
        Reads all raw EDF annotations and Trigger channels, then writes local
        outputs below _Data/eeg/event_source_inventory_v1.
    """

    started_at = datetime.now()
    repo_root = repo_root_from_script()
    try:
        dual_events, selected_events, file_inventory = run_inventory(
            repo_root
        )
        summary = build_summary(
            file_inventory, dual_events, selected_events, started_at
        )
        manifest = build_run_manifest(
            repo_root=repo_root,
            file_inventory=file_inventory,
            dual_events=dual_events,
            selected_events=selected_events,
            started_at=started_at,
        )
        write_outputs(
            repo_root=repo_root,
            dual_events=dual_events,
            selected_events=selected_events,
            file_inventory=file_inventory,
            summary=summary,
            manifest=manifest,
        )
    except Exception as error:
        raise SystemExit(f"FAIL dual-source event inventory: {error}") from error

    print(
        f"Inspected {len(file_inventory)} EDF files; "
        f"wrote versioned event-source evidence to {repo_root / OUTPUT_DIR}"
    )


if __name__ == "__main__":
    main()
