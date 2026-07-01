"""Compare raw EDF annotations with old-compatible event offsets.

This script is the first raw EEG event-linkage inspection step in the DEMI EEG
reanalysis. It reads raw EDF annotations from ``_Data/eeg/raw/`` with MNE and
compares the annotation-derived trial/event structure with the old-compatible
behavioural/tracing event-offset table at
``_Data/behavior/event_offsets/event_offsets_old_compatible.csv``.

Raw annotation comparison matters because the new EEG workflow is restarting
from raw recordings. Before any preprocessing, epoching, filtering, or event
repair can be planned, the repository needs a transparent account of how the
raw EDF annotations line up with the reconstructed behavioural/tracing timing
table. This script records observable annotation facts, reconstructs trial
numbers using the non-preprocessing parts of the historical event logic, and
surfaces mismatches for review.

The historical reference point is ``_Scripts/_functions/eeg.R::update_events``.
That old R function operated after earlier EEG import/conversion steps had
already turned numeric trigger annotations into semantic labels such as
``stim_on``, ``red_on``, ``trace_start``, and ``trace_end``. This script works
one stage earlier, directly from raw EDF annotation descriptions such as
``28``, ``30``, ``44``, and ``46``. The numeric-to-semantic mapping is taken
from the historical Python EDF-to-BIDS helper in
``external/DEMI_EEG_Pipeline/edf2bids.py`` and is written into the audit output
as a derived, inspectable layer. The original raw annotation descriptions are
kept unchanged.

Old logic ported here:

- trial counting from ``stim_on`` annotations and bad-start-like ``red_on``
  annotations that are not immediately preceded by ``stim_on``;
- practice-trial-like detection from the old figure-duration rule
  (``red_on - stim_on > 3000 ms``);
- the old ordering of practice handling, trial renumbering, and bad-start
  handling in ``update_events``;
- per-trial checks for missing ``stim_on``, ``red_on``, ``trace_start``, and
  ``trace_end`` events;
- right-join-style comparison against old ``epoch_offsets`` keys through the
  old-compatible ``id`` and ``trial_count`` fields;
- old duration warnings where the raw annotations provide enough information:
  ``trace_end`` versus ``trace_start + end_trigger`` with a 150 ms threshold,
  and ``red_on - stim_on`` versus task ``stimulus_mt`` with a 250 ms threshold.

Work intentionally deferred:

- EEG filtering, epoching, bad-channel work, ICA, interpolation, CSD, or any
  other preprocessing;
- writing modified EDF/FIF/BIDS files;
- concatenating split EDFs or trimming the concatenated ID 5 recording;
- applying the historical preprocessing-side duplicate-trigger repairs;
- adding ``real_trace_start`` or ``real_trace_end`` annotations to EEG files;
- participant-retention decisions or scientific interpretation.

Expected local inputs:

- raw EDF files directly below ``_Data/eeg/raw/``;
- filenames that match the manifest script's DEMI EDF pattern, such as
  ``demi_54_1 Data.edf`` for split files;
- ``_Data/behavior/event_offsets/event_offsets_old_compatible.csv``;
- a Python environment with ``pandas`` and ``mne``.

Generated local outputs:

- ``_Data/eeg/event_alignment/raw_annotation_events.csv``: one row per raw EDF
  annotation, with original descriptions plus derived numeric/semantic labels;
- ``_Data/eeg/event_alignment/annotation_trial_summary.csv``: one row per
  reconstructed trial-like annotation group, plus factual file rows when no
  trial-like group can be reconstructed;
- ``_Data/eeg/event_alignment/offset_annotation_join_audit.csv``: a full outer
  key audit between raw annotation-derived trial rows and old-compatible
  behavioural/tracing offsets;
- ``_Data/eeg/event_alignment/raw_annotation_alignment_summary.md``: concise
  local summary of file counts, annotation-count flags, split/concatenated
  cases, join counts, duration-mismatch counts, and open parity TODOs.

Safety boundaries:

- EDFs are opened with ``mne.io.read_raw_edf(..., preload=False)``.
- Raw signal arrays are not loaded.
- Raw EDF files are not edited.
- Behavioural/event-offset inputs are read only.
- Outputs are written only below ``_Data/eeg/event_alignment/``.
- One problematic EDF records an error row and does not stop the run.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/01_compare_raw_annotations_to_offsets.py

Paths are resolved relative to this file, so the command also works from
another current working directory.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

# MNE imports matplotlib in some environments. Direct its cache to a writable
# temp path before importing MNE so this read-only audit does not depend on the
# user's global matplotlib cache directory.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "demi_matplotlib"))

import mne
import numpy as np
import pandas as pd


RAW_EEG_DIR = Path("_Data") / "eeg" / "raw"
OLD_COMPATIBLE_OFFSETS_PATH = Path("_Data") / "behavior" / "event_offsets" / "event_offsets_old_compatible.csv"
OUTPUT_DIR = Path("_Data") / "eeg" / "event_alignment"

RAW_ANNOTATION_EVENTS_FILENAME = "raw_annotation_events.csv"
ANNOTATION_TRIAL_SUMMARY_FILENAME = "annotation_trial_summary.csv"
OFFSET_ANNOTATION_JOIN_AUDIT_FILENAME = "offset_annotation_join_audit.csv"
ALIGNMENT_SUMMARY_FILENAME = "raw_annotation_alignment_summary.md"

LOW_ANNOTATION_COUNT_THRESHOLD = 600
PRACTICE_FIG_DURATION_SECONDS_THRESHOLD = 3.0
TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD = 0.150
STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD = 0.250

WATCHLIST_IDS = {5, 11, 14, 54, 56, 65, 86, 89, 94, 100}

DOMINANT_TASK_EVENT_CODES = (28, 30, 44, 46, 60)
EXPECTED_TRIAL_EVENT_NAMES = ("stim_on", "red_on", "trace_start", "trace_end")
TASK_EVENT_NAMES = (
    "stim_on",
    "red_on",
    "trace_start",
    "trace_end",
    "accuracy_submit",
    "vividness_submit",
)

# The raw EDF descriptions are numeric strings. The mapping below mirrors the
# historical Python EDF-to-BIDS helper that prepared files consumed by the old
# R EEG import path. Keeping this mapping explicit prevents the raw annotation
# description from being silently overwritten by a derived label.
EVENT_CODE_TO_NAME = {
    28: "stim_on",
    30: "red_on",
    44: "trace_start",
    46: "trace_end",
    60: "accuracy_submit",
    62: "vividness_submit",
}
SPECIAL_DESCRIPTION_TO_NAME = {"file start": "file start"}

# Reused from analysis/eeg_mne/00_make_raw_eeg_manifest.py so both raw EEG
# scripts agree on participant ID and split-file parsing.
EDF_FILENAME_RE = re.compile(
    r"^demi_(?P<participant_id>\d{1,3})(?:_(?P<split_part>\d+))?(?P<label>.*)\.edf$",
    flags=re.IGNORECASE,
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Args:
        None.

    Returns:
        Absolute path to the repository root.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def relative_to_repo(path: Path, repo_root: Path) -> str:
    """Render a path relative to the repository root when possible.

    Args:
        path: Path to render.
        repo_root: Repository root used as the preferred base.

    Returns:
        POSIX-style path text.

    Side effects:
        None.
    """

    try:
        return path.relative_to(repo_root).as_posix()
    except ValueError:
        return path.as_posix()


def discover_edf_files(raw_dir: Path) -> list[Path]:
    """Find raw EDF files directly below the raw EEG directory.

    Args:
        raw_dir: Directory expected to contain DEMI raw EDF files.

    Returns:
        Sorted list of EDF paths. The search is case-insensitive for the
        ``.edf`` suffix and is not recursive.

    Side effects:
        Reads directory entries only.
    """

    if not raw_dir.exists():
        return []
    return sorted(path for path in raw_dir.iterdir() if path.is_file() and path.suffix.lower() == ".edf")


def parse_edf_filename(edf_path: Path) -> dict[str, Any]:
    """Parse DEMI participant and file-role fields from an EDF filename.

    Args:
        edf_path: EDF path. Only ``edf_path.name`` is parsed.

    Returns:
        Dictionary with participant ID, padded ID, file role, split-part
        number, and parse-warning text. Unparsed filenames receive null
        participant fields and role ``unknown``.

    Side effects:
        None.
    """

    match = EDF_FILENAME_RE.match(edf_path.name)
    if match is None:
        return {
            "participant_id": pd.NA,
            "participant_id_padded": "",
            "file_role": "unknown",
            "split_part": pd.NA,
            "filename_parse_warning": "filename_does_not_match_expected_demi_pattern",
        }

    participant_id = int(match.group("participant_id"))
    split_part_text = match.group("split_part")
    split_part = int(split_part_text) if split_part_text is not None else pd.NA
    label_text = match.group("label").lower()

    if split_part_text is not None:
        file_role = "split_part"
    elif "concatenated" in label_text:
        file_role = "concatenated"
    else:
        file_role = "single"

    return {
        "participant_id": participant_id,
        "participant_id_padded": f"{participant_id:03d}",
        "file_role": file_role,
        "split_part": split_part,
        "filename_parse_warning": "",
    }


def numeric_event_code(description: Any) -> int | None:
    """Parse an annotation description as an integer event code when possible.

    Args:
        description: Raw MNE annotation description.

    Returns:
        Integer code for plain whole-number descriptions, otherwise ``None``.

    Side effects:
        None.
    """

    text = str(description).strip()
    if re.fullmatch(r"[+-]?\d+", text) is None:
        return None
    return int(text)


def event_name_from_description(description: Any) -> str:
    """Map a raw annotation description to a semantic event name when known.

    Args:
        description: Raw MNE annotation description.

    Returns:
        Semantic event name from the historical mapping, or an empty string
        when the description is not part of that mapping.

    Side effects:
        None.
    """

    text = str(description).strip()
    lower_text = text.lower()
    if lower_text in SPECIAL_DESCRIPTION_TO_NAME:
        return SPECIAL_DESCRIPTION_TO_NAME[lower_text]

    code = numeric_event_code(text)
    if code is None:
        return ""
    return EVENT_CODE_TO_NAME.get(code, "")


def is_watchlist_id(participant_id: Any) -> bool:
    """Return whether a participant ID is on the event-alignment watchlist.

    Args:
        participant_id: Participant ID value that may be missing.

    Returns:
        ``True`` for IDs named in the current raw-annotation watchlist.

    Side effects:
        None.
    """

    if pd.isna(participant_id):
        return False
    return int(participant_id) in WATCHLIST_IDS


def empty_annotation_frame() -> pd.DataFrame:
    """Create an empty raw-annotation event table with stable columns.

    Args:
        None.

    Returns:
        Empty DataFrame using the same columns as non-empty raw annotation
        output.

    Side effects:
        None.
    """

    return pd.DataFrame(
        columns=[
            "participant_id",
            "participant_id_padded",
            "source_filename",
            "file_path",
            "file_role",
            "split_part",
            "watchlist_id",
            "read_status",
            "sampling_frequency_hz",
            "recording_duration_seconds",
            "mne_orig_time",
            "event_row_order",
            "onset_seconds",
            "duration_seconds",
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
        ]
    )


def read_annotations_for_edf(edf_path: Path, repo_root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read raw EDF annotations from one file without preloading signal data.

    Args:
        edf_path: Raw EDF file path.
        repo_root: Repository root used for relative provenance paths.

    Returns:
        ``(events, file_row)``. ``events`` has one row per annotation and may
        be empty for zero-annotation or read-error files. ``file_row`` contains
        file-level counts, MNE read status, and error details.

    Side effects:
        Opens the EDF through MNE with ``preload=False``. The Raw object is
        closed before returning when MNE exposes ``close()``.
    """

    parsed = parse_edf_filename(edf_path)
    participant_id = parsed["participant_id"]
    base_file_row: dict[str, Any] = {
        **parsed,
        "source_filename": edf_path.name,
        "file_path": relative_to_repo(edf_path, repo_root),
        "watchlist_id": is_watchlist_id(participant_id),
        "read_status": "not_read",
        "mne_read_error_type": "",
        "mne_read_error_message": "",
        "sampling_frequency_hz": np.nan,
        "recording_duration_seconds": np.nan,
        "mne_orig_time": "",
        "annotation_count": 0,
        "mapped_task_event_count": 0,
        "trial_numbering_candidate_count": 0,
        "unmapped_annotation_count": 0,
        "dominant_event_code_pattern": "",
        "annotation_description_counts_json": "{}",
        "file_issue_codes": "",
    }

    raw: mne.io.BaseRaw | None = None
    try:
        # MNE stores EDF annotations on the Raw object header side. With
        # preload=False, this reads metadata and annotation timing without
        # loading the EEG signal array into memory.
        raw = mne.io.read_raw_edf(edf_path, preload=False, verbose="ERROR")
        sampling_frequency_hz = float(raw.info["sfreq"])
        recording_duration_seconds = float(raw.n_times / sampling_frequency_hz) if sampling_frequency_hz else np.nan
        orig_time = raw.annotations.orig_time.isoformat() if raw.annotations.orig_time is not None else ""

        rows: list[dict[str, Any]] = []
        for event_row_order, annotation in enumerate(raw.annotations, start=1):
            description = str(annotation["description"])
            event_code = numeric_event_code(description)
            event_name = event_name_from_description(description)
            mapped_task_event = event_name in TASK_EVENT_NAMES
            rows.append(
                {
                    "participant_id": participant_id,
                    "participant_id_padded": parsed["participant_id_padded"],
                    "source_filename": edf_path.name,
                    "file_path": relative_to_repo(edf_path, repo_root),
                    "file_role": parsed["file_role"],
                    "split_part": parsed["split_part"],
                    "watchlist_id": is_watchlist_id(participant_id),
                    "read_status": "ok",
                    "sampling_frequency_hz": sampling_frequency_hz,
                    "recording_duration_seconds": recording_duration_seconds,
                    "mne_orig_time": orig_time,
                    "event_row_order": event_row_order,
                    "onset_seconds": float(annotation["onset"]),
                    "duration_seconds": float(annotation["duration"]),
                    "annotation_description": description,
                    "event_code": event_code if event_code is not None else pd.NA,
                    "event_name": event_name,
                    "mapped_task_event": mapped_task_event,
                    # The old R function computed trial groups from semantic
                    # trigger labels. For raw EDFs, only mapped task events are
                    # candidates for that reconstruction; annotations such as
                    # impedance checks remain in the event table but do not
                    # affect trial numbering.
                    "trial_numbering_candidate": mapped_task_event,
                    "raw_trial_sequence": pd.NA,
                    "old_update_trial_count": pd.NA,
                    "join_trial_count": pd.NA,
                    "old_update_removed_reason": "",
                    "fig_duration_seconds": np.nan,
                    "bad_start_like_event": False,
                    "practice_like_group": False,
                    "bad_start_like_group": False,
                }
            )

        events = pd.DataFrame(rows, columns=empty_annotation_frame().columns)
        description_counts = Counter(events["annotation_description"]) if not events.empty else Counter()
        event_name_counts = Counter(events["event_name"]) if not events.empty else Counter()
        dominant_event_code_pattern = format_event_pattern(event_name_counts)

        file_row = {
            **base_file_row,
            "read_status": "ok",
            "sampling_frequency_hz": sampling_frequency_hz,
            "recording_duration_seconds": recording_duration_seconds,
            "mne_orig_time": orig_time,
            "annotation_count": int(len(events)),
            "mapped_task_event_count": int(events["mapped_task_event"].sum()) if not events.empty else 0,
            "trial_numbering_candidate_count": int(events["trial_numbering_candidate"].sum()) if not events.empty else 0,
            "unmapped_annotation_count": int((events["event_name"].eq("")).sum()) if not events.empty else 0,
            "dominant_event_code_pattern": dominant_event_code_pattern,
            "annotation_description_counts_json": json.dumps(dict(sorted(description_counts.items())), sort_keys=True),
        }
        file_row["file_issue_codes"] = ";".join(file_issue_codes(file_row))
        return events, file_row

    except Exception as error:  # noqa: BLE001 - one EDF should not stop the audit.
        file_row = {
            **base_file_row,
            "read_status": "mne_read_error",
            "mne_read_error_type": type(error).__name__,
            "mne_read_error_message": str(error),
        }
        file_row["file_issue_codes"] = ";".join(file_issue_codes(file_row))
        return empty_annotation_frame(), file_row
    finally:
        if raw is not None and hasattr(raw, "close"):
            raw.close()


def format_event_pattern(event_name_counts: Counter[str]) -> str:
    """Format key mapped-event counts for compact file-level reporting.

    Args:
        event_name_counts: Counter keyed by derived event name.

    Returns:
        Semicolon-delimited count text for mapped task events and other mapped
        annotations.

    Side effects:
        None.
    """

    parts = [f"{event_name}={event_name_counts.get(event_name, 0)}" for event_name in TASK_EVENT_NAMES]
    special_count = event_name_counts.get("file start", 0)
    if special_count:
        parts.append(f"file start={special_count}")
    return ";".join(parts)


def file_issue_codes(file_row: dict[str, Any]) -> list[str]:
    """Create factual file-level audit codes.

    Args:
        file_row: File-level annotation/read metadata.

    Returns:
        List of short issue/status codes. Codes report observable facts and do
        not make participant-retention decisions.

    Side effects:
        None.
    """

    codes: list[str] = []
    if file_row["read_status"] == "mne_read_error":
        codes.append("mne_read_error")
    if file_row["annotation_count"] == 0:
        codes.append("zero_annotations")
    elif file_row["annotation_count"] < LOW_ANNOTATION_COUNT_THRESHOLD:
        codes.append("low_annotation_count")
    if file_row["file_role"] == "split_part":
        codes.append("split_file_not_concatenated")
    if file_row["file_role"] == "concatenated":
        codes.append("concatenated_file_not_repaired")
    if file_row["watchlist_id"]:
        codes.append("watchlist_id")
    if file_row["filename_parse_warning"]:
        codes.append("filename_parse_warning")

    description_counts = json.loads(file_row.get("annotation_description_counts_json", "{}"))
    for code in DOMINANT_TASK_EVENT_CODES:
        if str(code) not in description_counts:
            codes.append(f"missing_dominant_code_{code}")

    return codes


def read_all_raw_annotations(edf_files: list[Path], repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read annotations from all discovered EDF files.

    Args:
        edf_files: EDF paths to inspect.
        repo_root: Repository root used for relative paths.

    Returns:
        ``(raw_events, file_summary)``. ``raw_events`` has one row per raw
        annotation. ``file_summary`` has one row per EDF file.

    Side effects:
        Opens each EDF through MNE without loading raw signal arrays.
    """

    event_frames: list[pd.DataFrame] = []
    file_rows: list[dict[str, Any]] = []
    for edf_path in edf_files:
        events, file_row = read_annotations_for_edf(edf_path, repo_root)
        event_frames.append(events)
        file_rows.append(file_row)

    if event_frames:
        raw_events = pd.concat(event_frames, ignore_index=True, sort=False)
    else:
        raw_events = empty_annotation_frame()
    file_summary = pd.DataFrame(file_rows)
    return raw_events, file_summary


def first_nonmissing(values: pd.Series) -> Any:
    """Return the first nonmissing value from a Series.

    Args:
        values: Series of candidate values.

    Returns:
        First nonmissing value, or ``NaN`` when no value is present.

    Side effects:
        None.
    """

    nonmissing = values.dropna()
    if nonmissing.empty:
        return np.nan
    return nonmissing.iloc[0]


def is_true(value: Any) -> bool:
    """Return ``True`` only for explicit true-like scalar values.

    Args:
        value: Scalar value from a row or Series.

    Returns:
        ``True`` for explicit true-like values. Missing values return
        ``False``.

    Side effects:
        None.
    """

    if pd.isna(value):
        return False
    return bool(value)


def unique_joined_text(values: pd.Series) -> str:
    """Return unique non-empty values joined by semicolons.

    Args:
        values: Series of values that may be missing or empty.

    Returns:
        Semicolon-delimited text in first-seen order.

    Side effects:
        None.
    """

    out: list[str] = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value)
        if text and text not in out:
            out.append(text)
    return ";".join(out)


def assign_old_update_trial_numbers(raw_events: pd.DataFrame) -> pd.DataFrame:
    """Apply old ``update_events`` trial-numbering logic to raw annotations.

    Args:
        raw_events: One-row-per-annotation table produced by
            ``read_all_raw_annotations``.

    Returns:
        Copy of ``raw_events`` with raw trial sequence, old post-practice trial
        number, join trial number, practice-like flags, and bad-start-like
        flags filled for mapped task-event rows.

    Side effects:
        None.
    """

    if raw_events.empty:
        return raw_events.copy()

    out = raw_events.copy()
    group_columns = ["participant_id", "source_filename"]

    for _, index in out.groupby(group_columns, dropna=False, sort=False).groups.items():
        file_events = out.loc[index].sort_values("event_row_order").copy()
        candidate = file_events[file_events["trial_numbering_candidate"].fillna(False)].copy()
        if candidate.empty:
            continue

        candidate["original_event_index"] = candidate.index
        candidate["previous_task_event_name"] = candidate["event_name"].shift()
        candidate["previous_task_onset_seconds"] = candidate["onset_seconds"].shift()

        # Old R marked a bad start when red_on appeared without an immediately
        # preceding stim_on. That event still starts a trial group, preserving
        # the old trial counter's awareness of a missing stim trigger.
        candidate["bad_start_like_event"] = (
            candidate["event_name"].eq("red_on") & ~candidate["previous_task_event_name"].eq("stim_on")
        )
        starts_new_trial = candidate["event_name"].eq("stim_on") | candidate["bad_start_like_event"]
        candidate["raw_trial_sequence"] = starts_new_trial.cumsum()

        candidate["duration_since_previous_task_event_seconds"] = (
            candidate["onset_seconds"] - candidate["previous_task_onset_seconds"]
        )
        candidate.loc[candidate["bad_start_like_event"], "duration_since_previous_task_event_seconds"] = np.nan

        trial_rows: list[dict[str, Any]] = []
        for raw_trial_sequence, trial_events in candidate.groupby("raw_trial_sequence", sort=True):
            red_durations = trial_events.loc[
                trial_events["event_name"].eq("red_on"), "duration_since_previous_task_event_seconds"
            ].dropna()
            fig_duration_seconds = float(red_durations.iloc[0]) if not red_durations.empty else np.nan
            practice_like_group = (
                not pd.isna(fig_duration_seconds)
                and fig_duration_seconds > PRACTICE_FIG_DURATION_SECONDS_THRESHOLD
            )
            bad_start_like_group = bool(trial_events["bad_start_like_event"].any())
            trial_rows.append(
                {
                    "raw_trial_sequence": int(raw_trial_sequence),
                    "fig_duration_seconds": fig_duration_seconds,
                    "practice_like_group": practice_like_group,
                    "bad_start_like_group": bad_start_like_group,
                }
            )

        trial_info = pd.DataFrame(trial_rows)

        # This mirrors the order in update_events: practice-like groups are
        # removed first, then the remaining groups are renumbered, then
        # bad-start-like groups are removed from the event stream. The join
        # trial count is therefore missing for bad-start groups, but later
        # groups keep the old post-practice numbering.
        nonpractice = trial_info[~trial_info["practice_like_group"]].copy()
        nonpractice["old_update_trial_count"] = np.arange(1, len(nonpractice) + 1, dtype="int64")
        trial_info = trial_info.merge(
            nonpractice[["raw_trial_sequence", "old_update_trial_count"]],
            on="raw_trial_sequence",
            how="left",
        )
        trial_info["join_trial_count"] = trial_info["old_update_trial_count"]
        trial_info.loc[trial_info["bad_start_like_group"], "join_trial_count"] = np.nan

        trial_info["old_update_removed_reason"] = ""
        trial_info.loc[trial_info["practice_like_group"], "old_update_removed_reason"] = "practice_like_group"
        trial_info.loc[trial_info["bad_start_like_group"], "old_update_removed_reason"] = "bad_start_like_group"
        trial_info.loc[
            trial_info["practice_like_group"] & trial_info["bad_start_like_group"],
            "old_update_removed_reason",
        ] = "practice_like_group;bad_start_like_group"

        original_index = candidate["original_event_index"]
        trial_info_by_sequence = trial_info.set_index("raw_trial_sequence")
        for column in [
            "raw_trial_sequence",
            "old_update_trial_count",
            "join_trial_count",
            "old_update_removed_reason",
            "fig_duration_seconds",
            "practice_like_group",
            "bad_start_like_group",
            "bad_start_like_event",
        ]:
            if column in {"raw_trial_sequence", "bad_start_like_event"}:
                values = candidate[column]
            else:
                values = candidate["raw_trial_sequence"].map(trial_info_by_sequence[column])
            out.loc[original_index, column] = values.to_numpy()

    return coerce_nullable_integer_columns(out)


def coerce_nullable_integer_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert known integer-like audit columns to pandas nullable integers.

    Args:
        frame: DataFrame with integer-like columns that may have missing values.

    Returns:
        Copy with stable nullable integer dtypes where columns are present.

    Side effects:
        None.
    """

    out = frame.copy()
    for column in ["participant_id", "split_part", "event_code", "raw_trial_sequence", "old_update_trial_count", "join_trial_count"]:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
    return out


def summarize_trial_events(trial_events: pd.DataFrame) -> dict[str, Any]:
    """Summarize one raw annotation-derived trial group.

    Args:
        trial_events: Mapped task-event rows from one EDF and raw trial group.

    Returns:
        Dictionary containing per-trial event counts, key event onsets, missing
        event flags, and duration fields.

    Side effects:
        None.
    """

    counts = Counter(trial_events["event_name"])
    row: dict[str, Any] = {
        "summary_row_kind": "trial_group",
        "participant_id": first_nonmissing(trial_events["participant_id"]),
        "participant_id_padded": first_nonmissing(trial_events["participant_id_padded"]),
        "source_filename": first_nonmissing(trial_events["source_filename"]),
        "file_path": first_nonmissing(trial_events["file_path"]),
        "file_role": first_nonmissing(trial_events["file_role"]),
        "split_part": first_nonmissing(trial_events["split_part"]),
        "watchlist_id": is_true(first_nonmissing(trial_events["watchlist_id"])),
        "read_status": first_nonmissing(trial_events["read_status"]),
        "sampling_frequency_hz": first_nonmissing(trial_events["sampling_frequency_hz"]),
        "recording_duration_seconds": first_nonmissing(trial_events["recording_duration_seconds"]),
        "annotation_count": np.nan,
        "mapped_task_event_count": int(len(trial_events)),
        "dominant_event_code_pattern": "",
        "raw_trial_sequence": first_nonmissing(trial_events["raw_trial_sequence"]),
        "old_update_trial_count": first_nonmissing(trial_events["old_update_trial_count"]),
        "join_trial_count": first_nonmissing(trial_events["join_trial_count"]),
        "old_update_removed_reason": unique_joined_text(trial_events["old_update_removed_reason"]),
        "practice_like_group": is_true(first_nonmissing(trial_events["practice_like_group"])),
        "bad_start_like_group": is_true(first_nonmissing(trial_events["bad_start_like_group"])),
        "fig_duration_seconds": first_nonmissing(trial_events["fig_duration_seconds"]),
        "trial_event_count": int(len(trial_events)),
    }

    for event_name in TASK_EVENT_NAMES:
        event_rows = trial_events[trial_events["event_name"].eq(event_name)]
        row[f"{event_name}_count"] = int(counts.get(event_name, 0))
        row[f"{event_name}_onset_seconds"] = first_nonmissing(event_rows["onset_seconds"]) if not event_rows.empty else np.nan

    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        row[f"missing_{event_name}"] = row[f"{event_name}_count"] == 0

    missing_expected = [event_name for event_name in EXPECTED_TRIAL_EVENT_NAMES if row[f"missing_{event_name}"]]
    row["missing_expected_event_names"] = ";".join(missing_expected)

    trace_start = row["trace_start_onset_seconds"]
    trace_end = row["trace_end_onset_seconds"]
    if pd.notna(trace_start) and pd.notna(trace_end):
        row["trace_epoch_duration_seconds"] = float(trace_end) - float(trace_start)
    else:
        row["trace_epoch_duration_seconds"] = np.nan

    row["trial_issue_codes"] = ";".join(trial_issue_codes(row))
    return row


def file_without_trials_summary_row(file_row: pd.Series) -> dict[str, Any]:
    """Create a trial-summary row for files with no reconstructed trial group.

    Args:
        file_row: One row from the file-level annotation summary.

    Returns:
        Dictionary matching the trial-summary schema as closely as possible,
        with trial fields left missing.

    Side effects:
        None.
    """

    row: dict[str, Any] = {
        "summary_row_kind": "file_without_reconstructed_trial_groups",
        "participant_id": file_row["participant_id"],
        "participant_id_padded": file_row["participant_id_padded"],
        "source_filename": file_row["source_filename"],
        "file_path": file_row["file_path"],
        "file_role": file_row["file_role"],
        "split_part": file_row["split_part"],
        "watchlist_id": is_true(file_row["watchlist_id"]),
        "read_status": file_row["read_status"],
        "sampling_frequency_hz": file_row["sampling_frequency_hz"],
        "recording_duration_seconds": file_row["recording_duration_seconds"],
        "annotation_count": file_row["annotation_count"],
        "mapped_task_event_count": file_row["mapped_task_event_count"],
        "dominant_event_code_pattern": file_row["dominant_event_code_pattern"],
        "raw_trial_sequence": pd.NA,
        "old_update_trial_count": pd.NA,
        "join_trial_count": pd.NA,
        "old_update_removed_reason": "",
        "practice_like_group": False,
        "bad_start_like_group": False,
        "fig_duration_seconds": np.nan,
        "trial_event_count": 0,
        "trace_epoch_duration_seconds": np.nan,
        "missing_expected_event_names": ";".join(EXPECTED_TRIAL_EVENT_NAMES),
        "trial_issue_codes": file_row["file_issue_codes"],
    }
    for event_name in TASK_EVENT_NAMES:
        row[f"{event_name}_count"] = 0
        row[f"{event_name}_onset_seconds"] = np.nan
    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        row[f"missing_{event_name}"] = True
    return row


def trial_issue_codes(row: dict[str, Any]) -> list[str]:
    """Create per-trial annotation audit codes.

    Args:
        row: Trial-summary row dictionary.

    Returns:
        List of short issue/status codes.

    Side effects:
        None.
    """

    codes: list[str] = []
    if is_true(row.get("practice_like_group", False)):
        codes.append("practice_like_group")
    if is_true(row.get("bad_start_like_group", False)):
        codes.append("bad_start_like_group")
    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        if row.get(f"missing_{event_name}", False):
            codes.append(f"missing_{event_name}")
    if row.get("file_role") == "split_part":
        codes.append("split_file_not_concatenated")
    if row.get("file_role") == "concatenated":
        codes.append("concatenated_file_not_repaired")
    if is_true(row.get("watchlist_id", False)):
        codes.append("watchlist_id")
    return codes


def build_annotation_trial_summary(raw_events: pd.DataFrame, file_summary: pd.DataFrame) -> pd.DataFrame:
    """Build one trial-level annotation summary table.

    Args:
        raw_events: Raw annotation event table after old trial-number
            assignment.
        file_summary: One-row-per-EDF file summary.

    Returns:
        DataFrame with one row per reconstructed raw trial group and a
        placeholder file row for EDFs without such groups.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    if not raw_events.empty:
        candidates = raw_events[raw_events["trial_numbering_candidate"].fillna(False)].copy()
        candidates = candidates[candidates["raw_trial_sequence"].notna()]
        for _, trial_events in candidates.groupby(
            ["participant_id", "source_filename", "raw_trial_sequence"],
            dropna=False,
            sort=True,
        ):
            rows.append(summarize_trial_events(trial_events.sort_values("event_row_order")))

    files_with_trial_rows = {row["source_filename"] for row in rows}
    for _, file_row in file_summary.iterrows():
        if file_row["source_filename"] not in files_with_trial_rows:
            rows.append(file_without_trials_summary_row(file_row))

    trial_summary = pd.DataFrame(rows)
    if trial_summary.empty:
        return trial_summary

    for _, file_row in file_summary.iterrows():
        mask = trial_summary["source_filename"].eq(file_row["source_filename"])
        for column in [
            "annotation_count",
            "mapped_task_event_count",
            "dominant_event_code_pattern",
            "recording_duration_seconds",
            "sampling_frequency_hz",
        ]:
            trial_summary.loc[mask, column] = file_row[column]

    sort_columns = ["participant_id", "source_filename", "raw_trial_sequence"]
    trial_summary = trial_summary.sort_values(sort_columns, kind="mergesort", na_position="last")
    return coerce_nullable_integer_columns(trial_summary.reset_index(drop=True))


def read_old_compatible_offsets(path: Path) -> pd.DataFrame:
    """Read and validate the old-compatible behavioural/tracing offset table.

    Args:
        path: CSV path for ``event_offsets_old_compatible.csv``.

    Returns:
        DataFrame with old ``epoch_offsets`` key and timing fields.

    Side effects:
        Reads the CSV from disk. Raises ``RuntimeError`` if required columns
        are missing.
    """

    required_columns = {
        "id",
        "physical",
        "it",
        "mt",
        "stimulus_mt",
        "trial_count",
        "trace_onset",
        "trace_end",
        "start_shift",
        "real_start",
        "real_end",
        "end_trigger",
    }
    if not path.exists():
        raise RuntimeError(f"required old-compatible offset file is missing: {path.as_posix()}")

    offsets = pd.read_csv(path)
    missing_columns = sorted(required_columns.difference(offsets.columns))
    if missing_columns:
        raise RuntimeError(f"old-compatible offset file is missing column(s): {', '.join(missing_columns)}")

    offsets = offsets.copy()
    offsets["id"] = pd.to_numeric(offsets["id"], errors="coerce").astype("Int64")
    offsets["trial_count"] = pd.to_numeric(offsets["trial_count"], errors="coerce").astype("Int64")
    duplicate_keys = offsets.duplicated(["id", "trial_count"], keep=False)
    if duplicate_keys.any():
        duplicate_count = int(duplicate_keys.sum())
        raise RuntimeError(f"old-compatible offsets contain {duplicate_count} duplicate id/trial_count row(s)")

    return offsets.sort_values(["id", "trial_count"], kind="mergesort").reset_index(drop=True)


def participant_file_fact_maps(file_summary: pd.DataFrame) -> dict[str, set[int]]:
    """Create participant-level file fact maps for join-audit issue codes.

    Args:
        file_summary: One-row-per-EDF file summary.

    Returns:
        Dictionary mapping fact names to participant ID sets.

    Side effects:
        None.
    """

    maps = {
        "has_raw_edf_file": set(),
        "has_zero_annotation_file": set(),
        "has_low_annotation_file": set(),
        "has_split_file": set(),
        "has_concatenated_file": set(),
    }
    for _, row in file_summary.iterrows():
        participant_id = row["participant_id"]
        if pd.isna(participant_id):
            continue
        pid = int(participant_id)
        maps["has_raw_edf_file"].add(pid)
        if int(row["annotation_count"]) == 0:
            maps["has_zero_annotation_file"].add(pid)
        if 0 < int(row["annotation_count"]) < LOW_ANNOTATION_COUNT_THRESHOLD:
            maps["has_low_annotation_file"].add(pid)
        if row["file_role"] == "split_part":
            maps["has_split_file"].add(pid)
        if row["file_role"] == "concatenated":
            maps["has_concatenated_file"].add(pid)
    return maps


def build_offset_annotation_join_audit(
    trial_summary: pd.DataFrame,
    offsets: pd.DataFrame,
    file_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Compare raw annotation-derived trials with old-compatible offsets.

    Args:
        trial_summary: Per-trial annotation summary.
        offsets: Old-compatible behavioural/tracing offsets.
        file_summary: One-row-per-EDF file summary used for participant-level
            file facts on offset-only rows.

    Returns:
        Full outer join audit keyed by participant ID and trial count.

    Side effects:
        None.
    """

    raw_joinable = trial_summary[
        trial_summary["summary_row_kind"].eq("trial_group") & trial_summary["join_trial_count"].notna()
    ].copy()

    raw_joinable["participant_id"] = pd.to_numeric(raw_joinable["participant_id"], errors="coerce").astype("Int64")
    raw_joinable["join_trial_count"] = pd.to_numeric(raw_joinable["join_trial_count"], errors="coerce").astype("Int64")
    raw_joinable["raw_join_key_duplicate"] = raw_joinable.duplicated(
        ["participant_id", "join_trial_count"], keep=False
    )

    raw_columns = [
        "participant_id",
        "participant_id_padded",
        "source_filename",
        "file_path",
        "file_role",
        "split_part",
        "watchlist_id",
        "summary_row_kind",
        "raw_trial_sequence",
        "old_update_trial_count",
        "join_trial_count",
        "raw_join_key_duplicate",
        "practice_like_group",
        "bad_start_like_group",
        "fig_duration_seconds",
        "trial_event_count",
        "stim_on_count",
        "red_on_count",
        "trace_start_count",
        "trace_end_count",
        "accuracy_submit_count",
        "vividness_submit_count",
        "stim_on_onset_seconds",
        "red_on_onset_seconds",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
        "accuracy_submit_onset_seconds",
        "vividness_submit_onset_seconds",
        "trace_epoch_duration_seconds",
        "missing_expected_event_names",
        "trial_issue_codes",
    ]
    raw_joinable = raw_joinable[[column for column in raw_columns if column in raw_joinable.columns]].copy()

    offsets_for_join = offsets.rename(columns={"id": "offset_id", "trial_count": "offset_trial_count"})

    # The old R update used offset keys as the authority surface after trial
    # reconstruction. A full outer merge keeps both failure directions visible:
    # raw annotation-derived trials that have no old-compatible offset, and
    # old-compatible offsets that have no raw annotation-derived trial.
    joined = raw_joinable.merge(
        offsets_for_join,
        left_on=["participant_id", "join_trial_count"],
        right_on=["offset_id", "offset_trial_count"],
        how="outer",
        validate="many_to_one",
        indicator=True,
    )

    joined["audit_participant_id"] = joined["participant_id"].combine_first(joined["offset_id"]).astype("Int64")
    joined["audit_trial_count"] = joined["join_trial_count"].combine_first(joined["offset_trial_count"]).astype("Int64")
    joined["join_status"] = joined["_merge"].map(
        {
            "both": "raw_annotation_and_offset",
            "left_only": "raw_annotation_without_offset",
            "right_only": "offset_without_raw_annotation_trial",
        }
    )

    # Old update_events warned when the observed trace_end trigger was not close
    # to trace_start + end_trigger. Work in seconds here because MNE annotation
    # onsets are already seconds; the threshold is the old 150 ms rule.
    joined["trace_end_expected_from_offset_seconds"] = joined["trace_start_onset_seconds"] + joined["end_trigger"]
    joined["trace_epoch_end_difference_seconds"] = (
        joined["trace_end_onset_seconds"] - joined["trace_end_expected_from_offset_seconds"]
    ).abs()
    joined["old_trace_epoch_duration_mismatch"] = (
        joined["trace_epoch_end_difference_seconds"] > TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD
    ).fillna(False)

    # The old stimulus-duration warning compared red_on - stim_on with the task
    # stimulus_mt field. This catches event-sequence drift before any EEG
    # preprocessing decisions are made.
    joined["stimulus_duration_difference_seconds"] = joined["fig_duration_seconds"] - joined["stimulus_mt"]
    joined["old_stimulus_duration_mismatch"] = (
        joined["stimulus_duration_difference_seconds"].abs() > STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD
    ).fillna(False)

    fact_maps = participant_file_fact_maps(file_summary)
    joined["participant_has_raw_edf_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_raw_edf_file"] if pd.notna(value) else False
    )
    joined["participant_has_zero_annotation_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_zero_annotation_file"] if pd.notna(value) else False
    )
    joined["participant_has_low_annotation_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_low_annotation_file"] if pd.notna(value) else False
    )
    joined["participant_has_split_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_split_file"] if pd.notna(value) else False
    )
    joined["participant_has_concatenated_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_concatenated_file"] if pd.notna(value) else False
    )
    joined["watchlist_id"] = joined["audit_participant_id"].map(is_watchlist_id)
    joined["join_audit_issue_codes"] = joined.apply(join_issue_codes, axis=1)

    preferred = [
        "audit_participant_id",
        "audit_trial_count",
        "join_status",
        "join_audit_issue_codes",
        "source_filename",
        "file_role",
        "split_part",
        "watchlist_id",
        "participant_has_raw_edf_file",
        "participant_has_zero_annotation_file",
        "participant_has_low_annotation_file",
        "participant_has_split_file",
        "participant_has_concatenated_file",
        "raw_trial_sequence",
        "old_update_trial_count",
        "join_trial_count",
        "offset_id",
        "offset_trial_count",
        "raw_join_key_duplicate",
        "fig_duration_seconds",
        "stimulus_mt",
        "stimulus_duration_difference_seconds",
        "old_stimulus_duration_mismatch",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
        "trace_epoch_duration_seconds",
        "end_trigger",
        "trace_end_expected_from_offset_seconds",
        "trace_epoch_end_difference_seconds",
        "old_trace_epoch_duration_mismatch",
        "physical",
        "it",
        "mt",
        "trace_onset",
        "trace_end",
        "start_shift",
        "real_start",
        "real_end",
        "trial_event_count",
        "stim_on_count",
        "red_on_count",
        "trace_start_count",
        "trace_end_count",
        "accuracy_submit_count",
        "vividness_submit_count",
        "missing_expected_event_names",
        "trial_issue_codes",
    ]
    ordered = [column for column in preferred if column in joined.columns] + [
        column for column in joined.columns if column not in preferred and column != "_merge"
    ]
    return coerce_nullable_integer_columns(
        joined[ordered].sort_values(
            ["audit_participant_id", "audit_trial_count", "source_filename"],
            kind="mergesort",
            na_position="last",
        )
    ).reset_index(drop=True)


def join_issue_codes(row: pd.Series) -> str:
    """Build semicolon-delimited issue codes for one raw/offset join row.

    Args:
        row: Row from the full outer join audit.

    Returns:
        Semicolon-delimited issue/status codes.

    Side effects:
        None.
    """

    codes: list[str] = []
    if row.get("join_status") == "raw_annotation_without_offset":
        codes.append("raw_annotation_without_offset")
    if row.get("join_status") == "offset_without_raw_annotation_trial":
        codes.append("offset_without_raw_annotation_trial")
    if is_true(row.get("raw_join_key_duplicate", False)):
        codes.append("duplicate_raw_annotation_join_key")
    if is_true(row.get("old_trace_epoch_duration_mismatch", False)):
        codes.append("old_trace_epoch_duration_mismatch")
    if is_true(row.get("old_stimulus_duration_mismatch", False)):
        codes.append("old_stimulus_duration_mismatch")

    missing_expected = row.get("missing_expected_event_names", "")
    if isinstance(missing_expected, str) and missing_expected:
        for event_name in missing_expected.split(";"):
            codes.append(f"missing_{event_name}")

    if is_true(row.get("participant_has_zero_annotation_file", False)):
        codes.append("participant_has_zero_annotation_file")
    if is_true(row.get("participant_has_low_annotation_file", False)):
        codes.append("participant_has_low_annotation_file")
    if is_true(row.get("participant_has_split_file", False)):
        codes.append("split_file_not_concatenated")
    if is_true(row.get("participant_has_concatenated_file", False)):
        codes.append("concatenated_file_not_repaired")
    if not is_true(row.get("participant_has_raw_edf_file", True)):
        codes.append("participant_has_no_raw_edf_file")
    if is_true(row.get("watchlist_id", False)):
        codes.append("watchlist_id")

    return ";".join(dict.fromkeys(code for code in codes if code))


def count_issue_codes(values: pd.Series) -> Counter[str]:
    """Count semicolon-delimited issue codes.

    Args:
        values: Series of issue-code strings.

    Returns:
        Counter keyed by issue code.

    Side effects:
        None.
    """

    counts: Counter[str] = Counter()
    for value in values.fillna(""):
        for code in str(value).split(";"):
            if code:
                counts[code] += 1
    return counts


def format_count_table(counts: pd.Series) -> list[str]:
    """Format a pandas count Series as Markdown bullets.

    Args:
        counts: Count Series indexed by label.

    Returns:
        Markdown bullet lines.

    Side effects:
        None.
    """

    if counts.empty:
        return ["- none"]
    return [f"- {label}: {int(count)}" for label, count in counts.items()]


def format_id_list(values: pd.Series | list[Any] | set[Any]) -> str:
    """Format participant IDs as compact comma-delimited text.

    Args:
        values: Values that may contain missing data.

    Returns:
        Sorted comma-delimited integer ID text, or ``none``.

    Side effects:
        None.
    """

    ids = sorted({int(value) for value in values if pd.notna(value)})
    return ", ".join(str(value) for value in ids) if ids else "none"


def build_alignment_summary(
    raw_events: pd.DataFrame,
    file_summary: pd.DataFrame,
    trial_summary: pd.DataFrame,
    join_audit: pd.DataFrame,
    offsets: pd.DataFrame,
    raw_dir: Path,
    offsets_path: Path,
    started_at: datetime,
) -> str:
    """Build the local Markdown summary for raw annotation alignment.

    Args:
        raw_events: One-row-per-annotation table.
        file_summary: One-row-per-EDF summary.
        trial_summary: Per-trial annotation summary.
        join_audit: Full outer raw/offset join audit.
        offsets: Old-compatible event-offset table.
        raw_dir: Raw EDF input directory.
        offsets_path: Old-compatible offset CSV path.
        started_at: Timestamp captured at run start.

    Returns:
        Markdown summary text.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    zero_files = file_summary[file_summary["annotation_count"].eq(0)]
    low_files = file_summary[
        file_summary["annotation_count"].gt(0) & file_summary["annotation_count"].lt(LOW_ANNOTATION_COUNT_THRESHOLD)
    ]
    split_files = file_summary[file_summary["file_role"].eq("split_part")]
    concatenated_files = file_summary[file_summary["file_role"].eq("concatenated")]
    watchlist_files = file_summary[file_summary["watchlist_id"]]

    trial_groups = trial_summary[trial_summary["summary_row_kind"].eq("trial_group")]
    practice_groups = trial_groups[trial_groups["practice_like_group"]]
    bad_start_groups = trial_groups[trial_groups["bad_start_like_group"]]
    join_counts = join_audit["join_status"].value_counts(dropna=False).sort_index()
    file_role_counts = file_summary["file_role"].value_counts(dropna=False).sort_index()
    read_status_counts = file_summary["read_status"].value_counts(dropna=False).sort_index()
    join_issue_counts = count_issue_codes(join_audit["join_audit_issue_codes"])
    trial_issue_counts = count_issue_codes(trial_summary["trial_issue_codes"])

    dominant_code_missing_files = file_summary[
        file_summary["file_issue_codes"].fillna("").str.contains("missing_dominant_code", regex=False)
    ]

    lines = [
        "# Raw Annotation Alignment Summary",
        "",
        f"Generated: {finished_at.isoformat(timespec='seconds')}",
        f"Started: {started_at.isoformat(timespec='seconds')}",
        f"Duration seconds: {(finished_at - started_at).total_seconds():.1f}",
        "",
        "## Inputs",
        "",
        f"- Raw EDF directory: `{raw_dir.as_posix()}`",
        f"- Old-compatible event offsets: `{offsets_path.as_posix()}`",
        "- EDF read mode: `mne.io.read_raw_edf(..., preload=False)`",
        "- Trigger-name mapping source: `external/DEMI_EEG_Pipeline/edf2bids.py`",
        "",
        "## Files inspected",
        "",
        f"- EDF files inspected: {len(file_summary)}",
        f"- Parsed participants represented by EDFs: {file_summary['participant_id'].dropna().nunique()}",
        f"- Raw annotation rows written: {len(raw_events)}",
        f"- Old-compatible offset rows read: {len(offsets)}",
        f"- Old-compatible offset participants: {offsets['id'].dropna().nunique()}",
        "",
        "MNE read status:",
        "",
        *format_count_table(read_status_counts),
        "",
        "File roles:",
        "",
        *format_count_table(file_role_counts),
        "",
        "## Annotation-count flags",
        "",
        f"- Zero-annotation EDFs: {len(zero_files)} ({format_id_list(zero_files['participant_id'])})",
        f"- Low nonzero annotation EDFs (< {LOW_ANNOTATION_COUNT_THRESHOLD}): {len(low_files)} "
        f"({format_id_list(low_files['participant_id'])})",
        f"- EDFs missing one or more dominant task codes ({', '.join(str(code) for code in DOMINANT_TASK_EVENT_CODES)}): "
        f"{len(dominant_code_missing_files)} ({format_id_list(dominant_code_missing_files['participant_id'])})",
        "",
        "## Trial reconstruction",
        "",
        f"- Raw trial-like annotation groups: {len(trial_groups)}",
        f"- Practice-like groups by old duration rule: {len(practice_groups)}",
        f"- Bad-start-like groups by old trigger-order rule: {len(bad_start_groups)}",
        f"- Trial groups with old-compatible join trial count: {trial_groups['join_trial_count'].notna().sum()}",
        "",
        "Trial issue counts:",
        "",
    ]

    if trial_issue_counts:
        for issue, count in sorted(trial_issue_counts.items()):
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Offset join audit",
            "",
            "Join status counts:",
            "",
            *format_count_table(join_counts),
            "",
            f"- Trace epoch duration mismatches (> {TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD * 1000:.0f} ms): "
            f"{int(join_audit['old_trace_epoch_duration_mismatch'].fillna(False).sum())}",
            f"- Stimulus duration mismatches (> {STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD * 1000:.0f} ms): "
            f"{int(join_audit['old_stimulus_duration_mismatch'].fillna(False).sum())}",
            "",
            "Join issue counts:",
            "",
        ]
    )

    if join_issue_counts:
        for issue, count in sorted(join_issue_counts.items()):
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Watchlist and special file facts",
            "",
            f"- Watchlist IDs present in EDF files: {format_id_list(watchlist_files['participant_id'])}",
            f"- Split EDF rows: {len(split_files)} ({format_id_list(split_files['participant_id'])})",
            f"- Concatenated EDF rows: {len(concatenated_files)} ({format_id_list(concatenated_files['participant_id'])})",
            "",
        ]
    )

    if not split_files.empty:
        lines.append("Split EDF files:")
        lines.append("")
        for _, row in split_files.sort_values(["participant_id", "source_filename"]).iterrows():
            lines.append(
                f"- ID {int(row['participant_id'])}, `{row['source_filename']}`: "
                f"{int(row['annotation_count'])} annotations"
            )
        lines.append("")

    if not concatenated_files.empty:
        lines.append("Concatenated EDF files:")
        lines.append("")
        for _, row in concatenated_files.sort_values(["participant_id", "source_filename"]).iterrows():
            lines.append(
                f"- ID {int(row['participant_id'])}, `{row['source_filename']}`: "
                f"{int(row['annotation_count'])} annotations"
            )
        lines.append("")

    lines.extend(
        [
            "## Deviations and TODOs from old update_events parity",
            "",
            "- Raw EDF descriptions are numeric, so this script adds a derived semantic label layer before applying old trial logic.",
            "- Non-task annotations remain in `raw_annotation_events.csv` but do not affect old-style trial numbering.",
            "- Split EDFs are inspected as separate files; no event-continuity repair or concatenation is attempted here.",
            "- The concatenated ID 5 EDF is inspected as found; no duplicate-recording trim is attempted here.",
            "- Historical preprocessing-side duplicate-trigger relabeling is not applied in this read-only script.",
            "- `real_trace_start` and `real_trace_end` are not written back to EEG; only old-style mismatch quantities are audited.",
            "- Participant-retention and preprocessing decisions remain pending review.",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(
    raw_events: pd.DataFrame,
    trial_summary: pd.DataFrame,
    join_audit: pd.DataFrame,
    summary_text: str,
    output_dir: Path,
) -> None:
    """Write local event-alignment audit outputs.

    Args:
        raw_events: One-row-per-annotation table.
        trial_summary: Per-trial annotation summary.
        join_audit: Full outer raw/offset join audit.
        summary_text: Markdown summary text.
        output_dir: Output directory below ``_Data/eeg/``.

    Returns:
        None.

    Side effects:
        Creates ``output_dir`` if needed and writes four local-only audit files.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    raw_events.to_csv(output_dir / RAW_ANNOTATION_EVENTS_FILENAME, index=False)
    trial_summary.to_csv(output_dir / ANNOTATION_TRIAL_SUMMARY_FILENAME, index=False)
    join_audit.to_csv(output_dir / OFFSET_ANNOTATION_JOIN_AUDIT_FILENAME, index=False)
    (output_dir / ALIGNMENT_SUMMARY_FILENAME).write_text(summary_text, encoding="utf-8")


def run_alignment_audit(repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the raw annotation versus event-offset audit in memory.

    Args:
        repo_root: Repository root.

    Returns:
        ``(raw_events, file_summary, trial_summary, join_audit)``.

    Side effects:
        Reads raw EDF annotations and the old-compatible offset CSV. Does not
        write files.
    """

    raw_dir = repo_root / RAW_EEG_DIR
    offsets_path = repo_root / OLD_COMPATIBLE_OFFSETS_PATH
    edf_files = discover_edf_files(raw_dir)
    raw_events, file_summary = read_all_raw_annotations(edf_files, repo_root)
    raw_events = assign_old_update_trial_numbers(raw_events)
    trial_summary = build_annotation_trial_summary(raw_events, file_summary)
    offsets = read_old_compatible_offsets(offsets_path)
    join_audit = build_offset_annotation_join_audit(trial_summary, offsets, file_summary)
    return raw_events, file_summary, trial_summary, join_audit


def main() -> None:
    """Run the raw annotation comparison workflow.

    Args:
        None. Paths are resolved from this script's repository location.

    Returns:
        None.

    Side effects:
        Opens raw EDF files through MNE with ``preload=False``, reads the
        old-compatible event-offset CSV, and writes local audit outputs under
        ``_Data/eeg/event_alignment/``.
    """

    started_at = datetime.now()
    repo_root = repo_root_from_script()
    raw_dir = repo_root / RAW_EEG_DIR
    offsets_path = repo_root / OLD_COMPATIBLE_OFFSETS_PATH
    output_dir = repo_root / OUTPUT_DIR

    edf_files = discover_edf_files(raw_dir)
    raw_events, file_summary = read_all_raw_annotations(edf_files, repo_root)
    raw_events = assign_old_update_trial_numbers(raw_events)
    trial_summary = build_annotation_trial_summary(raw_events, file_summary)
    offsets = read_old_compatible_offsets(offsets_path)
    join_audit = build_offset_annotation_join_audit(trial_summary, offsets, file_summary)
    summary_text = build_alignment_summary(
        raw_events=raw_events,
        file_summary=file_summary,
        trial_summary=trial_summary,
        join_audit=join_audit,
        offsets=offsets,
        raw_dir=raw_dir,
        offsets_path=offsets_path,
        started_at=started_at,
    )
    write_outputs(raw_events, trial_summary, join_audit, summary_text, output_dir)

    print(f"Inspected {len(file_summary)} EDF files from {raw_dir}")
    print(f"Wrote event-alignment audit outputs to {output_dir}")


if __name__ == "__main__":
    main()
