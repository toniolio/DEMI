"""Audit focused DEMI EEG event-alignment special cases.

This script is a read-only follow-up to
``analysis/eeg_mne/02_audit_raw_event_sequences.py``. It does not build a new
event surface. Instead, it treats the proposed labels, trial numbers, join
statuses, and nearby timing candidates from script 02 as the current best
inspection surface and asks a narrower set of questions about the remaining
event-alignment blockers.

The focused audit is needed because the second raw event-sequence audit greatly
reduced the main timing mismatch surface, but it left several issues that still
need factual review before event-driven epoch construction can be designed:

- split-file continuity for IDs 54, 56, and 65;
- the concatenated/duplicate-file structure for ID 5;
- zero-annotation or no-offset-target cases 11, 14, 89, 94, and 100;
- classification of raw-only rows outside those known special cases;
- the clean-timing count mismatch between the Markdown summary and the
  ``clean_timing_row`` flag in the join CSV.

Inputs read from local paths:

- ``_Data/eeg/event_sequence_audit/proposed_annotation_events.csv``;
- ``_Data/eeg/event_sequence_audit/proposed_trial_summary.csv``;
- ``_Data/eeg/event_sequence_audit/proposed_offset_join_audit.csv``;
- ``_Data/eeg/event_sequence_audit/first_mismatch_by_file.csv``;
- ``_Data/eeg/event_sequence_audit/nearby_trial_candidate_audit.csv``;
- ``_Data/eeg/event_sequence_audit/event_sequence_audit_summary.md``;
- ``_Data/behavior/event_offsets/event_offsets_old_compatible.csv``.

Outputs written below the local-only directory
``_Data/eeg/event_special_case_audit/``:

- ``split_file_continuity_audit.csv``;
- ``id5_concatenated_file_audit.csv``;
- ``zero_no_offset_case_audit.csv``;
- ``raw_only_row_classification.csv``;
- ``clean_timing_definition_check.csv``;
- ``event_special_case_audit_summary.md``.

What this script explicitly does not do:

- It does not open, repair, rewrite, or save raw EEG files.
- It does not concatenate split EDFs.
- It does not crop, split, or otherwise alter the ID 5 EDF.
- It does not load EEG signal data.
- It does not preprocess, filter, epoch, interpolate, run ICA, or compute
  time-frequency features.
- It does not make participant-level or event-handling decisions. The output is
  an inspection aid for Tony's review.

Safety boundaries:

- All generated files are local outputs under ``_Data/``.
- The script uses ``pathlib`` and ``pandas`` and does not import MNE.
- Existing audit CSVs and the old-compatible offset table are read only.
- Error messages name missing files or missing columns before any partial audit
  output is written.
- Categories in the raw-only classifier are conservative labels for observable
  row patterns, not analytic decisions.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/03_audit_event_alignment_special_cases.py

The script also works from another working directory because paths are resolved
relative to this file's location in the repository.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


EVENT_SEQUENCE_AUDIT_DIR = Path("_Data") / "eeg" / "event_sequence_audit"
EVENT_SPECIAL_CASE_AUDIT_DIR = Path("_Data") / "eeg" / "event_special_case_audit"
OLD_COMPATIBLE_OFFSETS_PATH = Path("_Data") / "behavior" / "event_offsets" / "event_offsets_old_compatible.csv"

PROPOSED_ANNOTATION_EVENTS_FILENAME = "proposed_annotation_events.csv"
PROPOSED_TRIAL_SUMMARY_FILENAME = "proposed_trial_summary.csv"
PROPOSED_OFFSET_JOIN_AUDIT_FILENAME = "proposed_offset_join_audit.csv"
FIRST_MISMATCH_BY_FILE_FILENAME = "first_mismatch_by_file.csv"
NEARBY_TRIAL_CANDIDATE_AUDIT_FILENAME = "nearby_trial_candidate_audit.csv"
EVENT_SEQUENCE_AUDIT_SUMMARY_FILENAME = "event_sequence_audit_summary.md"

SPLIT_FILE_IDS = {54, 56, 65}
CONCATENATED_FILE_ID = 5
ZERO_NO_OFFSET_CASE_IDS = {11, 14, 89, 94, 100}
KNOWN_CONCATENATED_OR_SPLIT_IDS = SPLIT_FILE_IDS | {CONCATENATED_FILE_ID}

RAW_AND_OFFSET_STATUS = "raw_annotation_and_offset"
RAW_ONLY_STATUS = "raw_annotation_without_offset"
OFFSET_ONLY_STATUS = "offset_without_raw_annotation_trial"

TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD = 0.150
STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD = 0.250

SPLIT_FILE_CONTINUITY_AUDIT_FILENAME = "split_file_continuity_audit.csv"
ID5_CONCATENATED_FILE_AUDIT_FILENAME = "id5_concatenated_file_audit.csv"
ZERO_NO_OFFSET_CASE_AUDIT_FILENAME = "zero_no_offset_case_audit.csv"
RAW_ONLY_ROW_CLASSIFICATION_FILENAME = "raw_only_row_classification.csv"
CLEAN_TIMING_DEFINITION_CHECK_FILENAME = "clean_timing_definition_check.csv"
EVENT_SPECIAL_CASE_AUDIT_SUMMARY_FILENAME = "event_special_case_audit_summary.md"

INTEGER_COLUMNS = (
    "participant_id",
    "audit_participant_id",
    "audit_trial_count",
    "offset_id",
    "offset_trial_count",
    "trial_count",
    "split_part",
    "event_row_order",
    "raw_trial_sequence",
    "join_trial_count",
    "proposed_raw_trial_sequence",
    "proposed_join_trial_count",
    "candidate_offset_trial_count",
    "candidate_offset_delta",
    "candidate_timing_rank",
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script path.

    Args:
        None.

    Returns:
        Absolute repository root path.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def text_or_empty(value: Any) -> str:
    """Return stripped text while preserving missing values as empty text.

    Args:
        value: Scalar value from a CSV cell or derived table.

    Returns:
        Empty string for missing values, otherwise stripped text.

    Side effects:
        None.
    """

    if pd.isna(value):
        return ""
    return str(value).strip()


def is_true(value: Any) -> bool:
    """Interpret bool-like values read from CSVs.

    Args:
        value: Scalar value that may be boolean, numeric, text, or missing.

    Returns:
        ``True`` for true-like values and ``False`` otherwise.

    Side effects:
        None.
    """

    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def true_mask(frame: pd.DataFrame, column: str) -> pd.Series:
    """Return a boolean Series for a bool-like column.

    Args:
        frame: DataFrame to inspect.
        column: Column name to convert.

    Returns:
        Boolean Series aligned to ``frame.index``. Missing columns yield all
        ``False`` values.

    Side effects:
        None.
    """

    if column not in frame.columns:
        return pd.Series(False, index=frame.index)
    return frame[column].map(is_true)


def coerce_integer_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Coerce known integer-like columns to pandas nullable integers.

    Args:
        frame: DataFrame with audit columns that may be read as floats.

    Returns:
        Copy with known integer-like columns converted where present.

    Side effects:
        None.
    """

    out = frame.copy()
    for column in INTEGER_COLUMNS:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
    return out


def read_csv_checked(path: Path, required_columns: set[str], label: str) -> pd.DataFrame:
    """Read a CSV file and verify that required columns are present.

    Args:
        path: CSV path to read.
        required_columns: Column names required by this audit.
        label: Human-readable input label for failure messages.

    Returns:
        DataFrame read from ``path``.

    Side effects:
        Reads one CSV file. Raises ``RuntimeError`` when the file or required
        columns are missing.
    """

    if not path.exists():
        raise RuntimeError(f"required {label} is missing: {path.as_posix()}")

    frame = pd.read_csv(path, low_memory=False)
    missing = sorted(required_columns.difference(frame.columns))
    if missing:
        raise RuntimeError(f"{label} is missing required column(s): {', '.join(missing)}")
    return coerce_integer_columns(frame)


def read_text_checked(path: Path, label: str) -> str:
    """Read a UTF-8 text file with a clear missing-file error.

    Args:
        path: Text file path.
        label: Human-readable input label for failure messages.

    Returns:
        File text.

    Side effects:
        Reads one text file. Raises ``RuntimeError`` if the file is missing.
    """

    if not path.exists():
        raise RuntimeError(f"required {label} is missing: {path.as_posix()}")
    return path.read_text(encoding="utf-8")


def read_audit_inputs(repo_root: Path) -> dict[str, Any]:
    """Read all event-alignment special-case audit inputs.

    Args:
        repo_root: Absolute repository root path.

    Returns:
        Dictionary containing proposed events, proposed trial summary, proposed
        join audit, first mismatch table, nearby-candidate table, Markdown
        summary text, and old-compatible offsets.

    Side effects:
        Reads local CSV and Markdown files below ``_Data/``.
    """

    sequence_dir = repo_root / EVENT_SEQUENCE_AUDIT_DIR
    events = read_csv_checked(
        sequence_dir / PROPOSED_ANNOTATION_EVENTS_FILENAME,
        {
            "participant_id",
            "source_filename",
            "file_role",
            "split_part",
            "event_row_order",
            "onset_seconds",
            "annotation_description",
            "proposed_event_name",
            "proposed_sequence_action",
            "proposed_mapped_task_event",
            "proposed_raw_trial_sequence",
            "proposed_join_trial_count",
        },
        "proposed annotation-event CSV from script 02",
    )
    trials = read_csv_checked(
        sequence_dir / PROPOSED_TRIAL_SUMMARY_FILENAME,
        {
            "summary_row_kind",
            "participant_id",
            "source_filename",
            "file_role",
            "split_part",
            "annotation_count",
            "mapped_task_event_count",
            "proposed_raw_trial_sequence",
            "proposed_join_trial_count",
            "practice_like_group",
            "bad_start_like_group",
            "trial_event_count",
            "missing_expected_event_names",
            "extra_expected_event_names",
            "trial_issue_codes",
        },
        "proposed trial-summary CSV from script 02",
    )
    join = read_csv_checked(
        sequence_dir / PROPOSED_OFFSET_JOIN_AUDIT_FILENAME,
        {
            "audit_participant_id",
            "audit_trial_count",
            "join_status",
            "clean_timing_row",
            "source_filename",
            "file_role",
            "split_part",
            "raw_trial_sequence",
            "join_trial_count",
            "raw_join_key_duplicate",
            "old_trace_epoch_duration_mismatch",
            "old_stimulus_duration_mismatch",
            "missing_expected_event_names",
            "extra_expected_event_names",
            "alignment_problem_codes",
            "trial_issue_codes",
        },
        "proposed offset-join audit CSV from script 02",
    )
    first_mismatch = read_csv_checked(
        sequence_dir / FIRST_MISMATCH_BY_FILE_FILENAME,
        {
            "participant_id",
            "source_filename",
            "file_role",
            "split_part",
            "has_alignment_problem",
            "first_problem_codes",
            "file_join_rows",
            "file_problem_rows",
            "possible_cascade_after_first_problem",
        },
        "first-mismatch CSV from script 02",
    )
    nearby = read_csv_checked(
        sequence_dir / NEARBY_TRIAL_CANDIDATE_AUDIT_FILENAME,
        {
            "problem_row_number",
            "audit_participant_id",
            "source_filename",
            "file_role",
            "split_part",
            "raw_trial_sequence",
            "join_trial_count",
            "candidate_offset_trial_count",
            "candidate_offset_delta",
            "candidate_timing_score_seconds",
            "candidate_trace_within_old_threshold",
            "candidate_stimulus_within_old_threshold",
            "candidate_clean_by_timing",
            "nominal_candidate_timing_score_seconds",
            "candidate_better_than_nominal",
            "candidate_timing_rank",
        },
        "nearby trial-candidate CSV from script 02",
    )
    offsets = read_csv_checked(
        repo_root / OLD_COMPATIBLE_OFFSETS_PATH,
        {"id", "trial_count", "stimulus_mt", "end_trigger"},
        "old-compatible event-offset CSV",
    )
    summary_text = read_text_checked(
        sequence_dir / EVENT_SEQUENCE_AUDIT_SUMMARY_FILENAME,
        "event-sequence audit Markdown summary from script 02",
    )

    return {
        "events": events,
        "trials": trials,
        "join": join,
        "first_mismatch": first_mismatch,
        "nearby": nearby,
        "offsets": offsets,
        "summary_text": summary_text,
    }


def numeric_range_text(values: pd.Series, digits: int = 0) -> str:
    """Format the numeric range of a Series as compact text.

    Args:
        values: Series containing numeric or missing values.
        digits: Number of decimal places for the output.

    Returns:
        ``"min-max"`` text, or an empty string when no numeric value exists.

    Side effects:
        None.
    """

    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return ""
    if digits == 0:
        return f"{int(numeric.min())}-{int(numeric.max())}"
    return f"{numeric.min():.{digits}f}-{numeric.max():.{digits}f}"


def joined_unique_text(values: pd.Series) -> str:
    """Join unique non-empty text values in first-seen order.

    Args:
        values: Series containing text-like values.

    Returns:
        Semicolon-delimited text.

    Side effects:
        None.
    """

    seen: dict[str, None] = {}
    for value in values:
        text = text_or_empty(value)
        if text:
            seen[text] = None
    return ";".join(seen.keys())


def format_int_values(values: pd.Series) -> str:
    """Format unique integer values as comma-delimited text.

    Args:
        values: Series containing integer-like values.

    Returns:
        Sorted comma-delimited integer text, or an empty string.

    Side effects:
        None.
    """

    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return ""
    return ",".join(str(int(value)) for value in sorted(numeric.unique()))


def value_counts_text(values: pd.Series) -> str:
    """Format value counts as semicolon-delimited ``value=count`` text.

    Args:
        values: Series to count after missing values are removed.

    Returns:
        Semicolon-delimited count text.

    Side effects:
        None.
    """

    cleaned = values.dropna().map(text_or_empty)
    cleaned = cleaned[cleaned.ne("")]
    if cleaned.empty:
        return ""
    counts = cleaned.value_counts(sort=False)
    return ";".join(f"{label}={int(count)}" for label, count in counts.items())


def split_issue_codes(value: Any) -> list[str]:
    """Split a semicolon-delimited issue-code cell.

    Args:
        value: Issue-code cell value.

    Returns:
        List of non-empty issue codes.

    Side effects:
        None.
    """

    return [code for code in text_or_empty(value).split(";") if code]


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
    for value in values:
        for code in split_issue_codes(value):
            counts[code] += 1
    return counts


def strict_clean_timing_mask(join: pd.DataFrame) -> pd.Series:
    """Recompute the stricter clean-timing definition for join rows.

    Args:
        join: Proposed offset-join audit table or subset.

    Returns:
        Boolean Series that is true only for matched raw/offset rows with no
        trace-duration mismatch, no stimulus-duration mismatch, no missing
        expected event, no extra expected event, and no duplicate raw join key.

    Side effects:
        None.
    """

    return (
        join["join_status"].eq(RAW_AND_OFFSET_STATUS)
        & ~true_mask(join, "old_trace_epoch_duration_mismatch")
        & ~true_mask(join, "old_stimulus_duration_mismatch")
        & join["missing_expected_event_names"].fillna("").eq("")
        & join["extra_expected_event_names"].fillna("").eq("")
        & ~true_mask(join, "raw_join_key_duplicate")
    )


def best_nearby_candidate_rows(nearby: pd.DataFrame) -> pd.DataFrame:
    """Select the best nearby offset candidate for each problem row.

    Args:
        nearby: Nearby trial-candidate audit from script 02.

    Returns:
        One row per problem row/key, carrying the lowest-ranked candidate and
        renamed ``best_candidate_*`` columns for easier merging.

    Side effects:
        None.
    """

    if nearby.empty:
        return pd.DataFrame(
            columns=[
                "audit_participant_id",
                "source_filename",
                "raw_trial_sequence",
                "join_trial_count",
                "best_candidate_offset_trial_count",
                "best_candidate_offset_delta",
                "best_candidate_timing_score_seconds",
                "best_candidate_trace_within_old_threshold",
                "best_candidate_stimulus_within_old_threshold",
                "best_candidate_clean_by_timing",
                "best_candidate_better_than_nominal",
            ]
        )

    sort_columns = [
        "audit_participant_id",
        "source_filename",
        "raw_trial_sequence",
        "join_trial_count",
        "problem_row_number",
        "candidate_timing_rank",
    ]
    ordered = nearby.sort_values(sort_columns, kind="mergesort", na_position="last")
    best = (
        ordered.groupby(
            ["audit_participant_id", "source_filename", "raw_trial_sequence", "join_trial_count"],
            dropna=False,
            as_index=False,
        )
        .first()
        .copy()
    )
    return best.rename(
        columns={
            "candidate_offset_trial_count": "best_candidate_offset_trial_count",
            "candidate_offset_delta": "best_candidate_offset_delta",
            "candidate_timing_score_seconds": "best_candidate_timing_score_seconds",
            "candidate_trace_within_old_threshold": "best_candidate_trace_within_old_threshold",
            "candidate_stimulus_within_old_threshold": "best_candidate_stimulus_within_old_threshold",
            "candidate_clean_by_timing": "best_candidate_clean_by_timing",
            "candidate_better_than_nominal": "best_candidate_better_than_nominal",
        }
    )


def summarize_best_nearby(best_rows: pd.DataFrame) -> dict[str, Any]:
    """Summarize best nearby-candidate diagnostics for a row group.

    Args:
        best_rows: One-best-candidate rows for a participant/file subset.

    Returns:
        Dictionary with problem-row counts and timing-candidate summaries.

    Side effects:
        None.
    """

    if best_rows.empty:
        return {
            "nearby_problem_rows": 0,
            "nearby_best_better_than_nominal_rows": 0,
            "nearby_best_clean_by_timing_rows": 0,
            "nearby_best_delta_values": "",
            "nearby_best_delta_mode": "",
            "nearby_best_timing_score_min_seconds": pd.NA,
            "nearby_best_timing_score_median_seconds": pd.NA,
        }

    best_delta = pd.to_numeric(best_rows["best_candidate_offset_delta"], errors="coerce").dropna()
    best_score = pd.to_numeric(best_rows["best_candidate_timing_score_seconds"], errors="coerce").dropna()
    if best_delta.empty:
        delta_mode = ""
        delta_values = ""
    else:
        delta_mode = str(int(best_delta.mode().iloc[0]))
        delta_values = ",".join(str(int(value)) for value in sorted(best_delta.unique()))

    return {
        "nearby_problem_rows": int(len(best_rows)),
        "nearby_best_better_than_nominal_rows": int(true_mask(best_rows, "best_candidate_better_than_nominal").sum()),
        "nearby_best_clean_by_timing_rows": int(true_mask(best_rows, "best_candidate_clean_by_timing").sum()),
        "nearby_best_delta_values": delta_values,
        "nearby_best_delta_mode": delta_mode,
        "nearby_best_timing_score_min_seconds": best_score.min() if not best_score.empty else pd.NA,
        "nearby_best_timing_score_median_seconds": best_score.median() if not best_score.empty else pd.NA,
    }


def continuity_observation(row: dict[str, Any]) -> str:
    """Describe whether split-file continuity is visible in audit facts.

    Args:
        row: Split-file summary row with join and nearby-candidate counts.

    Returns:
        Short factual observation for the ``continuous_alignment_observable``
        column.

    Side effects:
        None.
    """

    if row["nearby_best_better_than_nominal_rows"] > 0:
        return "candidate-level continuity signal; nominal counts still reset_or_duplicate"
    if row["nearby_best_clean_by_timing_rows"] > 0 and row["duplicate_raw_join_key_rows"] > 0:
        return "within-part timing signal; duplicate keys still unresolved"
    if row["duplicate_raw_join_key_rows"] > 0 or row["offset_only_rows"] > 0 or row["raw_only_rows"] > 0:
        return "not_resolved_on_nominal_join_surface"
    return "no_continuity_problem_visible_in_current_rows"


def aggregate_join_rows(join_rows: pd.DataFrame) -> dict[str, Any]:
    """Aggregate join-status, duplicate-key, and timing facts.

    Args:
        join_rows: Join-audit subset for one participant/file/segment.

    Returns:
        Dictionary of row counts and trial-count ranges.

    Side effects:
        None.
    """

    if join_rows.empty:
        return {
            "join_rows": 0,
            "matched_rows": 0,
            "raw_only_rows": 0,
            "offset_only_rows": 0,
            "duplicate_raw_join_key_rows": 0,
            "duplicate_join_trial_counts": "",
            "raw_only_trial_counts": "",
            "offset_only_trial_counts": "",
            "csv_clean_timing_rows": 0,
            "strict_clean_timing_rows": 0,
            "trace_mismatch_rows": 0,
            "stimulus_mismatch_rows": 0,
            "missing_expected_rows": 0,
            "extra_expected_rows": 0,
        }

    duplicate_mask = true_mask(join_rows, "raw_join_key_duplicate")
    raw_only_mask = join_rows["join_status"].eq(RAW_ONLY_STATUS)
    offset_only_mask = join_rows["join_status"].eq(OFFSET_ONLY_STATUS)
    missing_expected = join_rows["missing_expected_event_names"].fillna("").ne("")
    extra_expected = join_rows["extra_expected_event_names"].fillna("").ne("")

    return {
        "join_rows": int(len(join_rows)),
        "matched_rows": int(join_rows["join_status"].eq(RAW_AND_OFFSET_STATUS).sum()),
        "raw_only_rows": int(raw_only_mask.sum()),
        "offset_only_rows": int(offset_only_mask.sum()),
        "duplicate_raw_join_key_rows": int(duplicate_mask.sum()),
        "duplicate_join_trial_counts": format_int_values(join_rows.loc[duplicate_mask, "join_trial_count"]),
        "raw_only_trial_counts": format_int_values(join_rows.loc[raw_only_mask, "audit_trial_count"]),
        "offset_only_trial_counts": format_int_values(join_rows.loc[offset_only_mask, "audit_trial_count"]),
        "csv_clean_timing_rows": int(true_mask(join_rows, "clean_timing_row").sum()),
        "strict_clean_timing_rows": int(strict_clean_timing_mask(join_rows).sum()),
        "trace_mismatch_rows": int(true_mask(join_rows, "old_trace_epoch_duration_mismatch").sum()),
        "stimulus_mismatch_rows": int(true_mask(join_rows, "old_stimulus_duration_mismatch").sum()),
        "missing_expected_rows": int(missing_expected.sum()),
        "extra_expected_rows": int(extra_expected.sum()),
    }


def build_split_file_continuity_audit(
    events: pd.DataFrame,
    trials: pd.DataFrame,
    join: pd.DataFrame,
    best_nearby: pd.DataFrame,
) -> pd.DataFrame:
    """Build the split-file continuity audit for IDs 54, 56, and 65.

    Args:
        events: Proposed annotation-event table from script 02.
        trials: Proposed trial-summary table from script 02.
        join: Proposed offset-join audit from script 02.
        best_nearby: Best nearby-candidate rows derived from script 02.

    Returns:
        DataFrame with one row per split EDF part plus one participant-total
        row per split ID.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []

    for participant_id in sorted(SPLIT_FILE_IDS):
        participant_events = events[events["participant_id"].eq(participant_id)].copy()
        participant_trials = trials[trials["participant_id"].eq(participant_id)].copy()
        participant_join = join[join["audit_participant_id"].eq(participant_id)].copy()
        participant_best = best_nearby[best_nearby["audit_participant_id"].eq(participant_id)].copy()

        for source_filename, file_events in participant_events.groupby("source_filename", sort=True):
            file_trials = participant_trials[participant_trials["source_filename"].eq(source_filename)].copy()
            file_join = participant_join[participant_join["source_filename"].eq(source_filename)].copy()
            file_best = participant_best[participant_best["source_filename"].eq(source_filename)].copy()
            task_events = file_events[true_mask(file_events, "proposed_mapped_task_event")]
            join_counts = aggregate_join_rows(file_join)
            nearby_counts = summarize_best_nearby(file_best)

            row = {
                "participant_id": participant_id,
                "row_kind": "split_part",
                "source_filename": source_filename,
                "file_role": joined_unique_text(file_events["file_role"]),
                "split_part": file_events["split_part"].dropna().iloc[0] if file_events["split_part"].notna().any() else pd.NA,
                "event_rows": int(len(file_events)),
                "proposed_task_event_rows": int(len(task_events)),
                "event_row_order_range": numeric_range_text(file_events["event_row_order"]),
                "onset_seconds_range": numeric_range_text(file_events["onset_seconds"], digits=3),
                "proposed_raw_trial_range_from_events": numeric_range_text(task_events["proposed_raw_trial_sequence"]),
                "proposed_join_trial_range_from_events": numeric_range_text(task_events["proposed_join_trial_count"]),
                "trial_group_rows": int(file_trials["summary_row_kind"].eq("trial_group").sum()),
                "proposed_raw_trial_range_from_trials": numeric_range_text(file_trials["proposed_raw_trial_sequence"]),
                "proposed_join_trial_range_from_trials": numeric_range_text(file_trials["proposed_join_trial_count"]),
                "proposed_sequence_action_counts": value_counts_text(file_events["proposed_sequence_action"]),
                **join_counts,
                **nearby_counts,
            }
            row["continuous_alignment_observable"] = continuity_observation(row)
            row["factual_note"] = (
                "split EDF part audited as a separate file; no concatenation or repair performed"
            )
            rows.append(row)

        total_join_counts = aggregate_join_rows(participant_join)
        total_nearby_counts = summarize_best_nearby(participant_best)
        total_row = {
            "participant_id": participant_id,
            "row_kind": "participant_total",
            "source_filename": "all_split_parts",
            "file_role": joined_unique_text(participant_events["file_role"]),
            "split_part": pd.NA,
            "event_rows": int(len(participant_events)),
            "proposed_task_event_rows": int(true_mask(participant_events, "proposed_mapped_task_event").sum()),
            "event_row_order_range": "",
            "onset_seconds_range": numeric_range_text(participant_events["onset_seconds"], digits=3),
            "proposed_raw_trial_range_from_events": numeric_range_text(
                participant_events["proposed_raw_trial_sequence"]
            ),
            "proposed_join_trial_range_from_events": numeric_range_text(
                participant_events["proposed_join_trial_count"]
            ),
            "trial_group_rows": int(participant_trials["summary_row_kind"].eq("trial_group").sum()),
            "proposed_raw_trial_range_from_trials": numeric_range_text(
                participant_trials["proposed_raw_trial_sequence"]
            ),
            "proposed_join_trial_range_from_trials": numeric_range_text(
                participant_trials["proposed_join_trial_count"]
            ),
            "proposed_sequence_action_counts": value_counts_text(participant_events["proposed_sequence_action"]),
            **total_join_counts,
            **total_nearby_counts,
        }
        total_row["continuous_alignment_observable"] = continuity_observation(total_row)
        total_row["factual_note"] = (
            "participant total across split parts; offset-only rows have no source EDF on the join surface"
        )
        rows.append(total_row)

    return pd.DataFrame(rows)


def id5_event_segment(onset_seconds: Any, file_start_onset: float | None) -> str:
    """Classify an ID 5 annotation or trial by position around ``file start``.

    Args:
        onset_seconds: Annotation or trial onset in seconds.
        file_start_onset: Onset of the ``file start`` marker, or ``None``.

    Returns:
        ``before_file_start``, ``file_start_row``, ``after_file_start``, or
        ``unknown_segment``.

    Side effects:
        None.
    """

    onset = pd.to_numeric(pd.Series([onset_seconds]), errors="coerce").iloc[0]
    if pd.isna(onset) or file_start_onset is None:
        return "unknown_segment"
    if onset < file_start_onset:
        return "before_file_start"
    if onset > file_start_onset:
        return "after_file_start"
    return "file_start_row"


def first_available_trial_onset(row: pd.Series) -> Any:
    """Return the earliest available trial-event onset from a trial row.

    Args:
        row: One proposed trial-summary row.

    Returns:
        First non-missing onset among expected task-event onset columns, or
        ``pd.NA`` if none is present.

    Side effects:
        None.
    """

    for column in [
        "stim_on_onset_seconds",
        "red_on_onset_seconds",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
    ]:
        if column in row.index and pd.notna(row[column]):
            return row[column]
    return pd.NA


def build_id5_concatenated_file_audit(
    events: pd.DataFrame,
    trials: pd.DataFrame,
    join: pd.DataFrame,
) -> pd.DataFrame:
    """Build the ID 5 concatenated/duplicate-file audit.

    Args:
        events: Proposed annotation-event table from script 02.
        trials: Proposed trial-summary table from script 02.
        join: Proposed offset-join audit from script 02.

    Returns:
        DataFrame summarizing ID 5 event ranges before/after ``file start``,
        trial-count structure, join statuses, and matched timing cleanliness.

    Side effects:
        None.
    """

    id5_events = events[events["participant_id"].eq(CONCATENATED_FILE_ID)].copy()
    id5_trials = trials[trials["participant_id"].eq(CONCATENATED_FILE_ID)].copy()
    id5_join = join[join["audit_participant_id"].eq(CONCATENATED_FILE_ID)].copy()

    file_start_rows = id5_events[id5_events["annotation_description"].map(text_or_empty).eq("file start")]
    file_start_onset = None
    if not file_start_rows.empty:
        file_start_onset = float(pd.to_numeric(file_start_rows["onset_seconds"], errors="coerce").dropna().iloc[0])

    id5_events["id5_segment"] = id5_events["onset_seconds"].map(
        lambda value: id5_event_segment(value, file_start_onset)
    )

    trial_groups = id5_trials[id5_trials["summary_row_kind"].eq("trial_group")].copy()
    if not trial_groups.empty:
        trial_groups["id5_trial_onset_seconds"] = trial_groups.apply(first_available_trial_onset, axis=1)
        trial_groups["id5_segment"] = trial_groups["id5_trial_onset_seconds"].map(
            lambda value: id5_event_segment(value, file_start_onset)
        )
    else:
        trial_groups["id5_trial_onset_seconds"] = pd.Series(dtype="float64")
        trial_groups["id5_segment"] = pd.Series(dtype="object")

    segment_key = trial_groups[
        ["source_filename", "raw_trial_sequence", "id5_segment"]
    ].rename(columns={"id5_segment": "join_segment"})
    id5_join = id5_join.merge(segment_key, on=["source_filename", "raw_trial_sequence"], how="left")
    id5_join["join_segment"] = id5_join["join_segment"].fillna("offset_only_or_unmapped_segment")

    rows: list[dict[str, Any]] = []
    for segment in ["before_file_start", "file_start_row", "after_file_start"]:
        segment_events = id5_events[id5_events["id5_segment"].eq(segment)].copy()
        segment_trials = trial_groups[trial_groups["id5_segment"].eq(segment)].copy()
        segment_join = id5_join[id5_join["join_segment"].eq(segment)].copy()
        task_events = segment_events[true_mask(segment_events, "proposed_mapped_task_event")]
        join_counts = aggregate_join_rows(segment_join)
        matched_rows = join_counts["matched_rows"]
        strict_clean_rows = join_counts["strict_clean_timing_rows"]

        rows.append(
            {
                "participant_id": CONCATENATED_FILE_ID,
                "row_kind": "file_start_segment",
                "segment": segment,
                "source_filename": joined_unique_text(id5_events["source_filename"]),
                "file_start_onset_seconds": file_start_onset if file_start_onset is not None else pd.NA,
                "event_rows": int(len(segment_events)),
                "proposed_task_event_rows": int(len(task_events)),
                "event_onset_seconds_range": numeric_range_text(segment_events["onset_seconds"], digits=3),
                "proposed_raw_trial_range_from_events": numeric_range_text(
                    task_events["proposed_raw_trial_sequence"]
                ),
                "proposed_join_trial_range_from_events": numeric_range_text(
                    task_events["proposed_join_trial_count"]
                ),
                "trial_group_rows": int(len(segment_trials)),
                "practice_like_trial_groups": int(true_mask(segment_trials, "practice_like_group").sum()),
                "bad_start_like_trial_groups": int(true_mask(segment_trials, "bad_start_like_group").sum()),
                "proposed_raw_trial_range_from_trials": numeric_range_text(
                    segment_trials["proposed_raw_trial_sequence"]
                ),
                "proposed_join_trial_range_from_trials": numeric_range_text(
                    segment_trials["proposed_join_trial_count"]
                ),
                "matched_rows_timing_clean_under_strict_definition": (
                    "yes" if matched_rows > 0 and matched_rows == strict_clean_rows else "no"
                ),
                **join_counts,
                "factual_note": "ID 5 segment audited around file start; no crop or split performed",
            }
        )

    total_join_counts = aggregate_join_rows(id5_join)
    rows.append(
        {
            "participant_id": CONCATENATED_FILE_ID,
            "row_kind": "participant_total",
            "segment": "all_segments",
            "source_filename": joined_unique_text(id5_events["source_filename"]),
            "file_start_onset_seconds": file_start_onset if file_start_onset is not None else pd.NA,
            "event_rows": int(len(id5_events)),
            "proposed_task_event_rows": int(true_mask(id5_events, "proposed_mapped_task_event").sum()),
            "event_onset_seconds_range": numeric_range_text(id5_events["onset_seconds"], digits=3),
            "proposed_raw_trial_range_from_events": numeric_range_text(id5_events["proposed_raw_trial_sequence"]),
            "proposed_join_trial_range_from_events": numeric_range_text(id5_events["proposed_join_trial_count"]),
            "trial_group_rows": int(len(trial_groups)),
            "practice_like_trial_groups": int(true_mask(trial_groups, "practice_like_group").sum()),
            "bad_start_like_trial_groups": int(true_mask(trial_groups, "bad_start_like_group").sum()),
            "proposed_raw_trial_range_from_trials": numeric_range_text(
                trial_groups["proposed_raw_trial_sequence"]
            ),
            "proposed_join_trial_range_from_trials": numeric_range_text(
                trial_groups["proposed_join_trial_count"]
            ),
            "matched_rows_timing_clean_under_strict_definition": (
                "yes"
                if total_join_counts["matched_rows"] > 0
                and total_join_counts["matched_rows"] == total_join_counts["strict_clean_timing_rows"]
                else "no"
            ),
            **total_join_counts,
            "factual_note": "ID 5 total across both sides of file start; no crop or split performed",
        }
    )

    return pd.DataFrame(rows)


def build_zero_no_offset_case_audit(
    events: pd.DataFrame,
    trials: pd.DataFrame,
    join: pd.DataFrame,
    offsets: pd.DataFrame,
) -> pd.DataFrame:
    """Build the zero-annotation and no-offset special-case audit.

    Args:
        events: Proposed annotation-event table from script 02.
        trials: Proposed trial-summary table from script 02.
        join: Proposed offset-join audit from script 02.
        offsets: Old-compatible event-offset table.

    Returns:
        DataFrame with one factual row for each of IDs 11, 14, 89, 94, and 100.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    for participant_id in sorted(ZERO_NO_OFFSET_CASE_IDS):
        participant_events = events[events["participant_id"].eq(participant_id)].copy()
        participant_trials = trials[trials["participant_id"].eq(participant_id)].copy()
        participant_join = join[join["audit_participant_id"].eq(participant_id)].copy()
        participant_offsets = offsets[offsets["id"].eq(participant_id)].copy()
        trial_groups = participant_trials[participant_trials["summary_row_kind"].eq("trial_group")]

        annotation_count_values = pd.to_numeric(participant_trials["annotation_count"], errors="coerce").dropna()
        annotation_count = int(annotation_count_values.max()) if not annotation_count_values.empty else int(len(participant_events))
        offset_rows = int(len(participant_offsets))

        if annotation_count == 0 and offset_rows > 0:
            note = "zero annotations in script-02 surface; old-compatible offsets exist"
        elif annotation_count == 0 and offset_rows == 0:
            note = "zero annotations in script-02 surface; no old-compatible offset rows"
        elif annotation_count > 0 and offset_rows == 0:
            note = "raw annotation sequence present; no old-compatible offset rows"
        else:
            note = "raw annotation sequence and old-compatible offsets both present"

        rows.append(
            {
                "participant_id": participant_id,
                "source_filenames": joined_unique_text(participant_trials["source_filename"])
                or joined_unique_text(participant_events["source_filename"]),
                "file_roles": joined_unique_text(participant_trials["file_role"])
                or joined_unique_text(participant_events["file_role"]),
                "annotation_count_from_trial_summary": annotation_count,
                "annotation_event_rows": int(len(participant_events)),
                "proposed_task_event_rows": int(true_mask(participant_events, "proposed_mapped_task_event").sum()),
                "annotation_description_counts": value_counts_text(participant_events["annotation_description"]),
                "proposed_event_name_counts": value_counts_text(participant_events["proposed_event_name"]),
                "trial_group_rows": int(len(trial_groups)),
                "proposed_raw_trial_range": numeric_range_text(trial_groups["proposed_raw_trial_sequence"]),
                "proposed_join_trial_range": numeric_range_text(trial_groups["proposed_join_trial_count"]),
                "join_status_counts": value_counts_text(participant_join["join_status"]),
                "matched_rows": int(participant_join["join_status"].eq(RAW_AND_OFFSET_STATUS).sum()),
                "raw_only_rows": int(participant_join["join_status"].eq(RAW_ONLY_STATUS).sum()),
                "offset_only_rows": int(participant_join["join_status"].eq(OFFSET_ONLY_STATUS).sum()),
                "old_compatible_offset_rows": offset_rows,
                "old_compatible_offset_trial_range": numeric_range_text(participant_offsets["trial_count"]),
                "factual_note": note,
            }
        )

    return pd.DataFrame(rows)


def offset_trial_maps(offsets: pd.DataFrame) -> tuple[dict[int, set[int]], dict[int, tuple[int, int]]]:
    """Create participant-level offset trial lookup maps.

    Args:
        offsets: Old-compatible event-offset table.

    Returns:
        ``(trial_sets, trial_ranges)`` where keys are participant IDs.

    Side effects:
        None.
    """

    trial_sets: dict[int, set[int]] = {}
    trial_ranges: dict[int, tuple[int, int]] = {}
    for participant_id, group in offsets.groupby("id", dropna=True):
        trials = pd.to_numeric(group["trial_count"], errors="coerce").dropna().astype(int)
        pid = int(participant_id)
        trial_sets[pid] = set(trials)
        if trials.empty:
            continue
        trial_ranges[pid] = (int(trials.min()), int(trials.max()))
    return trial_sets, trial_ranges


def offset_range_position(trial_count: Any, offset_range: tuple[int, int] | None) -> str:
    """Describe where a raw-only trial count sits relative to offset counts.

    Args:
        trial_count: Raw-derived audit trial count.
        offset_range: ``(min_trial, max_trial)`` for the participant, or
            ``None`` when no offsets exist.

    Returns:
        ``before_offset_range``, ``after_offset_range``,
        ``inside_offset_range``, or ``no_offset_range``.

    Side effects:
        None.
    """

    if offset_range is None:
        return "no_offset_range"
    numeric = pd.to_numeric(pd.Series([trial_count]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return "unknown_trial_count"
    if int(numeric) < offset_range[0]:
        return "before_offset_range"
    if int(numeric) > offset_range[1]:
        return "after_offset_range"
    return "inside_offset_range"


def raw_only_category(row: pd.Series) -> tuple[str, str]:
    """Classify one raw-only join row with conservative factual rules.

    Args:
        row: Raw-only join row after participant offset facts and best nearby
            candidate facts have been added.

    Returns:
        ``(category, reason)`` pair.

    Side effects:
        None.
    """

    participant_id = int(row["audit_participant_id"]) if pd.notna(row["audit_participant_id"]) else None
    file_role = text_or_empty(row.get("file_role", ""))
    missing_expected = text_or_empty(row.get("missing_expected_event_names", ""))
    extra_expected = text_or_empty(row.get("extra_expected_event_names", ""))

    if participant_id in KNOWN_CONCATENATED_OR_SPLIT_IDS or file_role in {"concatenated", "split_part"}:
        return (
            "known_concatenated_or_split_case",
            "row belongs to ID 5 or split-file IDs 54/56/65; handled in focused special-case tables",
        )

    if missing_expected or extra_expected:
        return (
            "malformed_or_non_task_annotation",
            "raw-only trial group has missing or extra expected task-event structure",
        )

    if not is_true(row.get("participant_has_old_compatible_offsets", False)):
        return (
            "participant_has_no_old_compatible_offsets",
            "participant has no old-compatible offset rows in this comparison surface",
        )

    # A raw-only row with a clean nearby nonzero candidate can reflect trial
    # count drift or another sequence issue. This is not a repair rule; it is a
    # flag that the nearby timing diagnostics found a plausible alternative
    # count for Tony to review.
    if is_true(row.get("best_candidate_clean_by_timing", False)):
        return (
            "possible_event_sequence_problem",
            "best nearby offset candidate is clean by both timing checks at a non-nominal trial count",
        )

    if not is_true(row.get("offset_trial_count_present", False)):
        return (
            "likely_old_compatible_row_surface_difference",
            "raw-derived trial count is absent from the old-compatible offset row surface",
        )

    return (
        "needs_manual_review",
        "raw-only row did not match the conservative classifier rules",
    )


def build_raw_only_row_classification(
    join: pd.DataFrame,
    offsets: pd.DataFrame,
    best_nearby: pd.DataFrame,
) -> pd.DataFrame:
    """Classify raw-only rows from the proposed join audit.

    Args:
        join: Proposed offset-join audit from script 02.
        offsets: Old-compatible event-offset table.
        best_nearby: Best nearby-candidate rows derived from script 02.

    Returns:
        One row per raw-only join row with a conservative category, reason,
        offset availability facts, and best nearby timing candidate fields.

    Side effects:
        None.
    """

    raw_only = join[join["join_status"].eq(RAW_ONLY_STATUS)].copy()
    if raw_only.empty:
        return raw_only

    merge_columns = [
        "audit_participant_id",
        "source_filename",
        "raw_trial_sequence",
        "join_trial_count",
        "best_candidate_offset_trial_count",
        "best_candidate_offset_delta",
        "best_candidate_timing_score_seconds",
        "best_candidate_trace_within_old_threshold",
        "best_candidate_stimulus_within_old_threshold",
        "best_candidate_clean_by_timing",
        "best_candidate_better_than_nominal",
    ]
    raw_only = raw_only.merge(
        best_nearby[[column for column in merge_columns if column in best_nearby.columns]],
        on=["audit_participant_id", "source_filename", "raw_trial_sequence", "join_trial_count"],
        how="left",
    )

    trial_sets, trial_ranges = offset_trial_maps(offsets)
    offset_row_counts = offsets.groupby("id").size().to_dict()

    raw_only["participant_has_old_compatible_offsets"] = raw_only["audit_participant_id"].map(
        lambda value: int(value) in trial_sets if pd.notna(value) else False
    )
    raw_only["participant_old_compatible_offset_rows"] = raw_only["audit_participant_id"].map(
        lambda value: int(offset_row_counts.get(int(value), 0)) if pd.notna(value) else 0
    )
    raw_only["participant_old_compatible_offset_trial_range"] = raw_only["audit_participant_id"].map(
        lambda value: (
            f"{trial_ranges[int(value)][0]}-{trial_ranges[int(value)][1]}"
            if pd.notna(value) and int(value) in trial_ranges
            else ""
        )
    )
    raw_only["offset_range_position"] = raw_only.apply(
        lambda row: offset_range_position(
            row["audit_trial_count"],
            trial_ranges.get(int(row["audit_participant_id"]))
            if pd.notna(row["audit_participant_id"])
            else None,
        ),
        axis=1,
    )
    raw_only["offset_trial_count_present"] = raw_only.apply(
        lambda row: (
            int(row["audit_trial_count"]) in trial_sets.get(int(row["audit_participant_id"]), set())
            if pd.notna(row["audit_participant_id"]) and pd.notna(row["audit_trial_count"])
            else False
        ),
        axis=1,
    )

    categories = raw_only.apply(raw_only_category, axis=1)
    raw_only["raw_only_category"] = [category for category, _ in categories]
    raw_only["raw_only_category_reason"] = [reason for _, reason in categories]

    preferred_columns = [
        "audit_participant_id",
        "audit_trial_count",
        "source_filename",
        "file_role",
        "split_part",
        "raw_trial_sequence",
        "join_trial_count",
        "raw_only_category",
        "raw_only_category_reason",
        "participant_has_old_compatible_offsets",
        "participant_old_compatible_offset_rows",
        "participant_old_compatible_offset_trial_range",
        "offset_range_position",
        "offset_trial_count_present",
        "best_candidate_offset_trial_count",
        "best_candidate_offset_delta",
        "best_candidate_timing_score_seconds",
        "best_candidate_trace_within_old_threshold",
        "best_candidate_stimulus_within_old_threshold",
        "best_candidate_clean_by_timing",
        "best_candidate_better_than_nominal",
        "raw_join_key_duplicate",
        "missing_expected_event_names",
        "extra_expected_event_names",
        "alignment_problem_codes",
        "trial_issue_codes",
    ]
    present_columns = [column for column in preferred_columns if column in raw_only.columns]
    return raw_only[present_columns].sort_values(
        ["audit_participant_id", "audit_trial_count", "source_filename"],
        kind="mergesort",
        na_position="last",
    )


def parse_summary_after_clean_timing_count(summary_text: str) -> int | None:
    """Parse the ``after`` clean-timing count from script 02 Markdown.

    Args:
        summary_text: Text from ``event_sequence_audit_summary.md``.

    Returns:
        Parsed integer count, or ``None`` if the table row is not found.

    Side effects:
        None.
    """

    for line in summary_text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("| after |"):
            continue
        cells = [cell.strip().replace(",", "") for cell in stripped.strip("|").split("|")]
        if len(cells) < 5:
            continue
        try:
            return int(cells[4])
        except ValueError:
            return None
    return None


def build_clean_timing_definition_check(join: pd.DataFrame, summary_text: str) -> pd.DataFrame:
    """Compare clean-timing counts from the Markdown, CSV flag, and strict rule.

    Args:
        join: Proposed offset-join audit from script 02.
        summary_text: Text from script 02 Markdown summary.

    Returns:
        DataFrame with named count definitions and a patch recommendation.

    Side effects:
        None.
    """

    summary_count = parse_summary_after_clean_timing_count(summary_text)
    csv_flag = true_mask(join, "clean_timing_row")
    strict_mask = strict_clean_timing_mask(join)
    no_duplicate_check_mask = (
        join["join_status"].eq(RAW_AND_OFFSET_STATUS)
        & ~true_mask(join, "old_trace_epoch_duration_mismatch")
        & ~true_mask(join, "old_stimulus_duration_mismatch")
        & join["missing_expected_event_names"].fillna("").eq("")
        & join["extra_expected_event_names"].fillna("").eq("")
    )
    flagged_but_not_strict = join[csv_flag & ~strict_mask].copy()
    flagged_but_not_strict_duplicate = flagged_but_not_strict[
        true_mask(flagged_but_not_strict, "raw_join_key_duplicate")
    ]

    strict_count = int(strict_mask.sum())
    csv_count = int(csv_flag.sum())
    no_duplicate_count = int(no_duplicate_check_mask.sum())
    patch_needed = csv_count != strict_count or summary_count != strict_count

    rows = [
        {
            "definition_name": "markdown_after_clean_timing_count",
            "count": summary_count if summary_count is not None else pd.NA,
            "matches_strict_definition": summary_count == strict_count if summary_count is not None else False,
            "note": "count parsed from event_sequence_audit_summary.md before-vs-after table",
        },
        {
            "definition_name": "csv_clean_timing_row_true_count",
            "count": csv_count,
            "matches_strict_definition": csv_count == strict_count,
            "note": "count of true clean_timing_row values in proposed_offset_join_audit.csv",
        },
        {
            "definition_name": "stricter_clean_timing_recomputed",
            "count": strict_count,
            "matches_strict_definition": True,
            "note": "matched rows with no timing mismatch, no missing or extra expected event, and no duplicate raw join key",
        },
        {
            "definition_name": "clean_count_without_duplicate_key_check",
            "count": no_duplicate_count,
            "matches_strict_definition": no_duplicate_count == strict_count,
            "note": "same timing/event checks but duplicate raw join keys are not subtracted",
        },
        {
            "definition_name": "csv_true_but_not_strict_rows",
            "count": int(len(flagged_but_not_strict)),
            "matches_strict_definition": int(len(flagged_but_not_strict)) == 0,
            "note": "rows marked clean_timing_row in the CSV but failing the stricter recomputation",
        },
        {
            "definition_name": "csv_true_but_not_strict_duplicate_key_rows",
            "count": int(len(flagged_but_not_strict_duplicate)),
            "matches_strict_definition": int(len(flagged_but_not_strict_duplicate)) == 0,
            "note": "subset explained by duplicate raw join keys",
        },
        {
            "definition_name": "small_patch_recommended",
            "count": int(patch_needed),
            "matches_strict_definition": not patch_needed,
            "note": (
                "yes: compute clean_timing_row with robust boolean coercion for raw_join_key_duplicate"
                if patch_needed
                else "no mismatch detected"
            ),
        },
    ]
    return pd.DataFrame(rows)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> list[str]:
    """Format selected DataFrame columns as a Markdown table.

    Args:
        frame: DataFrame to format.
        columns: Column names to render in order.

    Returns:
        Markdown table lines.

    Side effects:
        None.
    """

    if frame.empty:
        return ["_No rows._"]
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    body: list[str] = []
    for _, row in frame.iterrows():
        cells: list[str] = []
        for column in columns:
            value = row.get(column, "")
            if pd.isna(value):
                cells.append("")
            elif isinstance(value, float):
                cells.append(f"{value:.3f}")
            else:
                cells.append(str(value))
        body.append("| " + " | ".join(cells) + " |")
    return [header, divider, *body]


def category_counts_frame(classified_raw_only: pd.DataFrame) -> pd.DataFrame:
    """Count raw-only rows by classifier category.

    Args:
        classified_raw_only: Output of ``build_raw_only_row_classification``.

    Returns:
        DataFrame with ``raw_only_category`` and ``rows`` columns.

    Side effects:
        None.
    """

    if classified_raw_only.empty:
        return pd.DataFrame(columns=["raw_only_category", "rows"])
    return (
        classified_raw_only["raw_only_category"]
        .value_counts()
        .rename_axis("raw_only_category")
        .reset_index(name="rows")
        .sort_values(["rows", "raw_only_category"], ascending=[False, True], kind="mergesort")
    )


def build_summary(
    split_audit: pd.DataFrame,
    id5_audit: pd.DataFrame,
    zero_no_offset_audit: pd.DataFrame,
    raw_only_classification: pd.DataFrame,
    clean_timing_check: pd.DataFrame,
    first_mismatch: pd.DataFrame,
    started_at: datetime,
) -> str:
    """Build the Markdown summary for the focused special-case audit.

    Args:
        split_audit: Split-file continuity audit table.
        id5_audit: ID 5 concatenated-file audit table.
        zero_no_offset_audit: Zero/no-offset case table.
        raw_only_classification: Raw-only row classifier output.
        clean_timing_check: Clean-timing definition check table.
        first_mismatch: First-mismatch-by-file table from script 02.
        started_at: Timestamp captured at audit start.

    Returns:
        Markdown summary text.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    split_totals = split_audit[split_audit["row_kind"].eq("participant_total")].copy()
    id5_total = id5_audit[id5_audit["row_kind"].eq("participant_total")].copy()
    category_counts = category_counts_frame(raw_only_classification)

    strict_count_row = clean_timing_check[
        clean_timing_check["definition_name"].eq("stricter_clean_timing_recomputed")
    ]
    csv_count_row = clean_timing_check[
        clean_timing_check["definition_name"].eq("csv_clean_timing_row_true_count")
    ]
    patch_row = clean_timing_check[clean_timing_check["definition_name"].eq("small_patch_recommended")]
    strict_count = int(strict_count_row["count"].iloc[0]) if not strict_count_row.empty else 0
    csv_count = int(csv_count_row["count"].iloc[0]) if not csv_count_row.empty else 0
    patch_needed = bool(int(patch_row["count"].iloc[0])) if not patch_row.empty else False

    cascade_files = first_mismatch[
        first_mismatch["possible_cascade_after_first_problem"].map(is_true)
    ].copy()

    lines = [
        "# Event Alignment Special-Case Audit Summary",
        "",
        f"Generated: {finished_at.isoformat(timespec='seconds')}",
        f"Started: {started_at.isoformat(timespec='seconds')}",
        f"Duration seconds: {(finished_at - started_at).total_seconds():.1f}",
        "",
        "## Scope",
        "",
        "This local audit inspects remaining event-alignment blockers using the proposed read-only surface from script 02. It writes local CSV/Markdown outputs only and does not repair, concatenate, crop, preprocess, or epoch EEG data.",
        "",
        "## Outputs",
        "",
        f"- `{(EVENT_SPECIAL_CASE_AUDIT_DIR / SPLIT_FILE_CONTINUITY_AUDIT_FILENAME).as_posix()}`",
        f"- `{(EVENT_SPECIAL_CASE_AUDIT_DIR / ID5_CONCATENATED_FILE_AUDIT_FILENAME).as_posix()}`",
        f"- `{(EVENT_SPECIAL_CASE_AUDIT_DIR / ZERO_NO_OFFSET_CASE_AUDIT_FILENAME).as_posix()}`",
        f"- `{(EVENT_SPECIAL_CASE_AUDIT_DIR / RAW_ONLY_ROW_CLASSIFICATION_FILENAME).as_posix()}`",
        f"- `{(EVENT_SPECIAL_CASE_AUDIT_DIR / CLEAN_TIMING_DEFINITION_CHECK_FILENAME).as_posix()}`",
        "",
        "## Split Files",
        "",
        *markdown_table(
            split_totals,
            [
                "participant_id",
                "matched_rows",
                "raw_only_rows",
                "offset_only_rows",
                "duplicate_raw_join_key_rows",
                "strict_clean_timing_rows",
                "trace_mismatch_rows",
                "stimulus_mismatch_rows",
                "nearby_best_better_than_nominal_rows",
                "nearby_best_clean_by_timing_rows",
                "continuous_alignment_observable",
            ],
        ),
        "",
        "The split IDs still show duplicate nominal join keys, offset-only rows, and timing-candidate evidence that cross-part continuity is not represented by the current per-file trial counts. This audit does not concatenate the files.",
        "",
        "## ID 5",
        "",
        *markdown_table(
            id5_total,
            [
                "participant_id",
                "file_start_onset_seconds",
                "event_rows",
                "trial_group_rows",
                "matched_rows",
                "raw_only_rows",
                "strict_clean_timing_rows",
                "matched_rows_timing_clean_under_strict_definition",
            ],
        ),
        "",
        "ID 5's matched rows are timing-clean under the stricter definition, but the concatenated/duplicate-file structure and raw-only rows still need a reviewed event-handling policy. This audit does not crop or split ID 5.",
        "",
        "## Zero And No-Offset Cases",
        "",
        *markdown_table(
            zero_no_offset_audit,
            [
                "participant_id",
                "annotation_count_from_trial_summary",
                "trial_group_rows",
                "raw_only_rows",
                "offset_only_rows",
                "old_compatible_offset_rows",
                "old_compatible_offset_trial_range",
                "factual_note",
            ],
        ),
        "",
        "## Raw-Only Classifier",
        "",
        *markdown_table(category_counts, ["raw_only_category", "rows"]),
        "",
        "The classifier is row-level and conservative. It separates known split/concatenated rows, no-offset-target participants, structurally incomplete task-event rows, timing-candidate concerns, and likely old-compatible row-surface differences.",
        "",
        "## Clean Timing Definition",
        "",
        *markdown_table(
            clean_timing_check,
            ["definition_name", "count", "matches_strict_definition", "note"],
        ),
        "",
        "## Direct Answers",
        "",
        "- Are remaining mismatches mostly localized to known special cases? Severe timing and continuity problems are concentrated in the split-file IDs, and ID 5 remains a separate concatenated-file policy case. Raw-only rows are more widespread, but many classify as no-offset-target, known special-case, or likely old-compatible row-surface differences rather than direct duration mismatches.",
        "- Which cases need Tony decisions before epoch construction? IDs 54, 56, and 65 need a split-file continuity policy; ID 5 needs a concatenated/duplicate-file policy; IDs 11, 14, 89, 94, and 100 need explicit handling notes; raw-only categories need review before any event-driven epoching script depends on them.",
        f"- Is a small patch needed for clean-timing count reporting? {'Yes' if patch_needed else 'No'}. The CSV flag count is {csv_count:,}, while the stricter recomputation is {strict_count:,}.",
        "- Is high-level preprocessing planning still safe to begin? Yes, for planning that does not depend on final event repair or final epoch construction.",
        "- Is actual epoch construction still blocked? Yes. The split, concatenated, zero/no-offset, raw-only, and clean-timing-definition items above still need review.",
        "",
        "## Additional File-Level Context",
        "",
        f"- Files with possible cascade after first problem in script 02 first-mismatch table: {len(cascade_files):,}.",
        f"- Timing thresholds inherited from script 02: trace-end difference > {TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD * 1000:.0f} ms; stimulus-duration difference > {STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD * 1000:.0f} ms.",
        "",
        "## Safety Boundary",
        "",
        "- Proposed labels remain audit fields only.",
        "- No raw EDFs are read or modified.",
        "- No EEG signal data is loaded.",
        "- No preprocessing, filtering, or epoch construction is performed.",
        "- The generated files are local outputs under `_Data/`.",
        "",
    ]
    return "\n".join(lines)


def write_outputs(
    output_dir: Path,
    split_audit: pd.DataFrame,
    id5_audit: pd.DataFrame,
    zero_no_offset_audit: pd.DataFrame,
    raw_only_classification: pd.DataFrame,
    clean_timing_check: pd.DataFrame,
    summary_text: str,
) -> None:
    """Write focused special-case audit outputs.

    Args:
        output_dir: Local output directory below ``_Data/eeg``.
        split_audit: Split-file continuity audit table.
        id5_audit: ID 5 concatenated-file audit table.
        zero_no_offset_audit: Zero/no-offset case table.
        raw_only_classification: Raw-only row classifier table.
        clean_timing_check: Clean-timing definition check table.
        summary_text: Markdown summary text.

    Returns:
        None.

    Side effects:
        Creates ``output_dir`` when needed and writes CSV/Markdown audit files.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    split_audit.to_csv(output_dir / SPLIT_FILE_CONTINUITY_AUDIT_FILENAME, index=False)
    id5_audit.to_csv(output_dir / ID5_CONCATENATED_FILE_AUDIT_FILENAME, index=False)
    zero_no_offset_audit.to_csv(output_dir / ZERO_NO_OFFSET_CASE_AUDIT_FILENAME, index=False)
    raw_only_classification.to_csv(output_dir / RAW_ONLY_ROW_CLASSIFICATION_FILENAME, index=False)
    clean_timing_check.to_csv(output_dir / CLEAN_TIMING_DEFINITION_CHECK_FILENAME, index=False)
    (output_dir / EVENT_SPECIAL_CASE_AUDIT_SUMMARY_FILENAME).write_text(summary_text, encoding="utf-8")


def run_event_alignment_special_case_audit(repo_root: Path) -> dict[str, Any]:
    """Run the focused event-alignment special-case audit in memory.

    Args:
        repo_root: Absolute repository root path.

    Returns:
        Dictionary containing all output tables and summary text.

    Side effects:
        Reads local inputs below ``_Data/``. Does not write files.
    """

    started_at = datetime.now()
    inputs = read_audit_inputs(repo_root)
    events = inputs["events"]
    trials = inputs["trials"]
    join = inputs["join"]
    first_mismatch = inputs["first_mismatch"]
    nearby = inputs["nearby"]
    offsets = inputs["offsets"]
    summary_text = inputs["summary_text"]

    best_nearby = best_nearby_candidate_rows(nearby)
    split_audit = build_split_file_continuity_audit(events, trials, join, best_nearby)
    id5_audit = build_id5_concatenated_file_audit(events, trials, join)
    zero_no_offset_audit = build_zero_no_offset_case_audit(events, trials, join, offsets)
    raw_only_classification = build_raw_only_row_classification(join, offsets, best_nearby)
    clean_timing_check = build_clean_timing_definition_check(join, summary_text)
    special_case_summary = build_summary(
        split_audit=split_audit,
        id5_audit=id5_audit,
        zero_no_offset_audit=zero_no_offset_audit,
        raw_only_classification=raw_only_classification,
        clean_timing_check=clean_timing_check,
        first_mismatch=first_mismatch,
        started_at=started_at,
    )

    return {
        "split_audit": split_audit,
        "id5_audit": id5_audit,
        "zero_no_offset_audit": zero_no_offset_audit,
        "raw_only_classification": raw_only_classification,
        "clean_timing_check": clean_timing_check,
        "summary_text": special_case_summary,
    }


def main() -> None:
    """Run the focused event-alignment special-case audit and write outputs.

    Args:
        None. Paths are resolved from this script's repository location.

    Returns:
        None.

    Side effects:
        Reads script-02 CSV/Markdown outputs and old-compatible offsets,
        creates ``_Data/eeg/event_special_case_audit/`` if needed, writes
        local-only CSV/Markdown audit outputs, and prints a concise summary.
    """

    repo_root = repo_root_from_script()
    output_dir = repo_root / EVENT_SPECIAL_CASE_AUDIT_DIR

    try:
        outputs = run_event_alignment_special_case_audit(repo_root)
        write_outputs(
            output_dir=output_dir,
            split_audit=outputs["split_audit"],
            id5_audit=outputs["id5_audit"],
            zero_no_offset_audit=outputs["zero_no_offset_audit"],
            raw_only_classification=outputs["raw_only_classification"],
            clean_timing_check=outputs["clean_timing_check"],
            summary_text=outputs["summary_text"],
        )
    except Exception as error:
        raise SystemExit(f"FAIL event alignment special-case audit: {error}") from error

    raw_only_categories = category_counts_frame(outputs["raw_only_classification"])
    strict_row = outputs["clean_timing_check"][
        outputs["clean_timing_check"]["definition_name"].eq("stricter_clean_timing_recomputed")
    ]
    strict_count = int(strict_row["count"].iloc[0]) if not strict_row.empty else 0

    print(f"Wrote event-alignment special-case audit outputs to {output_dir}")
    print(f"Split-file audit rows: {len(outputs['split_audit']):,}")
    print(f"ID 5 audit rows: {len(outputs['id5_audit']):,}")
    print(f"Raw-only classified rows: {len(outputs['raw_only_classification']):,}")
    print(f"Raw-only categories: {len(raw_only_categories):,}")
    print(f"Stricter clean-timing rows: {strict_count:,}")


if __name__ == "__main__":
    main()
