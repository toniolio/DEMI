"""Reconstruct DEMI task-to-tracing event offsets before EEG alignment.

This script builds the behavioural/tracing timing table that the historical
DEMI EEG pipeline called ``epoch_offsets``. It reads local task text files and
the tracing-filter outputs produced by
``analysis/behavior/03_apply_tracing_filters.py``, joins task rows to tracing
summary rows by the old participant/session/block/trial logic, and constructs
the physical-tracing timing fields that later EEG trigger-alignment code needs.

Event offsets matter for the EEG reanalysis because raw EEG annotations should
not be compared to task rows or TraceLab samples until the behavioural/tracing
side has been reconstructed deterministically. The old EEG import code used
task timing, raw TraceLab trace onset/end, false-start shifts, and filtered
movement time to create real movement start and end times. This script ports
that pre-EEG timing step while keeping every row auditable.

The old authority source is ``_Scripts/04_import_eeg.R``. That script read the
post-cleanup behavioural table ``bdat2.rds``, raw ``tracings.rds``, and
``trace_filter_info.rds``; summarized raw trace onset/end; joined false-start
``start_shift`` values; computed ``trial_count = trial + (block - 1) * 20``;
and created ``real_start``, ``real_end``, and ``end_trigger``. The old
``_Scripts/_functions/eeg.R::update_events`` function later consumed those
offsets, but that EEG-trigger update is not ported here.

Old logic ported now:

- task TSV import with ``#`` comments ignored, matching the old readr import;
- old-R-style tracing join key from ``figure_file`` participant ID plus
  ``session_num``, ``block_num``, and ``trial_num``;
- preservation of task-file participant ID separately from the old tracing join
  ID, because old task filenames and task-row figure IDs can differ;
- trial-count reconstruction as ``trial + (block - 1) * 20``;
- raw trace onset/end fields from tracing-filter summaries;
- false-start ``start_shift`` recovery, defaulting missing shifts to zero as
  the old R code did;
- physical-trial ``real_start``, ``real_end``, and ``end_trigger`` formulas;
- imagery-trial offset formulas based on task movement time;
- row-level audit flags for no tracing payloads, tracing-filter removal, and
  task/tracing join gaps.

Work intentionally deferred:

- EEG annotation inspection and EEG preprocessing;
- the old ``update_events`` trigger repair/alignment logic;
- imagery mixture-model cleanup from ``_Scripts/02_imagery_clean.R``;
- final behavioural-analysis table reconstruction from
  ``_Scripts/03_behav_analysis.R``;
- accuracy/error metric reconstruction;
- participant-level decisions.

Expected local inputs:

- task text files below ``_Data/task/`` and ``_Data/task/incomplete/``;
- ``_Data/behavior/tracing_filters/tracing_trial_filter_summary.parquet``;
- ``_Data/behavior/tracing_filters/false_start_trials.parquet``;
- a Python environment with ``pandas`` and parquet support available.

Generated local outputs:

- ``_Data/behavior/event_offsets/event_offsets.csv``;
- ``_Data/behavior/event_offsets/event_offset_audit.csv``;
- ``_Data/behavior/event_offsets/task_file_read_audit.csv``;
- ``_Data/behavior/event_offsets/event_offset_summary.md``.

Safety boundaries:

- source task files and tracing-filter outputs are read-only inputs;
- outputs are written only below ``_Data/behavior/event_offsets/``;
- no raw EEG files or annotations are opened;
- no source data are modified;
- no participant-level decisions are made.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/04_build_event_offsets.py

Paths are resolved relative to this file, so the command also works from
another current working directory.
"""

from __future__ import annotations

import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


TASK_DIR = Path("_Data") / "task"
TRACING_FILTER_DIR = Path("_Data") / "behavior" / "tracing_filters"
EVENT_OFFSET_DIR = Path("_Data") / "behavior" / "event_offsets"

TRIALS_PER_BLOCK = 20

TASK_FILENAME_RE = re.compile(
    r"^p(?P<task_file_participant_id>\d{1,3})\.(?P<task_file_date>\d{4}-\d{2}-\d{2})(?P<label>.*)\.txt$",
    flags=re.IGNORECASE,
)
FIGURE_FILE_RE = re.compile(
    r"^p(?P<figure_file_participant_id>\d{1,3})_s(?P<figure_file_session>\d+)_"
    r"b(?P<figure_file_block>\d+)_t(?P<figure_file_trial>\d+)_"
    r"(?P<figure_file_date>\d{4}-\d{2}-\d{2})(?P<label>.*)\.tlf$",
    flags=re.IGNORECASE,
)

TASK_REQUIRED_COLUMNS = {
    "participant",
    "session_num",
    "block_num",
    "condition",
    "trial_num",
    "figure_file",
    "it",
    "mt",
    "stimulus_mt",
}

TRACE_REQUIRED_COLUMNS = {
    "participant_id",
    "session",
    "block",
    "trial",
    "tracings_parse_status",
    "tracings_row_count",
    "input_tracing_rows",
    "raw_trace_onset",
    "raw_trace_end",
    "raw_trace_duration",
    "final_tracing_rows",
    "mt_clip",
    "start_shift",
    "no_shape_trial",
    "incomplete_trial",
    "large_gap_trial",
    "hit_edge_trial",
    "no_tracing_payload",
    "has_reconstructed_tracing",
    "has_filtered_tracing",
}

SPECIAL_REVIEW_IDS = {
    "incomplete_task_marker": [2, 12, 20, 22, 27, 39, 88, 89, 99],
    "nonstandard_tracing_coverage": [2, 22, 27, 39, 88, 89, 99, 100, 999],
    "no_local_figure_zip_from_prior_review": [9, 10, 12, 13, 14, 20, 24, 26, 36, 78, 96],
}


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


def require_pyarrow() -> None:
    """Verify that parquet support is available in the active environment.

    Args:
        None.

    Returns:
        None.

    Side effects:
        Imports ``pyarrow``. Raises ``RuntimeError`` with setup context if the
        module is unavailable.
    """

    try:
        import pyarrow  # noqa: F401
    except ImportError as error:
        raise RuntimeError(
            "pyarrow is required to read tracing-filter parquet outputs. "
            "Run this script in the project environment."
        ) from error


def relative_to_repo(path: Path, repo_root: Path) -> str:
    """Return a POSIX path relative to the repository when possible.

    Args:
        path: Path to render.
        repo_root: Repository root.

    Returns:
        POSIX path text.

    Side effects:
        None.
    """

    try:
        return path.relative_to(repo_root).as_posix()
    except ValueError:
        return path.as_posix()


def discover_task_files(task_dir: Path) -> list[Path]:
    """Find task text files below the local task directory.

    Args:
        task_dir: Directory containing task text files and the incomplete
            subdirectory.

    Returns:
        Sorted list of ``.txt`` files.

    Side effects:
        Reads directory entries only.
    """

    if not task_dir.exists():
        return []
    return sorted(path for path in task_dir.rglob("*.txt") if path.is_file())


def parse_task_filename(path: Path, task_dir: Path) -> dict[str, Any]:
    """Parse task-file metadata from one local task path.

    Args:
        path: Task text path.
        task_dir: Root task directory.

    Returns:
        Dictionary with parsed file ID/date fields and incomplete-folder flags.

    Side effects:
        None.
    """

    relative_parts = path.relative_to(task_dir).parts
    in_incomplete_folder = "incomplete" in relative_parts[:-1]
    incomplete_style_name = "incomplete" in path.stem.lower()
    match = TASK_FILENAME_RE.match(path.name)

    if match is None:
        return {
            "task_file_participant_id": pd.NA,
            "task_file_date": "",
            "task_filename_parse_status": "task_filename_parse_error",
            "task_in_incomplete_folder": in_incomplete_folder,
            "task_incomplete_style_filename": incomplete_style_name,
            "task_incomplete_marker": in_incomplete_folder or incomplete_style_name,
        }

    return {
        "task_file_participant_id": int(match.group("task_file_participant_id")),
        "task_file_date": match.group("task_file_date"),
        "task_filename_parse_status": "ok",
        "task_in_incomplete_folder": in_incomplete_folder,
        "task_incomplete_style_filename": incomplete_style_name,
        "task_incomplete_marker": in_incomplete_folder or incomplete_style_name,
    }


def read_one_task_file(path: Path, task_dir: Path, repo_root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read one task TSV file and return row data plus file-level audit fields.

    Args:
        path: Task text file path.
        task_dir: Root task directory.
        repo_root: Repository root used for relative provenance paths.

    Returns:
        ``(rows, audit_row)``. ``rows`` may be empty when the file has only a
        header or cannot be read.

    Side effects:
        Opens the task file read-only. It does not modify the source file.
    """

    parsed = parse_task_filename(path, task_dir)
    audit_row: dict[str, Any] = {
        **parsed,
        "task_file_path": relative_to_repo(path, repo_root),
        "task_filename": path.name,
        "task_read_status": "not_read",
        "task_read_error_message": "",
        "task_row_count": 0,
        "task_column_count": 0,
        "task_columns": "",
    }

    try:
        frame = pd.read_csv(path, sep="\t", comment="#", dtype=str, engine="python")
    except Exception as error:  # noqa: BLE001 - keep one bad task file from stopping all file reads.
        audit_row["task_read_status"] = "task_read_error"
        audit_row["task_read_error_message"] = f"{type(error).__name__}: {error}"
        return pd.DataFrame(), audit_row

    audit_row["task_read_status"] = "ok"
    audit_row["task_row_count"] = int(len(frame))
    audit_row["task_column_count"] = int(len(frame.columns))
    audit_row["task_columns"] = ";".join(frame.columns)

    for column, value in audit_row.items():
        frame[column] = value
    frame["task_file_row_index"] = np.arange(1, len(frame) + 1, dtype="int64")
    return frame, audit_row


def read_task_rows(repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read all task rows and task-file audit metadata.

    Args:
        repo_root: Repository root.

    Returns:
        ``(task_rows, task_file_audit)``.

    Side effects:
        Opens local task text files read-only. Raises ``RuntimeError`` when no
        task files are found or required task columns are absent.
    """

    task_dir = repo_root / TASK_DIR
    task_files = discover_task_files(task_dir)
    if not task_files:
        raise RuntimeError(f"no task text files found under {relative_to_repo(task_dir, repo_root)}")

    row_frames: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    for path in task_files:
        rows, audit_row = read_one_task_file(path, task_dir, repo_root)
        audit_rows.append(audit_row)
        if not rows.empty:
            row_frames.append(rows)

    task_file_audit = pd.DataFrame(audit_rows)
    if not row_frames:
        raise RuntimeError("task files were found, but none contained task data rows")

    task_rows = pd.concat(row_frames, ignore_index=True, sort=False)
    missing_columns = sorted(TASK_REQUIRED_COLUMNS.difference(task_rows.columns))
    if missing_columns:
        raise RuntimeError(f"task rows are missing required column(s): {', '.join(missing_columns)}")

    return task_rows, task_file_audit


def parse_figure_file(value: Any) -> dict[str, Any]:
    """Parse old tracing join tokens from one task-row ``figure_file`` value.

    Args:
        value: Raw ``figure_file`` cell from a task TSV.

    Returns:
        Parsed participant/session/block/trial/date tokens and parse status.

    Side effects:
        None.
    """

    if pd.isna(value) or str(value).strip() == "":
        return {
            "figure_file_participant_id": pd.NA,
            "figure_file_session": pd.NA,
            "figure_file_block": pd.NA,
            "figure_file_trial": pd.NA,
            "figure_file_date": "",
            "figure_file_parse_status": "missing_figure_file",
        }

    text = str(value)
    match = FIGURE_FILE_RE.match(text)
    if match is None:
        return {
            "figure_file_participant_id": pd.NA,
            "figure_file_session": pd.NA,
            "figure_file_block": pd.NA,
            "figure_file_trial": pd.NA,
            "figure_file_date": "",
            "figure_file_parse_status": "figure_file_parse_error",
        }

    return {
        "figure_file_participant_id": int(match.group("figure_file_participant_id")),
        "figure_file_session": int(match.group("figure_file_session")),
        "figure_file_block": int(match.group("figure_file_block")),
        "figure_file_trial": int(match.group("figure_file_trial")),
        "figure_file_date": match.group("figure_file_date"),
        "figure_file_parse_status": "ok",
    }


def nullable_int(series: pd.Series) -> pd.Series:
    """Convert a series to pandas nullable integer dtype.

    Args:
        series: Input values.

    Returns:
        Numeric series with dtype ``Int64``. Values that cannot be parsed as
        whole numbers are converted to missing values so malformed local task
        rows remain visible in the audit output instead of stopping the run.

    Side effects:
        None.
    """

    numeric = pd.to_numeric(series, errors="coerce")
    whole_number = numeric.isna() | np.isclose(numeric % 1, 0)
    return numeric.where(whole_number).astype("Int64")


def numeric_float(series: pd.Series) -> pd.Series:
    """Convert a series to floating point values.

    Args:
        series: Input values.

    Returns:
        Floating point series with invalid values set to ``NaN``.

    Side effects:
        None.
    """

    return pd.to_numeric(series, errors="coerce").astype("float64")


def prepare_task_rows(task_rows: pd.DataFrame) -> pd.DataFrame:
    """Normalize task rows and create old-R-compatible join keys.

    Args:
        task_rows: Raw concatenated task rows from ``read_task_rows``.

    Returns:
        Task rows with parsed figure-file tokens, numeric timing fields, and
        explicit old-R join keys.

    Side effects:
        None.
    """

    out = task_rows.copy()
    parsed_figure = pd.DataFrame([parse_figure_file(value) for value in out["figure_file"]])
    out = pd.concat([out.reset_index(drop=True), parsed_figure], axis=1)

    # Old R inserted `db_id` from the task filename, but the tracing merge used
    # `fig_id` parsed from `figure_file`. Keep both so date-shifted file names
    # are visible while the old tracing join remains faithful.
    out["task_file_participant_id"] = nullable_int(out["task_file_participant_id"])
    out["old_epoch_id"] = nullable_int(out["participant"])
    out["trace_join_participant_id"] = nullable_int(out["figure_file_participant_id"])
    out["session_num"] = nullable_int(out["session_num"])
    out["block_num"] = nullable_int(out["block_num"])
    out["trial_num"] = nullable_int(out["trial_num"])
    out["figure_file_session"] = nullable_int(out["figure_file_session"])
    out["figure_file_block"] = nullable_int(out["figure_file_block"])
    out["figure_file_trial"] = nullable_int(out["figure_file_trial"])

    out["it"] = numeric_float(out["it"])
    out["task_mt"] = numeric_float(out["mt"])
    out["stimulus_mt"] = numeric_float(out["stimulus_mt"])

    out["physical"] = out["condition"].eq("physical")
    out["imagery"] = out["condition"].eq("imagery")
    out["trial_count"] = out["trial_num"] + ((out["block_num"] - 1) * TRIALS_PER_BLOCK)

    out["task_figure_key_matches_task_columns"] = (
        (out["figure_file_session"] == out["session_num"])
        & (out["figure_file_block"] == out["block_num"])
        & (out["figure_file_trial"] == out["trial_num"])
    )
    out.loc[out["figure_file_parse_status"] != "ok", "task_figure_key_matches_task_columns"] = False

    out["task_row_id"] = np.arange(1, len(out) + 1, dtype="int64")
    return out


def read_tracing_filter_summary(repo_root: Path) -> pd.DataFrame:
    """Read tracing-filter trial summaries for event-offset construction.

    Args:
        repo_root: Repository root.

    Returns:
        Trial-level tracing summary with old-R-compatible join-key columns.

    Side effects:
        Reads one parquet file below ``_Data/behavior/tracing_filters/``.
        Raises ``RuntimeError`` if the file or required columns are missing.
    """

    require_pyarrow()
    path = repo_root / TRACING_FILTER_DIR / "tracing_trial_filter_summary.parquet"
    if not path.exists():
        raise RuntimeError(
            f"required tracing filter summary is missing: {relative_to_repo(path, repo_root)}. "
            "Run analysis/behavior/03_apply_tracing_filters.py first."
        )

    try:
        summary = pd.read_parquet(path)
    except Exception as error:  # noqa: BLE001 - add path context.
        raise RuntimeError(f"failed to read {relative_to_repo(path, repo_root)}: {error}") from error

    missing_columns = sorted(TRACE_REQUIRED_COLUMNS.difference(summary.columns))
    if missing_columns:
        raise RuntimeError(f"tracing filter summary is missing column(s): {', '.join(missing_columns)}")

    out = summary.copy()
    out = out.rename(
        columns={
            "participant_id": "trace_join_participant_id",
            "session": "session_num",
            "block": "block_num",
            "trial": "trial_num",
            "raw_trace_onset": "trace_onset",
            "raw_trace_end": "trace_end",
        }
    )
    for column in ["trace_join_participant_id", "session_num", "block_num", "trial_num"]:
        out[column] = nullable_int(out[column])

    duplicate_count = len(out) - out[["trace_join_participant_id", "session_num", "block_num", "trial_num"]].drop_duplicates().shape[0]
    if duplicate_count:
        raise RuntimeError(f"tracing filter summary has {duplicate_count} duplicate trial join key(s)")

    return out


def join_task_to_tracing(task_rows: pd.DataFrame, tracing_summary: pd.DataFrame) -> pd.DataFrame:
    """Join task rows to tracing-filter trial summaries.

    Args:
        task_rows: Prepared task rows.
        tracing_summary: Trial-level tracing summaries.

    Returns:
        Joined DataFrame with a ``task_trace_join_status`` column.

    Side effects:
        None.
    """

    join_key = ["trace_join_participant_id", "session_num", "block_num", "trial_num"]
    joined = task_rows.merge(
        tracing_summary,
        on=join_key,
        how="left",
        validate="many_to_one",
        indicator="task_trace_join_status",
    )
    joined["task_trace_join_status"] = joined["task_trace_join_status"].map(
        {"both": "matched_tracing_summary", "left_only": "missing_tracing_summary", "right_only": "unexpected_right_only"}
    )
    return joined


def construct_old_offset_fields(joined: pd.DataFrame) -> pd.DataFrame:
    """Create old-compatible timing fields from task and tracing columns.

    Args:
        joined: Task rows joined to tracing-filter trial summaries.

    Returns:
        Copy with ``id``, ``session``, ``block``, ``trial``, ``mt``,
        ``real_start``, ``real_end``, and ``end_trigger`` fields.

    Side effects:
        None.
    """

    out = joined.copy()

    # `_Scripts/04_import_eeg.R` renamed task session/block/trial to
    # session/block/trial, then computed the sequential within-session trial
    # number. The old EEG updater later joined EEG trigger trial numbers to
    # `trial_count`.
    out["id"] = out["old_epoch_id"]
    out["session"] = out["session_num"]
    out["block"] = out["block_num"]
    out["trial"] = out["trial_num"]

    # Old R used filtered physical movement time (`mt_clip`) but task movement
    # time for imagery rows. This script does not rerun imagery cleanup.
    out["mt"] = np.where(out["physical"], out["mt_clip"], out["task_mt"])
    out["start_shift"] = out["start_shift"].fillna(0.0)

    physical_real_start = out["it"] + out["trace_onset"] + out["start_shift"]
    out["real_start"] = np.where(out["physical"], physical_real_start, 0.0)
    out["real_end"] = np.where(out["physical"], out["real_start"] + out["mt"], out["mt"])
    out["end_trigger"] = np.where(out["physical"], out["it"] + out["trace_end"], out["mt"])

    # Preserve the old final guard: when movement time is missing, the real
    # start should be missing even for imagery rows where it would otherwise be
    # zero.
    out.loc[out["mt"].isna(), "real_start"] = np.nan

    out["event_offset_constructed"] = out["real_start"].notna() & out["real_end"].notna()
    return out


def boolean_series(frame: pd.DataFrame, column: str) -> pd.Series:
    """Return a boolean series for an optional column.

    Args:
        frame: DataFrame that may or may not contain ``column``.
        column: Column name.

    Returns:
        Boolean Series aligned to ``frame``.

    Side effects:
        None.
    """

    if column not in frame.columns:
        return pd.Series(False, index=frame.index)
    return frame[column].fillna(False).astype(bool)


def issue_codes_for_row(row: pd.Series) -> str:
    """Build semicolon-delimited event-offset audit codes for one row.

    Args:
        row: Joined event-offset row.

    Returns:
        Semicolon-delimited issue/status code text. Empty text means no code.

    Side effects:
        None.
    """

    codes: list[str] = []

    if bool(row.get("task_incomplete_marker", False)):
        codes.append("task_incomplete_marker")
    if bool(row.get("task_in_incomplete_folder", False)):
        codes.append("task_file_outside_old_nonrecursive_task_import")
    if pd.isna(row.get("old_epoch_id")):
        codes.append("missing_task_participant")
    if row.get("figure_file_parse_status") != "ok":
        codes.append(str(row.get("figure_file_parse_status")))
    if row.get("figure_file_parse_status") == "ok" and not bool(row.get("task_figure_key_matches_task_columns", False)):
        codes.append("task_columns_differ_from_figure_file_tokens")
    if pd.notna(row.get("task_file_participant_id")) and pd.notna(row.get("trace_join_participant_id")):
        if int(row["task_file_participant_id"]) != int(row["trace_join_participant_id"]):
            codes.append("task_filename_id_differs_from_tracing_join_id")
    if row.get("task_trace_join_status") != "matched_tracing_summary":
        codes.append("missing_tracing_filter_summary")
    if bool(row.get("no_tracing_payload", False)):
        codes.append("no_tracing_payload")
    if bool(row.get("has_reconstructed_tracing", False)) and not bool(row.get("has_filtered_tracing", False)):
        codes.append("tracing_removed_by_filter_stage")
    for column, code in [
        ("no_shape_trial", "no_shape_trial"),
        ("incomplete_trial", "incomplete_tracing_trial"),
        ("large_gap_trial", "large_gap_trial"),
        ("hit_edge_trial", "edge_hit_trial"),
    ]:
        if bool(row.get(column, False)):
            codes.append(code)
    if row.get("condition") not in {"physical", "imagery"}:
        codes.append("unknown_condition")
    if bool(row.get("physical", False)):
        if pd.isna(row.get("it")):
            codes.append("physical_missing_it")
        if pd.isna(row.get("trace_onset")):
            codes.append("physical_missing_trace_onset")
        if pd.isna(row.get("trace_end")):
            codes.append("physical_missing_trace_end")
        if pd.isna(row.get("mt_clip")):
            codes.append("physical_missing_mt_clip")
    if bool(row.get("imagery", False)) and pd.isna(row.get("task_mt")):
        codes.append("imagery_missing_task_mt")
    if pd.notna(row.get("start_shift")) and float(row.get("start_shift")) != 0.0:
        codes.append("false_start_shift_present")
    if not bool(row.get("event_offset_constructed", False)):
        codes.append("event_offset_not_constructed")

    return ";".join(dict.fromkeys(code for code in codes if code))


def add_audit_fields(offsets: pd.DataFrame) -> pd.DataFrame:
    """Add row-level audit flags and issue-code text.

    Args:
        offsets: Event-offset rows with timing fields.

    Returns:
        Copy with audit columns.

    Side effects:
        None.
    """

    out = offsets.copy()
    out["has_tracing_summary_match"] = out["task_trace_join_status"].eq("matched_tracing_summary")
    out["tracing_trial_filtered_away"] = boolean_series(out, "has_reconstructed_tracing") & ~boolean_series(
        out, "has_filtered_tracing"
    )
    out["has_missing_tracing_data"] = (
        ~out["has_tracing_summary_match"]
        | boolean_series(out, "no_tracing_payload")
        | ~boolean_series(out, "has_reconstructed_tracing")
    )
    out["audit_issue_codes"] = out.apply(issue_codes_for_row, axis=1)
    return out


def event_offsets_columns(frame: pd.DataFrame) -> list[str]:
    """Return stable column order for the main event-offset output.

    Args:
        frame: Event-offset DataFrame.

    Returns:
        Existing columns in preferred order, followed by any remaining columns.

    Side effects:
        None.
    """

    preferred = [
        "id",
        "trial_count",
        "physical",
        "condition",
        "session",
        "block",
        "trial",
        "it",
        "mt",
        "stimulus_mt",
        "trace_onset",
        "trace_end",
        "start_shift",
        "real_start",
        "real_end",
        "end_trigger",
        "event_offset_constructed",
        "old_epoch_id",
        "trace_join_participant_id",
        "task_file_participant_id",
        "task_file_path",
        "task_file_row_index",
        "figure_file",
        "figure_file_parse_status",
        "task_trace_join_status",
        "has_missing_tracing_data",
        "tracing_trial_filtered_away",
        "no_tracing_payload",
        "has_reconstructed_tracing",
        "has_filtered_tracing",
        "mt_clip",
        "raw_trace_duration",
        "final_tracing_rows",
        "audit_issue_codes",
    ]
    return [column for column in preferred if column in frame.columns] + [
        column for column in frame.columns if column not in preferred
    ]


def audit_columns(frame: pd.DataFrame) -> list[str]:
    """Return stable column order for the row-level audit output.

    Args:
        frame: Event-offset DataFrame.

    Returns:
        Existing audit columns in preferred order.

    Side effects:
        None.
    """

    preferred = [
        "task_row_id",
        "task_file_path",
        "task_file_row_index",
        "task_incomplete_marker",
        "task_file_participant_id",
        "old_epoch_id",
        "trace_join_participant_id",
        "session_num",
        "block_num",
        "trial_num",
        "trial_count",
        "condition",
        "figure_file",
        "figure_file_parse_status",
        "task_figure_key_matches_task_columns",
        "task_trace_join_status",
        "has_missing_tracing_data",
        "no_tracing_payload",
        "has_reconstructed_tracing",
        "has_filtered_tracing",
        "tracing_trial_filtered_away",
        "no_shape_trial",
        "incomplete_trial",
        "large_gap_trial",
        "hit_edge_trial",
        "event_offset_constructed",
        "real_start",
        "real_end",
        "end_trigger",
        "audit_issue_codes",
    ]
    return [column for column in preferred if column in frame.columns]


def build_event_offsets(repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read inputs and build event-offset plus audit tables.

    Args:
        repo_root: Repository root.

    Returns:
        ``(event_offsets, event_offset_audit, task_file_audit)``.

    Side effects:
        Reads task files and tracing-filter parquet outputs. It does not write
        files.
    """

    task_rows_raw, task_file_audit = read_task_rows(repo_root)
    task_rows = prepare_task_rows(task_rows_raw)
    tracing_summary = read_tracing_filter_summary(repo_root)
    joined = join_task_to_tracing(task_rows, tracing_summary)
    offsets = construct_old_offset_fields(joined)
    audited = add_audit_fields(offsets)

    event_offsets = audited[event_offsets_columns(audited)].sort_values(
        ["id", "trial_count", "task_file_path", "task_file_row_index"],
        kind="mergesort",
        na_position="last",
    )
    event_offset_audit = audited[audit_columns(audited)].sort_values(
        ["old_epoch_id", "trial_count", "task_file_path", "task_file_row_index"],
        kind="mergesort",
        na_position="last",
    )
    return event_offsets.reset_index(drop=True), event_offset_audit.reset_index(drop=True), task_file_audit


def count_issue_codes(values: pd.Series) -> Counter[str]:
    """Count semicolon-delimited audit issue codes.

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


def format_id_list(values: list[int]) -> str:
    """Format participant IDs for Markdown.

    Args:
        values: Integer participant IDs.

    Returns:
        Comma-separated ID text or ``none``.

    Side effects:
        None.
    """

    return ", ".join(str(value) for value in sorted(set(values))) if values else "none"


def ids_present_in_offsets(offsets: pd.DataFrame, ids: list[int]) -> list[int]:
    """Find review-watchlist IDs present in event-offset rows.

    Args:
        offsets: Event-offset table.
        ids: Watchlist IDs.

    Returns:
        Sorted IDs present in any key-bearing participant column.

    Side effects:
        None.
    """

    present_values: set[int] = set()
    for column in ["old_epoch_id", "trace_join_participant_id", "task_file_participant_id"]:
        if column in offsets.columns:
            present_values.update(int(value) for value in offsets[column].dropna().astype(int).unique())
    return sorted(set(ids).intersection(present_values))


def build_summary(
    event_offsets: pd.DataFrame,
    event_offset_audit: pd.DataFrame,
    task_file_audit: pd.DataFrame,
    started_at: datetime,
) -> str:
    """Build the Markdown summary text for event-offset reconstruction.

    Args:
        event_offsets: Main event-offset output table.
        event_offset_audit: Row-level audit output table.
        task_file_audit: File-level task read audit table.
        started_at: Timestamp captured at run start.

    Returns:
        Markdown summary text.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    task_rows = len(event_offsets)
    top_level_rows = int((~event_offsets["task_in_incomplete_folder"].fillna(False).astype(bool)).sum())
    incomplete_folder_rows = int(event_offsets["task_in_incomplete_folder"].fillna(False).astype(bool).sum())
    joined_rows = int(event_offsets["has_tracing_summary_match"].sum())
    missing_tracing_rows = int(event_offsets["has_missing_tracing_data"].sum())
    filtered_away_rows = int(event_offsets["tracing_trial_filtered_away"].sum())
    constructed_rows = int(event_offsets["event_offset_constructed"].sum())
    physical_rows = int(event_offsets["physical"].fillna(False).astype(bool).sum())
    imagery_rows = int(event_offsets["imagery"].fillna(False).astype(bool).sum()) if "imagery" in event_offsets else int(
        event_offsets["condition"].eq("imagery").sum()
    )

    issue_counts = count_issue_codes(event_offset_audit["audit_issue_codes"])
    issue_lines = [f"- {code}: {count:,}" for code, count in sorted(issue_counts.items())] or ["- none"]

    condition_lines = []
    for condition, group in event_offsets.groupby("condition", dropna=False):
        condition_label = "missing" if pd.isna(condition) else str(condition)
        condition_lines.append(
            f"- {condition_label}: rows {len(group):,}; constructed offsets {int(group['event_offset_constructed'].sum()):,}"
        )

    watch_lines = []
    for label, ids in SPECIAL_REVIEW_IDS.items():
        watch_lines.append(f"- {label}: {format_id_list(ids_present_in_offsets(event_offsets, ids))}")

    task_file_status = task_file_audit["task_read_status"].value_counts(dropna=False).sort_index()
    task_file_lines = [f"- {status}: {count:,}" for status, count in task_file_status.items()]

    return f"""# Event Offset Reconstruction Summary

Generated: {finished_at.isoformat(timespec="seconds")}

Started: {started_at.isoformat(timespec="seconds")}

## Scope

This local summary reports reconstruction of task-to-tracing event offsets from
task TSV rows and tracing-filter outputs. It does not inspect raw EEG, update
EEG triggers, preprocess EEG, rerun imagery cleanup, or make participant-level
decisions.

## Task Files

- Task files read: {len(task_file_audit):,}
- Task rows read: {task_rows:,}
- Rows from top-level task files: {top_level_rows:,}
- Rows from `_Data/task/incomplete/`: {incomplete_folder_rows:,}

Task file read status:

{chr(10).join(task_file_lines)}

## Offset Rows

- Rows matched to tracing/filter summaries: {joined_rows:,}
- Rows with missing tracing data: {missing_tracing_rows:,}
- Rows whose tracing trial has reconstructed samples but no final filtered tracing rows: {filtered_away_rows:,}
- Rows with constructed `real_start` and `real_end`: {constructed_rows:,}
- Physical rows: {physical_rows:,}
- Imagery rows: {imagery_rows:,}

By task condition:

{chr(10).join(condition_lines)}

## Audit Issue Codes

{chr(10).join(issue_lines)}

## Known Special or Incomplete IDs

These IDs are surfaced from prior private reviews for linkage attention only.
They are not participant-level decisions.

{chr(10).join(watch_lines)}

## Old R Compatibility Notes

- The old tracing join used `fig_id` parsed from `figure_file`, not the task
  filename ID. This script preserves `task_file_participant_id`,
  `old_epoch_id`, and `trace_join_participant_id` separately.
- `trial_count` is reconstructed as `trial + (block - 1) * 20`, matching
  `_Scripts/04_import_eeg.R`.
- Physical rows use `it + trace_onset + start_shift` for `real_start`,
  filtered `mt_clip` for `mt`, and `it + trace_end` for `end_trigger`.
- Imagery rows use task `mt` for `real_end` and `end_trigger`, with
  `real_start = 0` when task `mt` is present.
- Missing false-start shifts are set to zero, matching the old R code.

## Deviations and TODOs

- The old `bdat2.rds` behavioural cleanup context is not rerun here.
- The imagery mixture-model cleanup from `_Scripts/02_imagery_clean.R` is not
  refit here.
- The old EEG trigger updater in `_Scripts/_functions/eeg.R::update_events` is
  not ported here.
- Event-offset rows with audit codes should be reviewed before raw EEG
  annotation comparison.
"""


def write_outputs(
    event_offsets: pd.DataFrame,
    event_offset_audit: pd.DataFrame,
    task_file_audit: pd.DataFrame,
    repo_root: Path,
    started_at: datetime,
) -> None:
    """Write event-offset outputs below the local output directory.

    Args:
        event_offsets: Main event-offset table.
        event_offset_audit: Row-level audit table.
        task_file_audit: File-level task read audit table.
        repo_root: Repository root.
        started_at: Timestamp captured at run start.

    Returns:
        None.

    Side effects:
        Creates ``_Data/behavior/event_offsets/`` if needed and writes CSV and
        Markdown files there.
    """

    output_dir = repo_root / EVENT_OFFSET_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    event_offsets_path = output_dir / "event_offsets.csv"
    event_offset_audit_path = output_dir / "event_offset_audit.csv"
    task_file_audit_path = output_dir / "task_file_read_audit.csv"
    summary_path = output_dir / "event_offset_summary.md"

    event_offsets.to_csv(event_offsets_path, index=False)
    event_offset_audit.to_csv(event_offset_audit_path, index=False)
    task_file_audit.to_csv(task_file_audit_path, index=False)
    summary_path.write_text(
        build_summary(event_offsets, event_offset_audit, task_file_audit, started_at),
        encoding="utf-8",
    )

    for path in [event_offsets_path, event_offset_audit_path, task_file_audit_path, summary_path]:
        print(f"Wrote {relative_to_repo(path, repo_root)}")


def main() -> int:
    """Run task-to-tracing event-offset reconstruction.

    Args:
        None.

    Returns:
        Process exit code: ``0`` on success, ``1`` on a handled failure.

    Side effects:
        Reads local task files and tracing-filter outputs, then writes local CSV
        and Markdown outputs under ``_Data/behavior/event_offsets/``.
    """

    started_at = datetime.now()
    repo_root = repo_root_from_script()

    try:
        event_offsets, event_offset_audit, task_file_audit = build_event_offsets(repo_root)
        write_outputs(event_offsets, event_offset_audit, task_file_audit, repo_root, started_at)
    except RuntimeError as error:
        print(f"FAIL event offsets: {error}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
