"""Inspect DEMI task files and figure zip filenames for linkage recovery.

This script performs a lightweight structural inspection of local DEMI task
text files and figure/tracing zip filenames. It is part of the inventory and
linkage-recovery phase of the EEG reanalysis. Its purpose is to make the
participant-level task/figure trial structure visible before any EEG event
reconstruction is attempted.

The previous behavioural file manifest,
``analysis/behavior/00_make_behavioral_file_manifest.py``, answers a first
question: which task files, figure directories, and zip counts are present
locally. This script takes the next narrow step. It checks whether each task
file has a recognizable tabular section, counts task data rows, records task
column names, and parses figure zip filename tokens that the old R import
script expected to represent participant ID, session, block, trial, and date.

The script records observable structure only:

- task ``.txt`` files below ``_Data/task/``, with files below
  ``_Data/task/incomplete/``;
- participant IDs and task dates parsed from task filenames;
- task line counts, comment-line counts, inferred header columns, and data-row
  counts after header inference;
- whether a task file can be read as a tab-delimited table at a structural
  level;
- recursive figure/tracing zip paths below ``_Data/figure/``;
- participant ID, session, block, trial, date token, learned-token status, and
  parent directory information inferred from zip filenames and paths;
- participant-level task/figure count summaries and factual linkage-risk
  flags.

This script does not parse behavioural outcome values, judge task performance,
open TraceLab zip contents, inspect EEG files, preprocess EEG, repair events,
or make participant-status decisions. Historical categories are carried as
compact labels only and are kept separate from directly observed file facts.

Expected inputs:

- task text files in ``_Data/task/`` and ``_Data/task/incomplete/``;
- figure/tracing zip files somewhere below ``_Data/figure/``;
- filenames following the observed DEMI pattern, for example
  ``p100_s1_b1_t10_2019-09-16.zip``;
- a Python environment with ``pandas`` available.

Generated local outputs:

- ``_Data/behavior/manifest/task_structure_manifest.csv``: one row per task
  file, with structural task metadata;
- ``_Data/behavior/manifest/figure_zip_filename_manifest.csv``: one row per
  zip file, with filename-token metadata;
- ``_Data/behavior/manifest/behavioral_linkage_summary.csv``: one row per
  expected or observed participant ID, with task/zip structural summaries and
  linkage-risk flags;
- ``_Data/behavior/manifest/behavioral_linkage_summary.md``: a concise local
  summary of structural counts and flagged linkage-risk groups.

Safety boundaries:

- task files are opened read-only;
- task values are loaded as strings solely so row and column structure can be
  counted;
- zip files are not opened;
- source files and directories are not modified;
- read or parse errors are recorded in output rows instead of stopping the
  whole run;
- generated outputs live under ``_Data/behavior/manifest/`` and should remain
  local-only.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/01_inspect_task_and_figure_linkage.py

The script also works when launched from another directory because paths are
resolved relative to this file's location in the repository.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


EXPECTED_PARTICIPANT_IDS = range(1, 101)
EXPECTED_TASK_ROWS_FROM_COMMON_PATTERN = 120
EXPECTED_NONLEARNED_ZIP_COUNT = 120
TRIALS_PER_BLOCK_FROM_TASK_COMMENTS = 20
DEV_TEST_PARTICIPANT_IDS = {999}

INCOMPLETE_TASK_MARKER_IDS = {2, 12, 20, 22, 27, 39, 88, 89, 99}
ZERO_ZIP_PARTICIPANT_IDS = {13}
FIGURE_WITHOUT_TASK_IDS = {100, 999}
TASK_WITHOUT_FIGURE_DIRECTORY_IDS = {9, 12, 14, 20, 26, 36, 78, 96}
DATE_MISMATCH_IDS_FROM_PRIVATE_REVIEW = {
    1,
    2,
    3,
    4,
    6,
    7,
    8,
    11,
    15,
    16,
    17,
    19,
    21,
    22,
    27,
    37,
    38,
    39,
    79,
    80,
    81,
    82,
    84,
    85,
    86,
    87,
    88,
    89,
    97,
    99,
}

# These compact labels mirror the first behavioural file manifest and keep
# historical notes separate from directly observed task/figure facts.
HISTORICAL_BEHAVIOURAL_SESSION_CATEGORIES = {
    9: "experimenter_error_or_wrong_protocol",
    10: "experimenter_error_or_wrong_protocol",
    12: "experimenter_error_or_wrong_protocol",
    14: "experimenter_error_or_wrong_protocol",
    20: "experimenter_error_or_wrong_protocol",
    24: "very_bad_eeg_session_note",
    89: "technical_noncompletion",
    96: "technical_noncompletion",
    100: "technical_noncompletion",
    13: "id_skipped",
    26: "id_skipped",
    36: "id_skipped",
    78: "id_skipped",
}

HISTORICAL_EEG_CATEGORIES = {
    6: "eeg_not_recorded_historical_note",
    7: "eeg_not_recorded_historical_note",
    11: "eeg_not_recorded_historical_note",
    24: "very_bad_eeg_historical_note",
    54: "triggers_missing_historical_note",
    56: "triggers_missing_historical_note",
    65: "software_crash_historical_note",
    86: "veol_or_veo_l_historical_note",
}

SPECIAL_CASE_CATEGORIES = {
    2: ["incomplete_task_marker"],
    3: ["eeg_quality_triage_note"],
    5: ["concatenated_or_duplicate_raw_eeg_review"],
    11: ["legacy_eeg_note_conflicts_with_local_raw_edf"],
    12: ["incomplete_task_marker"],
    20: ["incomplete_task_marker"],
    22: ["incomplete_task_marker"],
    27: ["incomplete_task_marker", "auxiliary_channel_or_quality_triage_note"],
    28: ["eeg_quality_triage_note"],
    30: ["eeg_quality_triage_note"],
    39: ["incomplete_task_marker"],
    54: ["split_raw_eeg_trigger_review"],
    56: ["split_raw_eeg_trigger_review"],
    61: ["auxiliary_channel_note"],
    64: ["eeg_quality_triage_note"],
    65: ["split_raw_eeg_software_crash_note"],
    69: ["auxiliary_channel_note"],
    84: ["auxiliary_channel_note"],
    86: ["auxiliary_channel_note"],
    88: ["incomplete_task_marker"],
    89: ["incomplete_task_marker", "partial_figure_zip_count"],
    94: ["raw_eeg_zero_annotation_review"],
    95: ["auxiliary_channel_note"],
    99: ["incomplete_task_marker"],
    100: ["partial_figure_zip_count"],
    999: ["dev_test_id"],
}

TASK_FILENAME_RE = re.compile(
    r"^p(?P<participant_id>\d{1,3})\.(?P<task_date>\d{4}-\d{2}-\d{2})(?P<label>.*)\.txt$",
    flags=re.IGNORECASE,
)
FIGURE_DIRECTORY_RE = re.compile(
    r"^p(?P<participant_id>\d{1,3})_(?P<figure_date>\d{4}-\d{2}-\d{2})(?P<label>.*)$",
    flags=re.IGNORECASE,
)
FIGURE_ZIP_FILENAME_RE = re.compile(
    r"^p(?P<participant_id>\d{1,3})_s(?P<session>\d+)_b(?P<block>\d+)_t(?P<trial>\d+)_"
    r"(?P<date_token>\d{4}-\d{2}-\d{2})(?P<extra>.*)\.zip$",
    flags=re.IGNORECASE,
)
SESSION_FOLDER_RE = re.compile(r"^session_(?P<session>\d+)$", flags=re.IGNORECASE)
DATE_TOKEN_RE = re.compile(r"\d{4}-\d{2}-\d{2}")
TIME_TOKEN_RE = re.compile(r"(?<!\d)\d{2}:\d{2}:\d{2}(?!\d)")


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Args:
        None.

    Returns:
        A ``Path`` pointing to the repository root. The script is expected to
        live at ``analysis/behavior/01_inspect_task_and_figure_linkage.py``.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def json_dump(value: Any) -> str:
    """Serialize a value into stable JSON text for CSV cells.

    Args:
        value: JSON-serializable Python value.

    Returns:
        Deterministic JSON text.

    Side effects:
        None.
    """

    return json.dumps(value, sort_keys=True)


def join_sorted_strings(values: list[Any]) -> str:
    """Join unique non-empty values in stable sorted order.

    Args:
        values: Values to clean and join.

    Returns:
        A semicolon-separated string. Missing values produce an empty string.

    Side effects:
        None.
    """

    cleaned = sorted({str(value) for value in values if value not in (None, "") and str(value) != "nan"})
    return ";".join(cleaned)


def join_sorted_ints(values: list[Any]) -> str:
    """Join unique integer values in numeric order.

    Args:
        values: Values that may contain integers, strings, or missing values.

    Returns:
        A semicolon-separated string of integers. Missing values produce an
        empty string.

    Side effects:
        None.
    """

    cleaned: set[int] = set()
    for value in values:
        if value in (None, ""):
            continue
        try:
            cleaned.add(int(value))
        except (TypeError, ValueError):
            continue
    return ";".join(str(value) for value in sorted(cleaned))


def discover_task_files(task_dir: Path) -> list[Path]:
    """Find task text files below the task-data directory.

    Args:
        task_dir: Directory expected to contain DEMI task text files.

    Returns:
        A sorted list of paths with a case-insensitive ``.txt`` suffix.

    Side effects:
        Reads directory entries only.
    """

    if not task_dir.exists():
        return []

    return sorted(path for path in task_dir.rglob("*") if path.is_file() and path.suffix.lower() == ".txt")


def parse_task_filename(task_path: Path, task_dir: Path) -> dict[str, Any]:
    """Parse task filename metadata and incomplete-marker fields.

    Args:
        task_path: Path to one task file.
        task_dir: Root task directory used to detect the incomplete folder.

    Returns:
        A dictionary with participant ID, task date, filename parse status, and
        incomplete-marker booleans.

    Side effects:
        None.
    """

    relative_parts = task_path.relative_to(task_dir).parts
    in_incomplete_folder = any(part.lower() == "incomplete" for part in relative_parts[:-1])
    incomplete_style_filename = "incomplete" in task_path.stem.lower()
    match = TASK_FILENAME_RE.match(task_path.name)

    if match is None:
        return {
            "participant_id": None,
            "participant_id_padded": "",
            "task_date": "",
            "task_filename_parse_status": "task_filename_parse_error",
            "task_in_incomplete_folder": in_incomplete_folder,
            "task_incomplete_style_filename": incomplete_style_filename,
            "incomplete_task_flag": in_incomplete_folder or incomplete_style_filename,
        }

    participant_id = int(match.group("participant_id"))
    return {
        "participant_id": participant_id,
        "participant_id_padded": f"{participant_id:03d}",
        "task_date": match.group("task_date"),
        "task_filename_parse_status": "ok",
        "task_in_incomplete_folder": in_incomplete_folder,
        "task_incomplete_style_filename": incomplete_style_filename,
        "incomplete_task_flag": in_incomplete_folder or incomplete_style_filename,
    }


def read_task_lines(task_path: Path) -> tuple[list[str], str, str]:
    """Read one task file as text with conservative error handling.

    Args:
        task_path: Path to one task file.

    Returns:
        ``(lines, read_status, read_error_message)``. On read error, ``lines``
        is empty and the error type/message are preserved.

    Side effects:
        Opens the task file read-only.
    """

    try:
        return task_path.read_text(encoding="utf-8", errors="replace").splitlines(), "ok", ""
    except Exception as error:  # noqa: BLE001 - one bad task file should not stop the run.
        return [], "task_read_error", f"{type(error).__name__}: {error}"


def infer_task_header(lines: list[str]) -> dict[str, Any]:
    """Infer the first tabular header line in a task file.

    Args:
        lines: Text lines from a task file.

    Returns:
        A dictionary with header line number, column names, comment-line count,
        blank-line count before header, and data-line count after the header.

    Side effects:
        None.
    """

    comment_line_count = 0
    blank_line_count_before_header = 0
    header_line_index: int | None = None

    # KLibs task files start with comment metadata. The first non-comment,
    # non-empty line is treated as the tabular header. This mirrors the old
    # readr call using comment="#", but records the inference explicitly.
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            if header_line_index is None:
                blank_line_count_before_header += 1
            continue
        if stripped.startswith("#"):
            if header_line_index is None:
                comment_line_count += 1
            continue
        header_line_index = index
        break

    if header_line_index is None:
        return {
            "comment_line_count_before_header": comment_line_count,
            "blank_line_count_before_header": blank_line_count_before_header,
            "header_line_number": None,
            "header_column_count": 0,
            "header_columns_json": "[]",
            "noncomment_data_line_count_after_header": 0,
            "header_inference_status": "no_header_line_found",
        }

    header_line = lines[header_line_index]
    columns = [column.strip() for column in header_line.split("\t")]
    data_lines_after_header = [
        line
        for line in lines[header_line_index + 1 :]
        if line.strip() and not line.lstrip().startswith("#")
    ]

    return {
        "comment_line_count_before_header": comment_line_count,
        "blank_line_count_before_header": blank_line_count_before_header,
        "header_line_number": header_line_index + 1,
        "header_column_count": len(columns),
        "header_columns_json": json_dump(columns),
        "noncomment_data_line_count_after_header": len(data_lines_after_header),
        "header_inference_status": "ok",
    }


def read_task_table_shape(task_path: Path) -> dict[str, Any]:
    """Read a task file as a tab-delimited table for structural shape only.

    Args:
        task_path: Path to one task file.

    Returns:
        A dictionary with pandas read status, row count, column count, and
        read-error message. Data values are read as strings and not interpreted.

    Side effects:
        Opens the task file read-only through pandas.
    """

    try:
        # The old R import used comment="#". Keep that behavior, but force
        # string dtype so this script does not classify or transform task
        # outcome values.
        task_df = pd.read_csv(task_path, sep="\t", comment="#", dtype=str, engine="python")
        return {
            "task_table_read_status": "ok",
            "task_data_row_count": int(len(task_df)),
            "task_table_column_count": int(len(task_df.columns)),
            "task_table_read_error_message": "",
        }
    except Exception as error:  # noqa: BLE001 - structural read errors are manifest facts.
        return {
            "task_table_read_status": "task_table_read_error",
            "task_data_row_count": None,
            "task_table_column_count": None,
            "task_table_read_error_message": f"{type(error).__name__}: {error}",
        }


def task_structure_issue_codes(row: dict[str, Any]) -> list[str]:
    """Create factual task-structure issue codes for one task row.

    Args:
        row: Task row after filename parsing and structural table inspection.

    Returns:
        A list of short issue codes.

    Side effects:
        None.
    """

    issue_codes: list[str] = []

    if row["task_filename_parse_status"] != "ok":
        issue_codes.append("task_filename_parse_error")
    if row["task_read_status"] != "ok":
        issue_codes.append("task_read_error")
    if row["header_inference_status"] != "ok":
        issue_codes.append("task_header_inference_issue")
    if row["task_table_read_status"] != "ok":
        issue_codes.append("task_table_read_issue")
    if row["incomplete_task_flag"]:
        issue_codes.append("incomplete_task_marker")

    task_rows = row.get("task_data_row_count")
    if task_rows == 0:
        issue_codes.append("zero_task_data_rows")
    elif task_rows is not None and task_rows != EXPECTED_TASK_ROWS_FROM_COMMON_PATTERN:
        issue_codes.append("nonstandard_task_data_row_count")

    return issue_codes


def task_row_for_file(task_path: Path, task_dir: Path, repo_root: Path) -> dict[str, Any]:
    """Build one task-structure manifest row.

    Args:
        task_path: Path to one task text file.
        task_dir: Root task directory used for incomplete-marker detection.
        repo_root: Repository root used for relative output paths.

    Returns:
        A dictionary containing filename metadata, header/row structure, and
        factual task-structure flags.

    Side effects:
        Opens the task file read-only.
    """

    parsed = parse_task_filename(task_path, task_dir)
    lines, read_status, read_error_message = read_task_lines(task_path)
    header_info = infer_task_header(lines)
    table_info = read_task_table_shape(task_path) if read_status == "ok" else {
        "task_table_read_status": "task_table_not_read",
        "task_data_row_count": None,
        "task_table_column_count": None,
        "task_table_read_error_message": "task text read failed",
    }

    row = {
        **parsed,
        "task_filename": task_path.name,
        "task_file_path": task_path.relative_to(repo_root).as_posix(),
        "task_line_count": len(lines) if read_status == "ok" else None,
        "task_read_status": read_status,
        "task_read_error_message": read_error_message,
        **header_info,
        **table_info,
    }
    row["task_structurally_parseable"] = (
        row["task_filename_parse_status"] == "ok"
        and row["task_read_status"] == "ok"
        and row["header_inference_status"] == "ok"
        and row["task_table_read_status"] == "ok"
        and row["header_column_count"] > 1
    )
    row["task_structure_issue_codes"] = ";".join(task_structure_issue_codes(row))
    return row


def build_task_structure_manifest(task_files: list[Path], task_dir: Path, repo_root: Path) -> pd.DataFrame:
    """Build the one-row-per-task-file structural manifest.

    Args:
        task_files: Sorted task text file paths.
        task_dir: Root task directory.
        repo_root: Repository root used for relative source paths.

    Returns:
        A pandas ``DataFrame`` with one row per task file.

    Side effects:
        Opens each task file read-only.
    """

    return pd.DataFrame([task_row_for_file(path, task_dir, repo_root) for path in task_files])


def discover_figure_zip_files(figure_dir: Path) -> list[Path]:
    """Find figure/tracing zip files below the figure directory.

    Args:
        figure_dir: Directory expected to contain DEMI figure/tracing zips.

    Returns:
        A sorted list of zip-file paths with case-insensitive suffix matching.

    Side effects:
        Reads directory entries only. Zip files are not opened.
    """

    if not figure_dir.exists():
        return []

    return sorted(path for path in figure_dir.rglob("*") if path.is_file() and path.suffix.lower() == ".zip")


def top_level_figure_directory(zip_path: Path, figure_dir: Path) -> str:
    """Return the first path component below ``_Data/figure/``.

    Args:
        zip_path: Path to one zip file.
        figure_dir: Root figure directory.

    Returns:
        The top-level directory name, or an empty string if the path cannot be
        represented relative to ``figure_dir``.

    Side effects:
        None.
    """

    try:
        return zip_path.relative_to(figure_dir).parts[0]
    except (ValueError, IndexError):
        return ""


def parse_figure_directory(directory_name: str) -> dict[str, Any]:
    """Parse participant ID and date from a top-level figure directory name.

    Args:
        directory_name: First path component below ``_Data/figure/``.

    Returns:
        A dictionary with directory participant ID, date, kind, and parse
        status.

    Side effects:
        None.
    """

    match = FIGURE_DIRECTORY_RE.match(directory_name)
    if match is None:
        return {
            "directory_participant_id": None,
            "directory_participant_id_padded": "",
            "figure_directory_date": "",
            "figure_directory_kind": "nonparticipant_directory",
            "figure_directory_parse_status": "figure_directory_parse_error",
        }

    participant_id = int(match.group("participant_id"))
    directory_kind = "dev_test_id" if participant_id in DEV_TEST_PARTICIPANT_IDS else "participant_directory"
    return {
        "directory_participant_id": participant_id,
        "directory_participant_id_padded": f"{participant_id:03d}",
        "figure_directory_date": match.group("figure_date"),
        "figure_directory_kind": directory_kind,
        "figure_directory_parse_status": "ok",
    }


def parse_session_folder(zip_path: Path) -> int | None:
    """Parse a ``session_<n>`` folder token from a zip path when present.

    Args:
        zip_path: Path to one figure/tracing zip file.

    Returns:
        The parsed session number from a parent folder, or ``None``.

    Side effects:
        None.
    """

    for part in zip_path.parts:
        match = SESSION_FOLDER_RE.match(part)
        if match is not None:
            return int(match.group("session"))
    return None


def old_import_numeric_stem(filename: str) -> str:
    """Recreate the old R import's numeric-stem simplification.

    Args:
        filename: Zip basename.

    Returns:
        Filename stem after removing lowercase letters and periods, mirroring
        ``gsub("[a-z\\.]", "", fname)`` from the old import script.

    Side effects:
        None.
    """

    stem = Path(filename).stem
    return re.sub(r"[a-z\.]", "", stem)


def parse_zip_filename(zip_path: Path, figure_dir: Path, repo_root: Path) -> dict[str, Any]:
    """Parse structural metadata from one figure/tracing zip filename.

    Args:
        zip_path: Path to one zip file.
        figure_dir: Root figure directory.
        repo_root: Repository root used for relative output paths.

    Returns:
        A dictionary containing parsed filename tokens, directory tokens, old
        import numeric-stem tokens, and parse issue codes.

    Side effects:
        None. Zip contents are not opened.
    """

    top_directory = top_level_figure_directory(zip_path, figure_dir)
    directory_info = parse_figure_directory(top_directory)
    filename = zip_path.name
    filename_lower = filename.lower()
    learned_flag = "learned" in filename_lower
    match = FIGURE_ZIP_FILENAME_RE.match(filename)

    date_tokens = DATE_TOKEN_RE.findall(filename)
    time_tokens = TIME_TOKEN_RE.findall(filename)
    numeric_stem = old_import_numeric_stem(filename)
    old_import_tokens = [token for token in numeric_stem.split("_") if token != ""]

    row: dict[str, Any] = {
        **directory_info,
        "zip_filename": filename,
        "zip_file_path": zip_path.relative_to(repo_root).as_posix(),
        "top_level_figure_directory": top_directory,
        "parent_session_folder": parse_session_folder(zip_path),
        "filename_has_learned": learned_flag,
        "filename_date_tokens": ";".join(date_tokens),
        "filename_time_tokens": ";".join(time_tokens),
        "old_import_numeric_stem": numeric_stem,
        "old_import_token_count": len(old_import_tokens),
        "old_import_tokens_json": json_dump(old_import_tokens),
    }

    if match is None:
        row.update(
            {
                "participant_id": None,
                "participant_id_padded": "",
                "session": None,
                "block": None,
                "trial": None,
                "zip_date_token": date_tokens[0] if date_tokens else "",
                "zip_filename_parse_status": "zip_filename_parse_error",
                "zip_filename_extra_text": "",
            }
        )
    else:
        participant_id = int(match.group("participant_id"))
        row.update(
            {
                "participant_id": participant_id,
                "participant_id_padded": f"{participant_id:03d}",
                "session": int(match.group("session")),
                "block": int(match.group("block")),
                "trial": int(match.group("trial")),
                "zip_date_token": match.group("date_token"),
                "zip_filename_parse_status": "ok",
                "zip_filename_extra_text": match.group("extra"),
            }
        )

    issue_codes: list[str] = []
    if row["zip_filename_parse_status"] != "ok":
        issue_codes.append("zip_filename_parse_error")
    if row["figure_directory_parse_status"] != "ok":
        issue_codes.append("figure_directory_parse_error")
    if row["participant_id"] is not None and row["directory_participant_id"] is not None:
        if row["participant_id"] != row["directory_participant_id"]:
            issue_codes.append("zip_directory_participant_id_mismatch")
    if row["parent_session_folder"] is not None and row["session"] is not None:
        if row["parent_session_folder"] != row["session"]:
            issue_codes.append("zip_parent_session_mismatch")
    if row["participant_id"] in DEV_TEST_PARTICIPANT_IDS or row["directory_participant_id"] in DEV_TEST_PARTICIPANT_IDS:
        issue_codes.append("dev_test_id")
    if learned_flag:
        issue_codes.append("learned_filename_token")

    row["zip_filename_issue_codes"] = ";".join(issue_codes)
    return row


def build_figure_zip_filename_manifest(zip_files: list[Path], figure_dir: Path, repo_root: Path) -> pd.DataFrame:
    """Build the one-row-per-zip filename manifest.

    Args:
        zip_files: Sorted figure/tracing zip paths.
        figure_dir: Root figure directory.
        repo_root: Repository root used for relative source paths.

    Returns:
        A pandas ``DataFrame`` with one row per zip file.

    Side effects:
        None. Zip contents are not opened.
    """

    return pd.DataFrame([parse_zip_filename(path, figure_dir, repo_root) for path in zip_files])


def historical_fields_for_participant(participant_id: int) -> dict[str, str]:
    """Return compact historical-note fields for one participant ID.

    Args:
        participant_id: Integer DEMI participant ID.

    Returns:
        A dictionary with historical behavioural/session category, historical
        EEG category, and special-case category. Missing categories are empty.

    Side effects:
        None.
    """

    return {
        "historical_behavioural_session_category": HISTORICAL_BEHAVIOURAL_SESSION_CATEGORIES.get(
            participant_id, ""
        ),
        "historical_eeg_category": HISTORICAL_EEG_CATEGORIES.get(participant_id, ""),
        "special_case_category": ";".join(SPECIAL_CASE_CATEGORIES.get(participant_id, [])),
    }


def aggregate_task_by_participant(task_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate task-structure rows to participant level.

    Args:
        task_df: One-row-per-task-file manifest.

    Returns:
        A pandas ``DataFrame`` with one row per parsed task participant ID.

    Side effects:
        None.
    """

    if task_df.empty:
        return pd.DataFrame()

    parsed_df = task_df.dropna(subset=["participant_id"]).copy()
    if parsed_df.empty:
        return pd.DataFrame()

    parsed_df["participant_id"] = parsed_df["participant_id"].astype(int)
    rows: list[dict[str, Any]] = []

    for participant_id, group in parsed_df.groupby("participant_id", sort=True):
        rows.append(
            {
                "participant_id": int(participant_id),
                "task_file_count": int(len(group)),
                "task_filename": join_sorted_strings(group["task_filename"].tolist()),
                "task_file_path": join_sorted_strings(group["task_file_path"].tolist()),
                "task_date": join_sorted_strings(group["task_date"].tolist()),
                "task_line_count": join_sorted_ints(group["task_line_count"].tolist()),
                "task_data_row_count": join_sorted_ints(group["task_data_row_count"].tolist()),
                "task_table_column_count": join_sorted_ints(group["task_table_column_count"].tolist()),
                "task_structurally_parseable": bool(group["task_structurally_parseable"].fillna(False).all()),
                "incomplete_task_flag": bool(group["incomplete_task_flag"].fillna(False).any()),
                "task_structure_issue_codes": join_issue_codes(group["task_structure_issue_codes"].tolist()),
            }
        )

    return pd.DataFrame(rows)


def aggregate_zips_by_participant(zip_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate zip filename rows to participant level.

    Args:
        zip_df: One-row-per-zip filename manifest.

    Returns:
        A pandas ``DataFrame`` with one row per parsed zip participant ID.

    Side effects:
        None.
    """

    if zip_df.empty:
        return pd.DataFrame()

    parsed_df = zip_df.dropna(subset=["participant_id"]).copy()
    if parsed_df.empty:
        return pd.DataFrame()

    parsed_df["participant_id"] = parsed_df["participant_id"].astype(int)
    rows: list[dict[str, Any]] = []

    for participant_id, group in parsed_df.groupby("participant_id", sort=True):
        nonlearned_group = group[~group["filename_has_learned"].fillna(False).astype(bool)]
        session_block_pairs = sorted(
            {
                f"{int(row.session)}:{int(row.block)}"
                for row in nonlearned_group.itertuples()
                if pd.notna(row.session) and pd.notna(row.block)
            }
        )
        session_block_trial_tuples = sorted(
            {
                f"{int(row.session)}:{int(row.block)}:{int(row.trial)}"
                for row in nonlearned_group.itertuples()
                if pd.notna(row.session) and pd.notna(row.block) and pd.notna(row.trial)
            }
        )

        rows.append(
            {
                "participant_id": int(participant_id),
                "figure_zip_count": int(len(group)),
                "nonlearned_figure_zip_count": int(len(nonlearned_group)),
                "learned_figure_zip_count": int(group["filename_has_learned"].fillna(False).astype(bool).sum()),
                "figure_directory_name": join_sorted_strings(group["top_level_figure_directory"].tolist()),
                "figure_directory_date": join_sorted_strings(group["figure_directory_date"].tolist()),
                "zip_date_tokens": join_sorted_strings(group["zip_date_token"].tolist()),
                "unique_sessions": join_sorted_ints(nonlearned_group["session"].tolist()),
                "unique_blocks": join_sorted_ints(nonlearned_group["block"].tolist()),
                "unique_trials": join_sorted_ints(nonlearned_group["trial"].tolist()),
                "unique_session_count": int(nonlearned_group["session"].dropna().nunique()),
                "unique_block_count": int(nonlearned_group["block"].dropna().nunique()),
                "unique_trial_count": int(nonlearned_group["trial"].dropna().nunique()),
                "unique_session_block_pair_count": len(session_block_pairs),
                "unique_session_block_trial_count": len(session_block_trial_tuples),
                "session_block_pairs": ";".join(session_block_pairs),
                "zip_filename_parse_error_count": int((group["zip_filename_parse_status"] != "ok").sum()),
                "zip_filename_issue_codes": join_issue_codes(group["zip_filename_issue_codes"].tolist()),
            }
        )

    return pd.DataFrame(rows)


def join_issue_codes(issue_values: list[Any]) -> str:
    """Join semicolon-separated issue-code cells into one stable string.

    Args:
        issue_values: Values that may contain semicolon-separated issue codes.

    Returns:
        A semicolon-separated string of unique issue codes.

    Side effects:
        None.
    """

    issue_codes: set[str] = set()
    for value in issue_values:
        if value in (None, "") or str(value) == "nan":
            continue
        issue_codes.update(part for part in str(value).split(";") if part)
    return ";".join(sorted(issue_codes))


def all_summary_participant_ids(task_df: pd.DataFrame, zip_df: pd.DataFrame) -> list[int]:
    """Return expected and observed participant IDs for the linkage summary.

    Args:
        task_df: One-row-per-task-file manifest.
        zip_df: One-row-per-zip filename manifest.

    Returns:
        Sorted IDs 1-100 plus any additional parsed IDs in task or zip rows.

    Side effects:
        None.
    """

    ids = set(EXPECTED_PARTICIPANT_IDS)
    for df in (task_df, zip_df):
        if df.empty or "participant_id" not in df.columns:
            continue
        for participant_id in df["participant_id"].dropna():
            ids.add(int(participant_id))
    return sorted(ids)


def first_int_from_summary_cell(value: Any) -> int | None:
    """Return the first integer from a semicolon-separated summary cell.

    Args:
        value: Cell value such as ``"120"`` or ``"8;10"``.

    Returns:
        The first parsed integer, or ``None``.

    Side effects:
        None.
    """

    if value in (None, "") or str(value) == "nan":
        return None
    first = str(value).split(";")[0]
    try:
        return int(first)
    except ValueError:
        return None


def linkage_flags_for_row(row: dict[str, Any]) -> list[str]:
    """Create participant-level structural linkage flags.

    Args:
        row: Participant-level summary row after task and zip aggregation.

    Returns:
        A list of factual linkage flags.

    Side effects:
        None.
    """

    flags: list[str] = []
    participant_id = int(row["participant_id"])
    task_count = int(row.get("task_file_count", 0) or 0)
    zip_count = int(row.get("figure_zip_count", 0) or 0)
    nonlearned_zip_count = int(row.get("nonlearned_figure_zip_count", 0) or 0)
    task_row_count = first_int_from_summary_cell(row.get("task_data_row_count", ""))

    if participant_id in INCOMPLETE_TASK_MARKER_IDS or row.get("incomplete_task_flag", False):
        flags.append("incomplete_task_marker")
    if task_count == 0:
        flags.append("missing_task_file")
    if zip_count == 0:
        flags.append("no_figure_zips")
    if task_count > 0 and zip_count == 0:
        flags.append("task_without_figure_zips")
    if task_count == 0 and zip_count > 0:
        flags.append("figure_zips_without_task")
    if participant_id in ZERO_ZIP_PARTICIPANT_IDS:
        flags.append("known_zero_zip_id")
    if participant_id in FIGURE_WITHOUT_TASK_IDS:
        flags.append("known_figure_without_task_id")
    if participant_id in TASK_WITHOUT_FIGURE_DIRECTORY_IDS:
        flags.append("known_task_without_figure_directory_id")
    if participant_id in DEV_TEST_PARTICIPANT_IDS:
        flags.append("dev_test_id")
    if participant_id in DATE_MISMATCH_IDS_FROM_PRIVATE_REVIEW:
        flags.append("date_mismatch_private_review")
    if not bool(row.get("task_structurally_parseable", False)) and task_count > 0:
        flags.append("task_structural_parse_issue")
    if nonlearned_zip_count not in (0, EXPECTED_NONLEARNED_ZIP_COUNT):
        flags.append("nonstandard_nonlearned_zip_count")
    if task_row_count is not None and task_row_count != EXPECTED_TASK_ROWS_FROM_COMMON_PATTERN:
        flags.append("nonstandard_task_data_row_count")
    if task_row_count is not None and nonlearned_zip_count > 0 and task_row_count != nonlearned_zip_count:
        flags.append("task_row_zip_count_mismatch")

    issue_cells = [
        row.get("task_structure_issue_codes", ""),
        row.get("zip_filename_issue_codes", ""),
    ]
    joined_existing = join_issue_codes(issue_cells)
    if joined_existing:
        flags.extend(joined_existing.split(";"))

    return sorted(set(flags))


def build_behavioral_linkage_summary(task_df: pd.DataFrame, zip_df: pd.DataFrame) -> pd.DataFrame:
    """Build the participant-level behavioural linkage summary.

    Args:
        task_df: One-row-per-task-file manifest.
        zip_df: One-row-per-zip filename manifest.

    Returns:
        A pandas ``DataFrame`` with one row per expected or observed
        participant ID.

    Side effects:
        None.
    """

    summary_df = pd.DataFrame({"participant_id": all_summary_participant_ids(task_df, zip_df)})
    task_by_id = aggregate_task_by_participant(task_df)
    zip_by_id = aggregate_zips_by_participant(zip_df)

    if not task_by_id.empty:
        summary_df = summary_df.merge(task_by_id, on="participant_id", how="left")
    if not zip_by_id.empty:
        summary_df = summary_df.merge(zip_by_id, on="participant_id", how="left")

    numeric_fill_zero_columns = [
        "task_file_count",
        "figure_zip_count",
        "nonlearned_figure_zip_count",
        "learned_figure_zip_count",
        "unique_session_count",
        "unique_block_count",
        "unique_trial_count",
        "unique_session_block_pair_count",
        "unique_session_block_trial_count",
        "zip_filename_parse_error_count",
    ]
    for column in numeric_fill_zero_columns:
        if column not in summary_df.columns:
            summary_df[column] = 0
        summary_df[column] = summary_df[column].fillna(0).astype(int)

    string_columns = [
        "task_filename",
        "task_file_path",
        "task_date",
        "task_line_count",
        "task_data_row_count",
        "task_table_column_count",
        "task_structure_issue_codes",
        "figure_directory_name",
        "figure_directory_date",
        "zip_date_tokens",
        "unique_sessions",
        "unique_blocks",
        "unique_trials",
        "session_block_pairs",
        "zip_filename_issue_codes",
    ]
    for column in string_columns:
        if column not in summary_df.columns:
            summary_df[column] = ""
        summary_df[column] = summary_df[column].fillna("")

    bool_columns = ["task_structurally_parseable", "incomplete_task_flag"]
    for column in bool_columns:
        if column not in summary_df.columns:
            summary_df[column] = False
        summary_df[column] = summary_df[column].fillna(False).astype(bool)

    summary_df["participant_id_padded"] = summary_df["participant_id"].map(lambda value: f"{int(value):03d}")
    historical_df = pd.DataFrame(
        [historical_fields_for_participant(int(pid)) for pid in summary_df["participant_id"]]
    )
    summary_df = pd.concat([summary_df.reset_index(drop=True), historical_df], axis=1)
    summary_df["linkage_risk_flags"] = summary_df.apply(
        lambda row: ";".join(linkage_flags_for_row(row.to_dict())),
        axis=1,
    )

    column_order = [
        "participant_id",
        "participant_id_padded",
        "task_file_count",
        "task_filename",
        "task_file_path",
        "task_date",
        "task_line_count",
        "task_data_row_count",
        "task_table_column_count",
        "task_structurally_parseable",
        "incomplete_task_flag",
        "figure_zip_count",
        "nonlearned_figure_zip_count",
        "learned_figure_zip_count",
        "figure_directory_name",
        "figure_directory_date",
        "zip_date_tokens",
        "unique_sessions",
        "unique_blocks",
        "unique_trials",
        "unique_session_count",
        "unique_block_count",
        "unique_trial_count",
        "unique_session_block_pair_count",
        "unique_session_block_trial_count",
        "session_block_pairs",
        "zip_filename_parse_error_count",
        "linkage_risk_flags",
        "task_structure_issue_codes",
        "zip_filename_issue_codes",
        "historical_behavioural_session_category",
        "historical_eeg_category",
        "special_case_category",
    ]
    return summary_df[column_order].sort_values("participant_id").reset_index(drop=True)


def count_semicolon_codes(values: pd.Series) -> Counter[str]:
    """Count semicolon-separated code strings from a pandas series.

    Args:
        values: Series containing semicolon-separated code cells.

    Returns:
        A ``Counter`` mapping each code to row count.

    Side effects:
        None.
    """

    counter: Counter[str] = Counter()
    for value in values.fillna(""):
        if not value:
            continue
        counter.update(part for part in str(value).split(";") if part)
    return counter


def participant_ids_with_flag(summary_df: pd.DataFrame, flag: str) -> list[int]:
    """Return participant IDs whose summary row contains a given risk flag.

    Args:
        summary_df: Participant-level linkage summary.
        flag: Risk flag to search for.

    Returns:
        A sorted list of participant IDs.

    Side effects:
        None.
    """

    ids: list[int] = []
    for _, row in summary_df.iterrows():
        flags = set(str(row["linkage_risk_flags"]).split(";")) if row["linkage_risk_flags"] else set()
        if flag in flags:
            ids.append(int(row["participant_id"]))
    return sorted(ids)


def format_id_list(ids: list[int]) -> str:
    """Format a list of participant IDs for Markdown output.

    Args:
        ids: Participant IDs.

    Returns:
        Comma-separated IDs, or ``none``.

    Side effects:
        None.
    """

    return ", ".join(str(value) for value in ids) if ids else "none"


def write_linkage_summary_markdown(
    summary_df: pd.DataFrame,
    task_df: pd.DataFrame,
    zip_df: pd.DataFrame,
    output_path: Path,
    task_dir: Path,
    figure_dir: Path,
) -> None:
    """Write a Markdown summary of behavioural linkage inspection outputs.

    Args:
        summary_df: Participant-level linkage summary.
        task_df: One-row-per-task-file manifest.
        zip_df: One-row-per-zip filename manifest.
        output_path: Destination Markdown path.
        task_dir: Task directory summarized by this run.
        figure_dir: Figure directory summarized by this run.

    Returns:
        ``None``.

    Side effects:
        Writes ``output_path``.
    """

    risk_counter = count_semicolon_codes(summary_df["linkage_risk_flags"])
    task_issue_counter = count_semicolon_codes(task_df["task_structure_issue_codes"]) if not task_df.empty else Counter()
    zip_issue_counter = count_semicolon_codes(zip_df["zip_filename_issue_codes"]) if not zip_df.empty else Counter()

    lines = [
        "# Behavioural Linkage Inspection Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Task directory: `{task_dir.as_posix()}`",
        f"Figure directory: `{figure_dir.as_posix()}`",
        f"Task files inspected: {len(task_df)}",
        f"Figure zip filenames inspected: {len(zip_df)}",
        f"Participant summary rows written: {len(summary_df)}",
        "",
        "## Participant-Level Risk Flag Counts",
        "",
    ]

    if risk_counter:
        for flag, count in sorted(risk_counter.items()):
            lines.append(f"- {flag}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Task Structure Issue Counts", ""])
    if task_issue_counter:
        for flag, count in sorted(task_issue_counter.items()):
            lines.append(f"- {flag}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Figure Filename Issue Counts", ""])
    if zip_issue_counter:
        for flag, count in sorted(zip_issue_counter.items()):
            lines.append(f"- {flag}: {count}")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Factual Linkage-Risk Groups",
            "",
            f"- Incomplete-task marker IDs: {format_id_list(participant_ids_with_flag(summary_df, 'incomplete_task_marker'))}",
            f"- Task file without figure zips: {format_id_list(participant_ids_with_flag(summary_df, 'task_without_figure_zips'))}",
            f"- Task file without figure directory from private review: {format_id_list(participant_ids_with_flag(summary_df, 'known_task_without_figure_directory_id'))}",
            f"- Figure zips without task file: {format_id_list(participant_ids_with_flag(summary_df, 'figure_zips_without_task'))}",
            f"- Zero figure zip IDs: {format_id_list(participant_ids_with_flag(summary_df, 'known_zero_zip_id'))}",
            f"- Nonstandard task row count IDs: {format_id_list(participant_ids_with_flag(summary_df, 'nonstandard_task_data_row_count'))}",
            f"- Nonstandard nonlearned zip count IDs: {format_id_list(participant_ids_with_flag(summary_df, 'nonstandard_nonlearned_zip_count'))}",
            f"- Task row / zip count mismatch IDs: {format_id_list(participant_ids_with_flag(summary_df, 'task_row_zip_count_mismatch'))}",
            f"- Date mismatch IDs from private review: {format_id_list(participant_ids_with_flag(summary_df, 'date_mismatch_private_review'))}",
            f"- Dev/test IDs: {format_id_list(participant_ids_with_flag(summary_df, 'dev_test_id'))}",
            "",
            "## Notes",
            "",
            "- This is a structural inspection, not a behavioural analysis.",
            "- Task table values are read as strings only to count rows and columns.",
            "- Zip files are not opened.",
            "- Risk flags are factual inspection flags, not participant-status decisions.",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(
    task_df: pd.DataFrame,
    zip_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_dir: Path,
    task_dir: Path,
    figure_dir: Path,
) -> None:
    """Write all linkage inspection outputs.

    Args:
        task_df: One-row-per-task-file manifest.
        zip_df: One-row-per-zip filename manifest.
        summary_df: Participant-level linkage summary.
        output_dir: Directory where outputs should be written.
        task_dir: Task directory summarized by this run.
        figure_dir: Figure directory summarized by this run.

    Returns:
        ``None``.

    Side effects:
        Creates ``output_dir`` when needed and writes four local-only output
        files.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    task_df.to_csv(output_dir / "task_structure_manifest.csv", index=False)
    zip_df.to_csv(output_dir / "figure_zip_filename_manifest.csv", index=False)
    summary_df.to_csv(output_dir / "behavioral_linkage_summary.csv", index=False)
    write_linkage_summary_markdown(
        summary_df=summary_df,
        task_df=task_df,
        zip_df=zip_df,
        output_path=output_dir / "behavioral_linkage_summary.md",
        task_dir=task_dir,
        figure_dir=figure_dir,
    )


def main() -> int:
    """Run the behavioural task/figure linkage inspection workflow.

    Args:
        None.

    Returns:
        ``0`` after local outputs are written.

    Side effects:
        Reads task files and figure zip paths, writes local manifest outputs
        under ``_Data/behavior/manifest/``, and prints a concise completion
        summary.
    """

    repo_root = repo_root_from_script()
    task_dir = repo_root / "_Data" / "task"
    figure_dir = repo_root / "_Data" / "figure"
    output_dir = repo_root / "_Data" / "behavior" / "manifest"

    task_files = discover_task_files(task_dir)
    zip_files = discover_figure_zip_files(figure_dir)

    task_df = build_task_structure_manifest(task_files, task_dir, repo_root)
    zip_df = build_figure_zip_filename_manifest(zip_files, figure_dir, repo_root)
    summary_df = build_behavioral_linkage_summary(task_df, zip_df)

    write_outputs(task_df, zip_df, summary_df, output_dir, task_dir, figure_dir)

    print(f"Wrote {output_dir / 'task_structure_manifest.csv'}")
    print(f"Wrote {output_dir / 'figure_zip_filename_manifest.csv'}")
    print(f"Wrote {output_dir / 'behavioral_linkage_summary.csv'}")
    print(f"Wrote {output_dir / 'behavioral_linkage_summary.md'}")
    print(f"Task files inspected: {len(task_df)}")
    print(f"Figure zip filenames inspected: {len(zip_df)}")
    print(f"Participant summary rows written: {len(summary_df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
