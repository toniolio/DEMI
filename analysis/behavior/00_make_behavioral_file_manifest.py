"""Create a local manifest of DEMI behavioural task and figure files.

This script inspects local behavioural task text files under ``_Data/task/``
and local figure/tracing zip files under ``_Data/figure/``. It writes local
manifest outputs under ``_Data/behavior/manifest/``. The manifest exists as an
early linkage-recovery step for the EEG reanalysis: before EEG annotations can
be linked back to task and tracing records, the project needs a plain inventory
of which behavioural source files are present locally, how participant IDs are
encoded in their filenames, and where obvious file-level irregularities appear.

The script records observable file facts only:

- task ``.txt`` files, with files currently stored in
  ``_Data/task/incomplete/``;
- participant IDs and task dates parsed from task filenames such as
  ``p1.2018-03-06.txt`` and ``p2.2018-03-16_incomplete.txt``;
- whether a task file is stored in the incomplete folder or has an
  incomplete-style filename;
- line counts for task files, without parsing behavioural outcomes;
- top-level figure directories under ``_Data/figure/``;
- participant IDs and dates parsed from figure directory names such as
  ``p100_2019-09-16 15:14:43``;
- recursive zip-file counts for each top-level figure directory;
- compact historical category labels recovered from private review documents,
  kept separate from directly observed file facts.

The behavioural/figure manifest exists because the current DEMI EEG reanalysis
is still in inventory and linkage-recovery mode. The old behavioural analysis
is treated as frozen. This script therefore does not re-run or revise the
published behavioural analysis. Its output is meant to help identify which
task and tracing files need later review before event reconstruction, not to
make participant-level analysis decisions.

The script does not parse full task outcomes, unzip or parse TraceLab contents,
read trial-level response traces, inspect EEG files, preprocess EEG, decide
participant status, or modify any source data. Missing files, incomplete-task
markers, nonstandard zip counts, dev/test IDs, and historical categories are
reported as inspection flags only.

Expected inputs:

- ``_Data/task/`` containing DEMI task ``.txt`` files;
- ``_Data/task/incomplete/`` optionally containing incomplete-task marker files;
- ``_Data/figure/`` containing participant-style figure directories and any
  non-participant local directories such as development or incomplete folders;
- a Python environment with ``pandas`` available.

Generated local outputs:

- ``_Data/behavior/manifest/behavioral_file_manifest.csv``: one row per
  expected or observed participant ID, with task and figure facts merged by ID;
- ``_Data/behavior/manifest/figure_zip_counts.csv``: one row per top-level
  figure directory, with non-participant/dev/test directories kept visible;
- ``_Data/behavior/manifest/behavioral_file_manifest_summary.md``: a short
  local summary of discovered files and factual inspection flags.

Safety boundaries:

- source task files and figure directories are read-only inputs;
- zip files are counted by path only and are not opened;
- task files are opened only to count text lines;
- read errors are recorded in manifest rows instead of stopping the run;
- outputs are written only under ``_Data/behavior/manifest/`` and should remain
  local-only;
- historical category labels are compact pointers for later human review, not
  current analysis decisions.

Run from the repository root with:

    python3 analysis/behavior/00_make_behavioral_file_manifest.py

The script also works when launched from another directory because paths are
resolved relative to this file's location in the repository.
"""

from __future__ import annotations

import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


# The behavioural study participant numbers are represented as IDs 1-100 in
# the local/private reconciliation reports. This range is used here only to
# make missing local task or figure files visible in the manifest. It is not a
# participant-retention rule.
EXPECTED_PARTICIPANT_IDS = range(1, 101)

# Most complete participant figure directories contain 120 zip files in the
# local inventory. Counts that differ from this value are factual review flags;
# the script does not infer whether a partial count can or cannot support later
# linkage work.
STANDARD_FIGURE_ZIP_COUNT = 120

# ID 999 appears in the local figure directory tree and is not treated as a
# study participant ID by the private inventory. Keep it explicit so it is easy
# to find in both the figure table and the merged participant table.
DEV_TEST_PARTICIPANT_IDS = {999}

# Historical behavioural/session categories recovered from old R-script
# comments and summarized in private inventory reports. These compact labels
# keep historical notes separate from observable file facts.
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

# Historical EEG categories recovered from the old EEG path and private
# decision-source review. These labels are historical context, not current MNE
# or behavioural decisions.
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

# Special-case categories summarize private inventory findings that should be
# visible during linkage recovery. Several IDs carry more than one compact
# label; values are joined with semicolons in the manifest.
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

# This set powers a simple issue code requested for known special-case IDs. ID
# 94 is listed even though the behavioural files themselves do not explain the
# raw EEG zero-annotation finding.
KNOWN_SPECIAL_CASE_IDS = set(SPECIAL_CASE_CATEGORIES)

TASK_FILENAME_RE = re.compile(
    r"^p(?P<participant_id>\d{1,3})\.(?P<task_date>\d{4}-\d{2}-\d{2})(?P<label>.*)\.txt$",
    flags=re.IGNORECASE,
)
FIGURE_DIRECTORY_RE = re.compile(
    r"^p(?P<participant_id>\d{1,3})_(?P<figure_date>\d{4}-\d{2}-\d{2})(?P<label>.*)$",
    flags=re.IGNORECASE,
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Args:
        None.

    Returns:
        A ``Path`` pointing to the repository root. The script is expected to
        live at ``analysis/behavior/00_make_behavioral_file_manifest.py``.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def discover_task_files(task_dir: Path) -> list[Path]:
    """Find task text files below the task-data directory.

    Args:
        task_dir: Directory expected to contain DEMI task ``.txt`` files and
            possibly an ``incomplete`` subdirectory.

    Returns:
        A sorted list of ``.txt`` files found recursively below ``task_dir``.
        The extension match is case-insensitive.

    Side effects:
        Reads the directory tree. It does not open or modify task files.
    """

    if not task_dir.exists():
        return []

    return sorted(path for path in task_dir.rglob("*") if path.is_file() and path.suffix.lower() == ".txt")


def parse_task_filename(task_path: Path, task_dir: Path) -> dict[str, Any]:
    """Parse participant ID and incomplete markers from a task filename.

    Args:
        task_path: Path to one task ``.txt`` file.
        task_dir: Root task directory. This is used to detect whether the file
            is under ``_Data/task/incomplete/``.

    Returns:
        A dictionary with parsed participant ID, task date, incomplete-folder
        flag, incomplete-style filename flag, and a parse-warning string.

    Side effects:
        None.
    """

    match = TASK_FILENAME_RE.match(task_path.name)
    relative_parts = task_path.relative_to(task_dir).parts

    # Incomplete-task status is represented in two observable ways locally:
    # storage under _Data/task/incomplete/ and an "_incomplete" filename token.
    # Keep both fields so later review can distinguish folder organization from
    # filename content.
    in_incomplete_folder = any(part.lower() == "incomplete" for part in relative_parts[:-1])
    incomplete_style_filename = "incomplete" in task_path.stem.lower()

    if match is None:
        return {
            "participant_id": None,
            "participant_id_padded": "",
            "task_date": "",
            "task_filename_parse_warning": "filename_does_not_match_expected_task_pattern",
            "task_in_incomplete_folder": in_incomplete_folder,
            "task_incomplete_style_filename": incomplete_style_filename,
            "incomplete_task_flag": in_incomplete_folder or incomplete_style_filename,
        }

    participant_id = int(match.group("participant_id"))
    return {
        "participant_id": participant_id,
        "participant_id_padded": f"{participant_id:03d}",
        "task_date": match.group("task_date"),
        "task_filename_parse_warning": "",
        "task_in_incomplete_folder": in_incomplete_folder,
        "task_incomplete_style_filename": incomplete_style_filename,
        "incomplete_task_flag": in_incomplete_folder or incomplete_style_filename,
    }


def count_text_lines(task_path: Path) -> tuple[int | None, str, str]:
    """Count lines in a task text file with conservative error handling.

    Args:
        task_path: Path to one task ``.txt`` file.

    Returns:
        A tuple of ``(line_count, read_status, read_error_message)``. A
        successfully read file returns ``read_status='ok'``. A read error
        returns ``line_count=None`` and preserves the exception type/message.

    Side effects:
        Opens the task file for reading. The file is not modified.
    """

    try:
        with task_path.open("r", encoding="utf-8", errors="replace") as file:
            return sum(1 for _ in file), "ok", ""
    except Exception as error:  # noqa: BLE001 - one bad text file should not stop the manifest.
        return None, "task_read_error", f"{type(error).__name__}: {error}"


def task_row_for_file(task_path: Path, task_dir: Path, repo_root: Path) -> dict[str, Any]:
    """Build one task-file manifest row.

    Args:
        task_path: Path to one task ``.txt`` file.
        task_dir: Root task directory used for incomplete-folder detection.
        repo_root: Repository root used to store readable relative paths.

    Returns:
        A dictionary containing filename-derived fields, line count, read
        status, and source path fields for one task file.

    Side effects:
        Opens the task file to count lines.
    """

    parsed = parse_task_filename(task_path, task_dir)
    line_count, read_status, read_error_message = count_text_lines(task_path)

    return {
        **parsed,
        "task_filename": task_path.name,
        "task_file_path": task_path.relative_to(repo_root).as_posix(),
        "task_line_count": line_count,
        "task_read_status": read_status,
        "task_read_error_message": read_error_message,
    }


def build_task_file_table(task_files: list[Path], task_dir: Path, repo_root: Path) -> pd.DataFrame:
    """Build the one-row-per-task-file table.

    Args:
        task_files: Sorted list of task text files.
        task_dir: Root task directory.
        repo_root: Repository root used for relative source paths.

    Returns:
        A pandas ``DataFrame`` with one row per discovered task file.

    Side effects:
        Opens each task file to count lines.
    """

    rows = [task_row_for_file(task_path, task_dir, repo_root) for task_path in task_files]
    return pd.DataFrame(rows)


def discover_figure_directories(figure_dir: Path) -> list[Path]:
    """Find top-level figure directories below ``_Data/figure/``.

    Args:
        figure_dir: Directory expected to contain participant-style figure
            directories and possible non-participant local folders.

    Returns:
        A sorted list of immediate child directories. The search is not
        recursive because each returned directory is counted as one top-level
        source container.

    Side effects:
        Reads the directory tree. It does not open zip files.
    """

    if not figure_dir.exists():
        return []

    return sorted(path for path in figure_dir.iterdir() if path.is_dir())


def parse_figure_directory_name(directory_path: Path) -> dict[str, Any]:
    """Parse participant ID and date from a top-level figure directory name.

    Args:
        directory_path: Top-level directory below ``_Data/figure/``.

    Returns:
        A dictionary with parsed participant ID, figure date, directory kind,
        and a parse-warning string. ID 999 is marked as ``dev_test_id`` even
        though its directory name follows the participant-style pattern.

    Side effects:
        None.
    """

    match = FIGURE_DIRECTORY_RE.match(directory_path.name)
    if match is None:
        return {
            "participant_id": None,
            "participant_id_padded": "",
            "figure_date": "",
            "figure_directory_kind": "nonparticipant_directory",
            "figure_directory_parse_warning": "directory_does_not_match_expected_participant_pattern",
        }

    participant_id = int(match.group("participant_id"))
    directory_kind = "dev_test_id" if participant_id in DEV_TEST_PARTICIPANT_IDS else "participant_directory"

    return {
        "participant_id": participant_id,
        "participant_id_padded": f"{participant_id:03d}",
        "figure_date": match.group("figure_date"),
        "figure_directory_kind": directory_kind,
        "figure_directory_parse_warning": "",
    }


def count_recursive_zip_files(directory_path: Path) -> int:
    """Count zip files recursively below one figure directory.

    Args:
        directory_path: Top-level figure directory to inspect.

    Returns:
        The number of files below ``directory_path`` with a case-insensitive
        ``.zip`` extension.

    Side effects:
        Reads the directory tree. Zip files are counted by path only and are
        not opened.
    """

    return sum(1 for path in directory_path.rglob("*") if path.is_file() and path.suffix.lower() == ".zip")


def figure_issue_codes(row: dict[str, Any]) -> list[str]:
    """Create factual inspection flags for one figure-directory row.

    Args:
        row: A figure-directory row after name parsing and zip counting.

    Returns:
        A list of short issue codes describing observable directory facts.

    Side effects:
        None.
    """

    issue_codes: list[str] = []
    zip_count = row["figure_zip_count"]

    if row["figure_directory_parse_warning"]:
        issue_codes.append("figure_directory_parse_warning")

    if row["figure_directory_kind"] == "dev_test_id":
        issue_codes.append("dev_test_id")
    elif row["figure_directory_kind"] == "nonparticipant_directory":
        issue_codes.append("nonparticipant_directory")

    # Keep zero-zip and nonstandard nonzero counts as separate facts. A zero
    # count usually means a directory placeholder; a nonzero partial count
    # points to a different review question.
    if zip_count == 0:
        issue_codes.append("zero_zip_files")
    elif zip_count != STANDARD_FIGURE_ZIP_COUNT:
        issue_codes.append("nonstandard_zip_count")

    participant_id = row["participant_id"]
    if participant_id in KNOWN_SPECIAL_CASE_IDS:
        issue_codes.append("known_special_case_id")

    return issue_codes


def figure_row_for_directory(directory_path: Path, repo_root: Path) -> dict[str, Any]:
    """Build one figure-directory count row.

    Args:
        directory_path: Top-level figure directory below ``_Data/figure/``.
        repo_root: Repository root used to store a readable relative path.

    Returns:
        A dictionary containing parsed directory fields, recursive zip count,
        and semicolon-separated issue codes.

    Side effects:
        Recursively reads paths below ``directory_path`` to count zip files.
        Zip files are not opened or modified.
    """

    parsed = parse_figure_directory_name(directory_path)
    row = {
        **parsed,
        "figure_directory_name": directory_path.name,
        "figure_directory_path": directory_path.relative_to(repo_root).as_posix(),
        "figure_zip_count": count_recursive_zip_files(directory_path),
    }
    row["figure_issue_codes"] = ";".join(figure_issue_codes(row))
    return row


def build_figure_directory_table(figure_directories: list[Path], repo_root: Path) -> pd.DataFrame:
    """Build the one-row-per-top-level-figure-directory table.

    Args:
        figure_directories: Sorted list of top-level figure directories.
        repo_root: Repository root used for relative source paths.

    Returns:
        A pandas ``DataFrame`` with one row per top-level figure directory.

    Side effects:
        Recursively reads each figure directory to count zip files. Zip files
        are not opened or modified.
    """

    rows = [figure_row_for_directory(directory_path, repo_root) for directory_path in figure_directories]
    return pd.DataFrame(rows)


def join_values(values: pd.Series, missing_value: str = "") -> str:
    """Join unique non-empty values from a pandas series.

    Args:
        values: Series containing values from one participant group.
        missing_value: Value returned when no non-empty values are present.

    Returns:
        A semicolon-separated string of unique values, preserving sorted order
        for stable CSV output.

    Side effects:
        None.
    """

    cleaned = sorted({str(value) for value in values.dropna() if str(value) not in {"", "nan", "None"}})
    return ";".join(cleaned) if cleaned else missing_value


def join_bool(values: pd.Series) -> bool:
    """Collapse a boolean-like pandas series with logical OR.

    Args:
        values: Series containing boolean-like values.

    Returns:
        ``True`` if any value is truthy, otherwise ``False``.

    Side effects:
        None.
    """

    return bool(values.fillna(False).astype(bool).any())


def join_int_values(values: pd.Series) -> str:
    """Join integer values from a pandas series for participant-level output.

    Args:
        values: Series containing line counts or zip counts.

    Returns:
        A semicolon-separated string of integer values. Missing values are
        omitted. A string is used because duplicate files for one participant
        would otherwise make a single numeric value misleading.

    Side effects:
        None.
    """

    cleaned: list[str] = []
    for value in values.dropna():
        cleaned.append(str(int(value)))
    return ";".join(cleaned)


def aggregate_task_files(task_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse task-file rows to one row per parsed participant ID.

    Args:
        task_df: One-row-per-task-file table.

    Returns:
        A participant-level task table. Multiple files for one ID are retained
        as semicolon-separated values and flagged later.

    Side effects:
        None.
    """

    if task_df.empty:
        return pd.DataFrame(
            columns=[
                "participant_id",
                "task_file_count",
                "task_file_status",
                "task_filename",
                "task_file_path",
                "task_date",
                "task_line_count",
                "task_in_incomplete_folder",
                "task_incomplete_style_filename",
                "incomplete_task_flag",
                "task_filename_parse_warning",
                "task_read_status",
                "task_read_error_message",
            ]
        )

    parsed_df = task_df.dropna(subset=["participant_id"]).copy()
    if parsed_df.empty:
        return pd.DataFrame()

    parsed_df["participant_id"] = parsed_df["participant_id"].astype(int)
    grouped = parsed_df.groupby("participant_id", as_index=False)

    aggregated = grouped.agg(
        task_file_count=("task_filename", "size"),
        task_filename=("task_filename", join_values),
        task_file_path=("task_file_path", join_values),
        task_date=("task_date", join_values),
        task_line_count=("task_line_count", join_int_values),
        task_in_incomplete_folder=("task_in_incomplete_folder", join_bool),
        task_incomplete_style_filename=("task_incomplete_style_filename", join_bool),
        incomplete_task_flag=("incomplete_task_flag", join_bool),
        task_filename_parse_warning=("task_filename_parse_warning", join_values),
        task_read_status=("task_read_status", join_values),
        task_read_error_message=("task_read_error_message", join_values),
    )
    aggregated["task_file_status"] = aggregated["task_file_count"].map(
        lambda count: "present" if count == 1 else "multiple_files"
    )
    return aggregated


def aggregate_figure_directories(figure_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse participant-style figure directory rows to one row per ID.

    Args:
        figure_df: One-row-per-top-level-figure-directory table.

    Returns:
        A participant-level figure table. Non-participant directories without
        parsed IDs are left out of this merged table and remain visible in
        ``figure_zip_counts.csv``.

    Side effects:
        None.
    """

    if figure_df.empty:
        return pd.DataFrame(
            columns=[
                "participant_id",
                "figure_directory_count",
                "figure_directory_status",
                "figure_directory_name",
                "figure_directory_path",
                "figure_date",
                "figure_zip_count",
                "figure_directory_kind",
                "figure_directory_parse_warning",
            ]
        )

    parsed_df = figure_df.dropna(subset=["participant_id"]).copy()
    if parsed_df.empty:
        return pd.DataFrame()

    parsed_df["participant_id"] = parsed_df["participant_id"].astype(int)
    grouped = parsed_df.groupby("participant_id", as_index=False)

    aggregated = grouped.agg(
        figure_directory_count=("figure_directory_name", "size"),
        figure_directory_name=("figure_directory_name", join_values),
        figure_directory_path=("figure_directory_path", join_values),
        figure_date=("figure_date", join_values),
        figure_zip_count=("figure_zip_count", join_int_values),
        figure_directory_kind=("figure_directory_kind", join_values),
        figure_directory_parse_warning=("figure_directory_parse_warning", join_values),
    )
    aggregated["figure_directory_status"] = aggregated["figure_directory_count"].map(
        lambda count: "present" if count == 1 else "multiple_directories"
    )
    return aggregated


def historical_fields_for_participant(participant_id: int) -> dict[str, str]:
    """Return compact historical-note fields for one participant ID.

    Args:
        participant_id: Integer DEMI participant ID.

    Returns:
        A dictionary containing historical behavioural/session category,
        historical EEG category, and special-case category. Empty strings mean
        no compact category was recovered for that ID.

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


def parse_first_date(date_text: str) -> str:
    """Return the first date from a semicolon-separated date field.

    Args:
        date_text: A date string or semicolon-separated date strings.

    Returns:
        The first non-empty date token, or an empty string if no date is
        available.

    Side effects:
        None.
    """

    for part in str(date_text).split(";"):
        if part:
            return part
    return ""


def observable_mismatch_flags(row: dict[str, Any]) -> list[str]:
    """Create observable task/figure mismatch flags for one merged row.

    Args:
        row: Participant-level manifest row containing task and figure status
            fields.

    Returns:
        A list of short mismatch codes. These codes describe file-presence,
        duplicate-file, and parsed-date differences only.

    Side effects:
        None.
    """

    mismatch_flags: list[str] = []

    if row["task_file_status"] == "missing" and row["figure_directory_status"] != "missing":
        mismatch_flags.append("figure_without_task_file")
    if row["task_file_status"] != "missing" and row["figure_directory_status"] == "missing":
        mismatch_flags.append("task_file_without_figure_directory")
    if row["task_file_status"] == "multiple_files":
        mismatch_flags.append("multiple_task_files")
    if row["figure_directory_status"] == "multiple_directories":
        mismatch_flags.append("multiple_figure_directories")

    # Dates are not used as join keys, but a task-date/figure-date difference
    # is an observable clue that a session was restarted, exported later, or
    # otherwise named differently. Keep it as a review flag only.
    task_date = parse_first_date(row.get("task_date", ""))
    figure_date = parse_first_date(row.get("figure_date", ""))
    if task_date and figure_date and task_date != figure_date:
        mismatch_flags.append("task_figure_date_mismatch")

    return mismatch_flags


def participant_issue_codes(row: dict[str, Any]) -> list[str]:
    """Create factual inspection flags for one participant-level row.

    Args:
        row: Participant-level manifest row after task/figure merging and
            historical category lookup.

    Returns:
        A list of short issue codes. These are file-inspection flags, not
        participant-level analysis decisions.

    Side effects:
        None.
    """

    issue_codes: list[str] = []
    participant_id = row["participant_id"]

    if row["task_file_status"] == "missing":
        issue_codes.append("missing_task_file")
    elif row["task_file_status"] == "multiple_files":
        issue_codes.append("multiple_task_files")

    if row["figure_directory_status"] == "missing":
        issue_codes.append("missing_figure_directory")
    elif row["figure_directory_status"] == "multiple_directories":
        issue_codes.append("multiple_figure_directories")

    if row["incomplete_task_flag"]:
        issue_codes.append("incomplete_task_marker")

    if participant_id in DEV_TEST_PARTICIPANT_IDS:
        issue_codes.append("dev_test_id")

    if participant_id in KNOWN_SPECIAL_CASE_IDS:
        issue_codes.append("known_special_case_id")

    if row.get("task_read_status", "") and row["task_read_status"] != "ok":
        issue_codes.append("task_read_issue")

    zip_count_text = str(row.get("figure_zip_count", ""))
    if zip_count_text:
        zip_counts = [int(part) for part in zip_count_text.split(";") if part]
        if any(count == 0 for count in zip_counts):
            issue_codes.append("zero_zip_files")
        if any(count != 0 and count != STANDARD_FIGURE_ZIP_COUNT for count in zip_counts):
            issue_codes.append("nonstandard_zip_count")

    issue_codes.extend(observable_mismatch_flags(row))
    return issue_codes


def all_manifest_participant_ids(task_df: pd.DataFrame, figure_df: pd.DataFrame) -> list[int]:
    """Return expected and observed participant IDs for the merged manifest.

    Args:
        task_df: One-row-per-task-file table.
        figure_df: One-row-per-figure-directory table.

    Returns:
        A sorted list containing IDs 1-100 plus any additional parsed IDs found
        in task or figure files, such as local dev/test ID 999.

    Side effects:
        None.
    """

    ids = set(EXPECTED_PARTICIPANT_IDS)

    for df in (task_df, figure_df):
        if df.empty or "participant_id" not in df.columns:
            continue
        for participant_id in df["participant_id"].dropna():
            ids.add(int(participant_id))

    return sorted(ids)


def build_participant_manifest(task_df: pd.DataFrame, figure_df: pd.DataFrame) -> pd.DataFrame:
    """Build the merged participant-level behavioural/figure manifest.

    Args:
        task_df: One-row-per-task-file table.
        figure_df: One-row-per-top-level-figure-directory table.

    Returns:
        A pandas ``DataFrame`` with one row per expected or observed participant
        ID, task-file fields, figure-directory fields, historical categories,
        observable mismatch flags, and factual inspection issue codes.

    Side effects:
        None.
    """

    task_by_id = aggregate_task_files(task_df)
    figure_by_id = aggregate_figure_directories(figure_df)
    manifest_df = pd.DataFrame({"participant_id": all_manifest_participant_ids(task_df, figure_df)})

    if not task_by_id.empty:
        manifest_df = manifest_df.merge(task_by_id, on="participant_id", how="left")
    if not figure_by_id.empty:
        manifest_df = manifest_df.merge(figure_by_id, on="participant_id", how="left")

    manifest_df["participant_id_padded"] = manifest_df["participant_id"].map(lambda value: f"{int(value):03d}")

    # Fill absent task/figure fields after merging. Missing status here means
    # no parsed local file/directory for the participant ID, not a study-level
    # decision.
    manifest_df["task_file_count"] = manifest_df.get("task_file_count", pd.Series(dtype="float")).fillna(0).astype(int)
    manifest_df["figure_directory_count"] = (
        manifest_df.get("figure_directory_count", pd.Series(dtype="float")).fillna(0).astype(int)
    )
    manifest_df["task_file_status"] = manifest_df["task_file_status"].fillna("missing")
    manifest_df["figure_directory_status"] = manifest_df["figure_directory_status"].fillna("missing")

    string_columns = [
        "task_filename",
        "task_file_path",
        "task_date",
        "task_line_count",
        "task_filename_parse_warning",
        "task_read_status",
        "task_read_error_message",
        "figure_directory_name",
        "figure_directory_path",
        "figure_date",
        "figure_zip_count",
        "figure_directory_kind",
        "figure_directory_parse_warning",
    ]
    for column in string_columns:
        if column not in manifest_df.columns:
            manifest_df[column] = ""
        manifest_df[column] = manifest_df[column].fillna("")

    bool_columns = [
        "task_in_incomplete_folder",
        "task_incomplete_style_filename",
        "incomplete_task_flag",
    ]
    for column in bool_columns:
        if column not in manifest_df.columns:
            manifest_df[column] = False
        manifest_df[column] = manifest_df[column].fillna(False).astype(bool)

    historical_rows = [historical_fields_for_participant(int(pid)) for pid in manifest_df["participant_id"]]
    historical_df = pd.DataFrame(historical_rows)
    manifest_df = pd.concat([manifest_df.reset_index(drop=True), historical_df], axis=1)

    manifest_df["observable_mismatch_flags"] = manifest_df.apply(
        lambda row: ";".join(observable_mismatch_flags(row.to_dict())),
        axis=1,
    )
    manifest_df["inspection_issue_codes"] = manifest_df.apply(
        lambda row: ";".join(participant_issue_codes(row.to_dict())),
        axis=1,
    )

    column_order = [
        "participant_id",
        "participant_id_padded",
        "task_file_status",
        "task_filename",
        "task_file_path",
        "task_date",
        "task_line_count",
        "task_in_incomplete_folder",
        "task_incomplete_style_filename",
        "incomplete_task_flag",
        "figure_directory_status",
        "figure_directory_name",
        "figure_directory_path",
        "figure_date",
        "figure_zip_count",
        "observable_mismatch_flags",
        "historical_behavioural_session_category",
        "historical_eeg_category",
        "special_case_category",
        "inspection_issue_codes",
        "task_file_count",
        "figure_directory_count",
        "task_filename_parse_warning",
        "task_read_status",
        "task_read_error_message",
        "figure_directory_kind",
        "figure_directory_parse_warning",
    ]
    return manifest_df[column_order].sort_values("participant_id").reset_index(drop=True)


def count_issue_codes(issue_series: pd.Series) -> Counter[str]:
    """Count semicolon-separated issue codes in a manifest series.

    Args:
        issue_series: Series containing semicolon-separated issue-code strings.

    Returns:
        A ``Counter`` mapping issue code to occurrence count.

    Side effects:
        None.
    """

    counter: Counter[str] = Counter()
    for issue_text in issue_series.fillna(""):
        if not issue_text:
            continue
        counter.update(part for part in str(issue_text).split(";") if part)
    return counter


def write_summary(
    manifest_df: pd.DataFrame,
    figure_df: pd.DataFrame,
    output_path: Path,
    task_dir: Path,
    figure_dir: Path,
) -> None:
    """Write a short Markdown summary of the behavioural file manifest.

    Args:
        manifest_df: Participant-level behavioural/figure manifest.
        figure_df: One-row-per-top-level-figure-directory table.
        output_path: Destination Markdown path.
        task_dir: Task-data directory summarized by this run.
        figure_dir: Figure-data directory summarized by this run.

    Returns:
        ``None``.

    Side effects:
        Writes ``output_path``.
    """

    issue_counter = count_issue_codes(manifest_df["inspection_issue_codes"])
    mismatch_counter = count_issue_codes(manifest_df["observable_mismatch_flags"])
    figure_issue_counter = count_issue_codes(figure_df["figure_issue_codes"]) if not figure_df.empty else Counter()

    incomplete_rows = manifest_df[manifest_df["incomplete_task_flag"]]
    nonstandard_zip_rows = manifest_df[
        manifest_df["inspection_issue_codes"].str.contains("nonstandard_zip_count", na=False)
    ]
    missing_task_rows = manifest_df[manifest_df["inspection_issue_codes"].str.contains("missing_task_file", na=False)]
    missing_figure_rows = manifest_df[
        manifest_df["inspection_issue_codes"].str.contains("missing_figure_directory", na=False)
    ]
    nonparticipant_figure_rows = figure_df[
        figure_df["figure_directory_kind"].isin(["nonparticipant_directory", "dev_test_id"])
    ] if not figure_df.empty else pd.DataFrame()

    lines = [
        "# Behavioural File Manifest Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Task directory: `{task_dir.as_posix()}`",
        f"Figure directory: `{figure_dir.as_posix()}`",
        f"Participant rows written: {len(manifest_df)}",
        f"Task files represented: {int(manifest_df['task_file_count'].sum())}",
        f"Top-level figure directories represented: {len(figure_df)}",
        "",
        "## Participant-level inspection issue counts",
        "",
    ]

    if issue_counter:
        for issue, count in sorted(issue_counter.items()):
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Observable mismatch flag counts", ""])
    if mismatch_counter:
        for issue, count in sorted(mismatch_counter.items()):
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Figure-directory issue counts", ""])
    if figure_issue_counter:
        for issue, count in sorted(figure_issue_counter.items()):
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Incomplete-task marker rows", ""])
    if incomplete_rows.empty:
        lines.append("- none")
    else:
        for _, row in incomplete_rows.iterrows():
            lines.append(f"- {int(row['participant_id'])}: `{row['task_filename']}`")

    lines.extend(["", "## Nonstandard nonzero figure zip counts", ""])
    if nonstandard_zip_rows.empty:
        lines.append("- none")
    else:
        for _, row in nonstandard_zip_rows.iterrows():
            lines.append(f"- {int(row['participant_id'])}: {row['figure_zip_count']} zips")

    lines.extend(["", "## Missing task file IDs", ""])
    if missing_task_rows.empty:
        lines.append("- none")
    else:
        lines.append("- " + ", ".join(str(int(pid)) for pid in missing_task_rows["participant_id"]))

    lines.extend(["", "## Missing figure directory IDs", ""])
    if missing_figure_rows.empty:
        lines.append("- none")
    else:
        lines.append("- " + ", ".join(str(int(pid)) for pid in missing_figure_rows["participant_id"]))

    lines.extend(["", "## Non-participant or dev/test figure directories", ""])
    if nonparticipant_figure_rows.empty:
        lines.append("- none")
    else:
        for _, row in nonparticipant_figure_rows.iterrows():
            lines.append(
                f"- `{row['figure_directory_name']}`: {row['figure_directory_kind']}, "
                f"{row['figure_zip_count']} zips"
            )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This manifest reports file facts and compact historical categories only.",
            "- Zip files are counted but not opened.",
            "- Task files are counted by line but behavioural outcomes are not parsed.",
            "- Missing files and nonstandard counts are review flags, not analysis decisions.",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(manifest_df: pd.DataFrame, figure_df: pd.DataFrame, output_dir: Path, task_dir: Path, figure_dir: Path) -> None:
    """Write all behavioural manifest outputs.

    Args:
        manifest_df: Participant-level behavioural/figure manifest.
        figure_df: One-row-per-top-level-figure-directory table.
        output_dir: Directory where manifest outputs should be written.
        task_dir: Task-data directory summarized by this run.
        figure_dir: Figure-data directory summarized by this run.

    Returns:
        ``None``.

    Side effects:
        Creates ``output_dir`` if needed and writes three local manifest files.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_df.to_csv(output_dir / "behavioral_file_manifest.csv", index=False)
    figure_df.to_csv(output_dir / "figure_zip_counts.csv", index=False)
    write_summary(
        manifest_df=manifest_df,
        figure_df=figure_df,
        output_path=output_dir / "behavioral_file_manifest_summary.md",
        task_dir=task_dir,
        figure_dir=figure_dir,
    )


def main() -> int:
    """Run the behavioural task and figure manifest workflow.

    Args:
        None.

    Returns:
        ``0`` after outputs are written successfully.

    Side effects:
        Reads local task and figure directory trees, counts task-file lines,
        counts zip files by path, creates ``_Data/behavior/manifest/`` if
        needed, writes CSV and Markdown manifest outputs, and prints a concise
        completion summary.
    """

    repo_root = repo_root_from_script()
    task_dir = repo_root / "_Data" / "task"
    figure_dir = repo_root / "_Data" / "figure"
    output_dir = repo_root / "_Data" / "behavior" / "manifest"

    task_files = discover_task_files(task_dir)
    figure_directories = discover_figure_directories(figure_dir)

    task_df = build_task_file_table(task_files, task_dir, repo_root)
    figure_df = build_figure_directory_table(figure_directories, repo_root)
    manifest_df = build_participant_manifest(task_df, figure_df)

    write_outputs(manifest_df, figure_df, output_dir, task_dir, figure_dir)

    print(f"Wrote {output_dir / 'behavioral_file_manifest.csv'}")
    print(f"Wrote {output_dir / 'figure_zip_counts.csv'}")
    print(f"Wrote {output_dir / 'behavioral_file_manifest_summary.md'}")
    print(f"Task files inspected: {len(task_files)}")
    print(f"Top-level figure directories inspected: {len(figure_directories)}")
    print(f"Participant rows written: {len(manifest_df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
