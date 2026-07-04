"""Dry-run validation for future continuous MNE preprocessing.

This script is the first continuous-raw preprocessing entry point for the new
DEMI EEG reanalysis path, but it deliberately does not preprocess signals. It
performs a dry run: it opens each raw EDF with MNE using ``preload=False``,
checks header-level facts, validates channel labels and channel-type mapping,
checks EEG channel names against the tracked BESA-81 coordinate table after
case normalization, and carries forward local event-audit context as QC facts.

A dry-run stage comes before actual continuous preprocessing because the active
MNE workflow is restarting from raw EDF recordings, while the final event and
epoch policy is still pending. Continuous header, channel, coordinate, file
role, annotation-count, and event-context checks can be made now without
filtering data, repairing triggers, constructing epochs, or making
participant-level analytic decisions. This script therefore creates a review
surface for Tony before any later script writes cleaned continuous FIF files.

The script reads these required local inputs:

- ``_Data/eeg/manifest/raw_eeg_file_manifest.csv``;
- ``_Data/eeg/manifest/raw_eeg_annotation_counts.csv``;
- ``_Data/eeg/BESA-81.csv``;
- raw EDF files discovered directly below ``_Data/eeg/raw/``.

When present, it also reads these local event-audit context tables:

- ``_Data/eeg/event_sequence_audit/proposed_offset_join_audit.csv``;
- ``_Data/eeg/event_special_case_audit/split_file_continuity_audit.csv``;
- ``_Data/eeg/event_special_case_audit/id5_concatenated_file_audit.csv``;
- ``_Data/eeg/event_special_case_audit/zero_no_offset_case_audit.csv``;
- ``_Data/eeg/event_special_case_audit/raw_only_row_classification.csv``.

Generated outputs are local-only dry-run QC files under
``_Data/eeg/mne_preprocessing/dry_run/``:

- ``preprocessing_dry_run_file_summary.csv``;
- ``channel_type_validation.csv``;
- ``montage_coordinate_validation.csv``;
- ``event_context_flags.csv``;
- ``preprocessing_dry_run_summary.md``.

What this script explicitly does not do:

- no EEG filtering;
- no ICA;
- no bad-channel interpolation;
- no CSD computation;
- no epoch construction;
- no cleaned FIF or EDF derivative writing;
- no event repair;
- no participant-retention decision.

Safety boundaries:

- EDFs are opened with ``mne.io.read_raw_edf(..., preload=False)``.
- Raw signal arrays are not loaded.
- Raw EDF files are not modified.
- One missing or problematic EDF records a factual row and does not stop the
  dry run.
- Optional event-audit tables are treated as context only. Their facts are not
  converted into event policy.
- All generated files are below ignored ``_Data/`` paths.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/04_preprocess_continuous_raw.py

Paths are resolved relative to this script, so the command also works when
launched from another current working directory.
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

# MNE may import matplotlib depending on the local environment. Keep any cache
# writes inside a temp directory so this header-only dry run does not depend on
# the user's global matplotlib configuration path.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "demi_matplotlib"))

import mne
import pandas as pd


RAW_EEG_DIR = Path("_Data") / "eeg" / "raw"
RAW_MANIFEST_PATH = Path("_Data") / "eeg" / "manifest" / "raw_eeg_file_manifest.csv"
ANNOTATION_COUNTS_PATH = Path("_Data") / "eeg" / "manifest" / "raw_eeg_annotation_counts.csv"
BESA_COORDINATES_PATH = Path("_Data") / "eeg" / "BESA-81.csv"
DRY_RUN_OUTPUT_DIR = Path("_Data") / "eeg" / "mne_preprocessing" / "dry_run"

OPTIONAL_EVENT_CONTEXT_PATHS = {
    "proposed_offset_join_audit": Path("_Data")
    / "eeg"
    / "event_sequence_audit"
    / "proposed_offset_join_audit.csv",
    "split_file_continuity_audit": Path("_Data")
    / "eeg"
    / "event_special_case_audit"
    / "split_file_continuity_audit.csv",
    "id5_concatenated_file_audit": Path("_Data")
    / "eeg"
    / "event_special_case_audit"
    / "id5_concatenated_file_audit.csv",
    "zero_no_offset_case_audit": Path("_Data")
    / "eeg"
    / "event_special_case_audit"
    / "zero_no_offset_case_audit.csv",
    "raw_only_row_classification": Path("_Data")
    / "eeg"
    / "event_special_case_audit"
    / "raw_only_row_classification.csv",
}

FILE_SUMMARY_FILENAME = "preprocessing_dry_run_file_summary.csv"
CHANNEL_TYPE_VALIDATION_FILENAME = "channel_type_validation.csv"
MONTAGE_COORDINATE_VALIDATION_FILENAME = "montage_coordinate_validation.csv"
EVENT_CONTEXT_FLAGS_FILENAME = "event_context_flags.csv"
DRY_RUN_SUMMARY_FILENAME = "preprocessing_dry_run_summary.md"

EXPECTED_SAMPLING_FREQUENCY_HZ = 1000.0
SAMPLING_FREQUENCY_TOLERANCE_HZ = 0.01
EXPECTED_CHANNEL_COUNT = 37
LOW_ANNOTATION_COUNT_THRESHOLD = 600

EXPECTED_EEG_CHANNELS = (
    "FP1",
    "FP2",
    "F7",
    "F3",
    "FZ",
    "F4",
    "F8",
    "FT7",
    "FC3",
    "FCZ",
    "FC4",
    "FT8",
    "T7",
    "C3",
    "CZ",
    "C4",
    "T8",
    "M1",
    "TP7",
    "CP3",
    "CPZ",
    "CP4",
    "TP8",
    "M2",
    "P7",
    "P3",
    "PZ",
    "P4",
    "P8",
    "O1",
    "OZ",
    "O2",
)
EXPECTED_EOG_CHANNELS = ("HEO", "VEO")
EXPECTED_EMG_CHANNELS = ("EMG-L", "EMG-A")
EXPECTED_TRIGGER_CHANNELS = ("TRIGGER",)
EXPECTED_ALL_CHANNELS = (
    EXPECTED_EEG_CHANNELS + EXPECTED_EOG_CHANNELS + EXPECTED_EMG_CHANNELS + EXPECTED_TRIGGER_CHANNELS
)

EXPECTED_CHANNEL_TYPE_BY_NAME = {
    **{channel: "eeg" for channel in EXPECTED_EEG_CHANNELS},
    **{channel: "eog" for channel in EXPECTED_EOG_CHANNELS},
    **{channel: "emg" for channel in EXPECTED_EMG_CHANNELS},
    **{channel: "stim" for channel in EXPECTED_TRIGGER_CHANNELS},
}

KNOWN_FILE_ROLES = {"single", "split_part", "concatenated", "other"}
SPLIT_FILE_IDS = {54, 56, 65}
CONCATENATED_FILE_IDS = {5}
ZERO_NO_OFFSET_CONTEXT_IDS = {11, 14, 89, 94, 100}

EDF_FILENAME_RE = re.compile(
    r"^demi_(?P<participant_id>\d{1,3})(?:_(?P<split_part>\d+))?(?P<label>.*)\.edf$",
    flags=re.IGNORECASE,
)

FILE_SUMMARY_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_path",
    "manifest_present",
    "discovered_in_raw_dir",
    "file_exists",
    "file_role",
    "split_part",
    "read_status",
    "mne_read_error_type",
    "mne_read_error_message",
    "sampling_frequency_hz",
    "sampling_frequency_approximately_1000_hz",
    "duration_seconds",
    "n_channels",
    "expected_37_channels_present",
    "expected_eeg_labels_present",
    "expected_eog_labels_present",
    "expected_emg_labels_present",
    "expected_trigger_label_present",
    "unexpected_channel_names",
    "missing_expected_channel_names",
    "channel_order_matches_expected",
    "channel_type_assignment_complete",
    "unmapped_channel_type_names",
    "eeg_channels_matched_to_besa",
    "besa_unmatched_eeg_channel_names",
    "manifest_annotation_count",
    "annotation_count_from_counts_csv",
    "mne_annotation_count",
    "annotation_count_sources_agree",
    "zero_annotation_fact",
    "low_annotation_fact",
    "file_role_context",
    "event_context_pending_codes",
    "dry_run_issue_codes",
    "no_preprocessing_performed",
]

CHANNEL_TYPE_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_role",
    "split_part",
    "read_status",
    "channel_index",
    "raw_channel_name",
    "normalized_channel_name",
    "mne_reported_type",
    "expected_dry_run_type",
    "channel_type_assignment_status",
    "channel_group",
    "channel_type_issue_codes",
]

MONTAGE_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_role",
    "split_part",
    "read_status",
    "raw_channel_name",
    "normalized_channel_name",
    "raw_eeg_label_present",
    "expected_eeg_channel",
    "besa_coordinate_match",
    "besa_channel_name",
    "besa_x",
    "besa_y",
    "besa_z",
    "montage_issue_codes",
]

EVENT_CONTEXT_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_role",
    "split_part",
    "annotation_count",
    "zero_annotation_fact",
    "low_annotation_fact",
    "split_file_context_present",
    "split_file_context_row_kinds",
    "split_file_strict_clean_timing_rows",
    "split_file_raw_only_rows",
    "split_file_offset_only_rows",
    "concatenated_context_present",
    "id5_file_start_onset_seconds",
    "id5_segments_with_strict_clean_timing",
    "zero_no_offset_context_present",
    "zero_no_offset_factual_note",
    "old_compatible_offset_rows",
    "proposed_join_status_counts",
    "raw_only_context_rows",
    "raw_only_category_counts",
    "event_context_pending_codes",
    "event_context_source_status",
]


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


def text_or_empty(value: Any) -> str:
    """Return stripped text while preserving missing values as empty text.

    Args:
        value: Scalar value from a DataFrame cell or Python object.

    Returns:
        Empty string for missing values, otherwise stripped string text.

    Side effects:
        None.
    """

    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_channel_name(value: Any) -> str:
    """Normalize a channel label for conservative case-insensitive matching.

    Args:
        value: Channel label from an EDF header, manifest cell, or coordinate
            table.

    Returns:
        Uppercase stripped channel label. The function intentionally does not
        remove hyphens or otherwise rewrite labels, because the current dry run
        is checking whether simple case normalization is enough.

    Side effects:
        None.
    """

    return text_or_empty(value).upper()


def semicolon_join(values: list[Any] | pd.Series) -> str:
    """Join unique non-empty values into deterministic semicolon text.

    Args:
        values: Values to convert to text.

    Returns:
        Semicolon-delimited text in first-seen order.

    Side effects:
        None.
    """

    seen: dict[str, None] = {}
    for value in values:
        text = text_or_empty(value)
        if text:
            seen[text] = None
    return ";".join(seen.keys())


def value_counts_text(values: pd.Series) -> str:
    """Format non-empty value counts as ``label=count`` text.

    Args:
        values: Series of values to count.

    Returns:
        Semicolon-delimited count text, or empty text when no value is present.

    Side effects:
        None.
    """

    cleaned = values.dropna().map(text_or_empty)
    cleaned = cleaned[cleaned.ne("")]
    if cleaned.empty:
        return ""
    counts = cleaned.value_counts(sort=False)
    return ";".join(f"{label}={int(count)}" for label, count in counts.items())


def safe_int(value: Any) -> int | None:
    """Convert a value to ``int`` when possible.

    Args:
        value: Scalar value that may be numeric, text, or missing.

    Returns:
        Integer value, or ``None`` when conversion is not possible.

    Side effects:
        None.
    """

    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return None
    return int(numeric)


def safe_float(value: Any) -> float | None:
    """Convert a value to ``float`` when possible.

    Args:
        value: Scalar value that may be numeric, text, or missing.

    Returns:
        Floating-point value, or ``None`` when conversion is not possible.

    Side effects:
        None.
    """

    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return None
    return float(numeric)


def parse_edf_filename(edf_path: Path) -> dict[str, Any]:
    """Parse participant ID and file-role hints from an EDF filename.

    Args:
        edf_path: EDF path. Only the filename is parsed.

    Returns:
        Dictionary with participant ID, padded ID, normalized file role,
        split-part number, and filename warning text. Unparsed filenames are
        assigned the ``other`` file role.

    Side effects:
        None.
    """

    match = EDF_FILENAME_RE.match(edf_path.name)
    if match is None:
        return {
            "participant_id": None,
            "participant_id_padded": "",
            "file_role": "other",
            "split_part": None,
            "filename_parse_warning": "filename_does_not_match_expected_demi_pattern",
        }

    participant_id = int(match.group("participant_id"))
    split_part_text = match.group("split_part")
    split_part = int(split_part_text) if split_part_text is not None else None
    label_text = match.group("label").lower()

    if split_part is not None:
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


def normalize_file_role(value: Any) -> str:
    """Convert a manifest or filename role into this script's role vocabulary.

    Args:
        value: Role value from the raw manifest or filename parser.

    Returns:
        One of ``single``, ``split_part``, ``concatenated``, or ``other``.

    Side effects:
        None.
    """

    role = text_or_empty(value).lower()
    if role in KNOWN_FILE_ROLES:
        return role
    if role in {"split", "split part", "split-part"}:
        return "split_part"
    if not role or role == "unknown":
        return "other"
    return "other"


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
        return path.resolve().relative_to(repo_root).as_posix()
    except ValueError:
        return path.as_posix()


def read_csv_checked(path: Path, required_columns: set[str], label: str) -> pd.DataFrame:
    """Read a required CSV and verify expected columns.

    Args:
        path: CSV path.
        required_columns: Columns required by this dry run.
        label: Human-readable input label for failure messages.

    Returns:
        DataFrame read from ``path``.

    Side effects:
        Reads a local CSV file. Raises ``RuntimeError`` when the file or a
        required column is missing.
    """

    if not path.exists():
        raise RuntimeError(f"required {label} is missing: {path.as_posix()}")

    frame = pd.read_csv(path, low_memory=False)
    missing = sorted(required_columns.difference(frame.columns))
    if missing:
        raise RuntimeError(f"{label} is missing required column(s): {', '.join(missing)}")
    return frame


def read_optional_event_context(repo_root: Path) -> dict[str, pd.DataFrame | None]:
    """Read optional event-audit context tables when they are present.

    Args:
        repo_root: Absolute repository root path.

    Returns:
        Dictionary keyed by context-table label. Values are DataFrames for
        present tables and ``None`` for absent tables.

    Side effects:
        Reads local CSV files below ``_Data/`` when they exist.
    """

    context: dict[str, pd.DataFrame | None] = {}
    for label, relative_path in OPTIONAL_EVENT_CONTEXT_PATHS.items():
        path = repo_root / relative_path
        context[label] = pd.read_csv(path, low_memory=False) if path.exists() else None
    return context


def optional_context_source_status(context: dict[str, pd.DataFrame | None]) -> str:
    """Summarize which optional event-context tables were present.

    Args:
        context: Optional context dictionary returned by
            ``read_optional_event_context``.

    Returns:
        Semicolon-delimited ``label=present`` or ``label=missing`` text.

    Side effects:
        None.
    """

    return ";".join(f"{label}={'present' if frame is not None else 'missing'}" for label, frame in context.items())


def discover_edf_files(raw_dir: Path) -> list[Path]:
    """Find EDF files directly below the raw EEG directory.

    Args:
        raw_dir: Directory expected to contain DEMI raw EDF files.

    Returns:
        Sorted list of EDF paths. The search is case-insensitive for the EDF
        suffix and is not recursive.

    Side effects:
        Reads directory entries only.
    """

    if not raw_dir.exists():
        return []
    return sorted(path for path in raw_dir.iterdir() if path.is_file() and path.suffix.lower() == ".edf")


def load_besa_coordinates(path: Path) -> pd.DataFrame:
    """Read BESA-81 coordinates with normalized channel labels.

    Args:
        path: Path to ``_Data/eeg/BESA-81.csv``.

    Returns:
        Coordinate DataFrame with ``normalized_channel_name`` added.

    Side effects:
        Reads the tracked BESA coordinate CSV.
    """

    coordinates = read_csv_checked(path, {"chan", "x", "y", "z"}, "BESA-81 coordinate CSV")
    coordinates = coordinates.copy()
    coordinates["normalized_channel_name"] = coordinates["chan"].map(normalize_channel_name)
    return coordinates


def aggregate_annotation_counts(annotation_counts: pd.DataFrame) -> pd.DataFrame:
    """Reduce the long annotation-count CSV to one row per EDF file.

    Args:
        annotation_counts: ``raw_eeg_annotation_counts.csv``.

    Returns:
        DataFrame with ``source_filename`` and total annotation count from the
        long count table.

    Side effects:
        None.
    """

    totals = (
        annotation_counts.groupby("source_filename", dropna=False, as_index=False)["annotation_count"]
        .sum()
        .rename(columns={"annotation_count": "annotation_count_from_counts_csv"})
    )
    return totals


def manifest_path_for_row(row: pd.Series, raw_dir: Path, repo_root: Path) -> Path:
    """Resolve the EDF path represented by a raw-manifest row.

    Args:
        row: Row from ``raw_eeg_file_manifest.csv``.
        raw_dir: Raw EDF directory.
        repo_root: Repository root.

    Returns:
        Best local path for the EDF. The manifest ``file_path`` is preferred
        when available; otherwise ``raw_dir / source_filename`` is used.

    Side effects:
        None.
    """

    file_path_text = text_or_empty(row.get("file_path"))
    if file_path_text:
        path = Path(file_path_text)
        return path if path.is_absolute() else repo_root / path
    return raw_dir / text_or_empty(row.get("source_filename"))


def build_file_inventory(
    manifest: pd.DataFrame,
    annotation_totals: pd.DataFrame,
    edf_files: list[Path],
    raw_dir: Path,
    repo_root: Path,
) -> pd.DataFrame:
    """Build the union of manifest rows and discovered raw EDF files.

    Args:
        manifest: Raw EEG file manifest.
        annotation_totals: Per-file totals from the long annotation-count CSV.
        edf_files: EDF files discovered below the raw directory.
        raw_dir: Raw EDF directory.
        repo_root: Repository root.

    Returns:
        One row per manifest or discovered EDF filename. Rows track whether the
        file was present in the manifest, discovered on disk, and currently
        exists at the resolved local path.

    Side effects:
        Checks file existence.
    """

    manifest_by_name = manifest.drop_duplicates("source_filename", keep="first").set_index("source_filename")
    totals_by_name = annotation_totals.set_index("source_filename")
    discovered_by_name = {path.name: path for path in edf_files}
    all_names = sorted(set(manifest_by_name.index).union(discovered_by_name))

    rows: list[dict[str, Any]] = []
    for source_filename in all_names:
        manifest_present = source_filename in manifest_by_name.index
        discovered_in_raw_dir = source_filename in discovered_by_name
        manifest_row = manifest_by_name.loc[source_filename] if manifest_present else pd.Series(dtype=object)
        parsed = parse_edf_filename(Path(source_filename))

        edf_path = discovered_by_name[source_filename] if discovered_in_raw_dir else manifest_path_for_row(
            manifest_row, raw_dir, repo_root
        )
        manifest_role = normalize_file_role(manifest_row.get("file_role")) if manifest_present else "other"
        parsed_role = normalize_file_role(parsed["file_role"])
        file_role = manifest_role if manifest_present and manifest_role != "other" else parsed_role

        participant_id = safe_int(manifest_row.get("participant_id")) if manifest_present else parsed["participant_id"]
        split_part = safe_int(manifest_row.get("split_part")) if manifest_present else parsed["split_part"]
        manifest_annotation_count = (
            safe_int(manifest_row.get("annotation_count")) if manifest_present else None
        )
        annotation_count_from_counts_csv = (
            safe_int(totals_by_name.loc[source_filename, "annotation_count_from_counts_csv"])
            if source_filename in totals_by_name.index
            else None
        )

        rows.append(
            {
                "participant_id": participant_id,
                "participant_id_padded": f"{participant_id:03d}" if participant_id is not None else "",
                "source_filename": source_filename,
                "file_path": relative_to_repo(edf_path, repo_root),
                "manifest_present": manifest_present,
                "discovered_in_raw_dir": discovered_in_raw_dir,
                "file_exists": edf_path.exists(),
                "resolved_edf_path": edf_path,
                "file_role": file_role,
                "split_part": split_part,
                "filename_parse_warning": parsed["filename_parse_warning"],
                "manifest_annotation_count": manifest_annotation_count,
                "annotation_count_from_counts_csv": annotation_count_from_counts_csv,
            }
        )

    inventory = pd.DataFrame(rows)
    if inventory.empty:
        return pd.DataFrame(columns=["source_filename"])
    return inventory.sort_values(
        ["participant_id", "split_part", "source_filename"], kind="mergesort", na_position="last"
    ).reset_index(drop=True)


def expected_channel_type(normalized_channel_name: str) -> str:
    """Return the dry-run channel type expected for a normalized channel name.

    Args:
        normalized_channel_name: Uppercase channel label.

    Returns:
        ``eeg``, ``eog``, ``emg``, ``stim``, or ``unmapped``.

    Side effects:
        None.
    """

    return EXPECTED_CHANNEL_TYPE_BY_NAME.get(normalized_channel_name, "unmapped")


def channel_group_for_type(channel_type: str) -> str:
    """Return a broad channel group label for dry-run output.

    Args:
        channel_type: Expected dry-run channel type.

    Returns:
        Group label used in the channel-type validation table.

    Side effects:
        None.
    """

    if channel_type == "eeg":
        return "scalp_eeg"
    if channel_type in {"eog", "emg"}:
        return "auxiliary"
    if channel_type == "stim":
        return "trigger"
    return "unmapped"


def inspect_edf_header(edf_path: Path) -> dict[str, Any]:
    """Read EDF header metadata through MNE without loading signal data.

    Args:
        edf_path: Path to one raw EDF file.

    Returns:
        Dictionary with read status, sampling frequency, duration, channel
        labels, MNE-reported channel types, and annotation count. Error rows
        keep the exception type and message.

    Side effects:
        Opens the EDF through MNE with ``preload=False``. The Raw object is
        closed before returning when MNE exposes ``close()``.
    """

    raw: mne.io.BaseRaw | None = None
    try:
        raw = mne.io.read_raw_edf(edf_path, preload=False, verbose="ERROR")
        sfreq = float(raw.info["sfreq"])
        duration_seconds = float(raw.n_times / sfreq) if sfreq else None
        channel_names = list(raw.ch_names)
        mne_channel_types = list(raw.get_channel_types())
        annotation_count = int(len(raw.annotations))
        annotation_description_counts = Counter(str(description) for description in raw.annotations.description)

        return {
            "read_status": "ok",
            "mne_read_error_type": "",
            "mne_read_error_message": "",
            "sampling_frequency_hz": sfreq,
            "duration_seconds": duration_seconds,
            "channel_names": channel_names,
            "mne_channel_types": mne_channel_types,
            "mne_annotation_count": annotation_count,
            "mne_annotation_description_counts_json": json.dumps(
                dict(sorted(annotation_description_counts.items())), sort_keys=True
            ),
        }
    except Exception as error:  # noqa: BLE001 - one problematic EDF should not stop the dry run.
        return {
            "read_status": "mne_read_error",
            "mne_read_error_type": type(error).__name__,
            "mne_read_error_message": str(error),
            "sampling_frequency_hz": None,
            "duration_seconds": None,
            "channel_names": [],
            "mne_channel_types": [],
            "mne_annotation_count": None,
            "mne_annotation_description_counts_json": "{}",
        }
    finally:
        if raw is not None and hasattr(raw, "close"):
            raw.close()


def build_channel_type_rows(file_row: dict[str, Any], header: dict[str, Any]) -> list[dict[str, Any]]:
    """Build one channel-type validation row per observed EDF channel.

    Args:
        file_row: File inventory row as a dictionary.
        header: Header metadata returned by ``inspect_edf_header``.

    Returns:
        List of channel validation rows for the output CSV. Read failures and
        missing files receive a single factual row with blank channel fields.

    Side effects:
        None.
    """

    if header["read_status"] != "ok":
        return [
            {
                **file_identity_fields(file_row),
                "read_status": header["read_status"],
                "channel_index": pd.NA,
                "raw_channel_name": "",
                "normalized_channel_name": "",
                "mne_reported_type": "",
                "expected_dry_run_type": "",
                "channel_type_assignment_status": "file_not_read",
                "channel_group": "",
                "channel_type_issue_codes": header["read_status"],
            }
        ]

    rows: list[dict[str, Any]] = []
    for index, (channel_name, mne_type) in enumerate(
        zip(header["channel_names"], header["mne_channel_types"]), start=1
    ):
        normalized_name = normalize_channel_name(channel_name)
        dry_run_type = expected_channel_type(normalized_name)

        # Channel typing is intentionally name-based here. EDF channel metadata
        # can be sparse, while the DEMI acquisition labels have stable names for
        # scalp EEG, eye, muscle, and trigger channels. This dry run checks that
        # a future in-memory ``set_channel_types`` mapping would be unambiguous.
        if dry_run_type == "unmapped":
            status = "unmapped_channel_label"
            issue_codes = "unmapped_channel_label"
        elif dry_run_type == mne_type:
            status = "matches_mne_reported_type"
            issue_codes = ""
        else:
            status = "dry_run_type_override_needed"
            issue_codes = "mne_reported_type_differs_from_dry_run_type"

        rows.append(
            {
                **file_identity_fields(file_row),
                "read_status": header["read_status"],
                "channel_index": index,
                "raw_channel_name": channel_name,
                "normalized_channel_name": normalized_name,
                "mne_reported_type": mne_type,
                "expected_dry_run_type": dry_run_type,
                "channel_type_assignment_status": status,
                "channel_group": channel_group_for_type(dry_run_type),
                "channel_type_issue_codes": issue_codes,
            }
        )

    return rows


def build_montage_coordinate_rows(
    file_row: dict[str, Any],
    header: dict[str, Any],
    besa_by_normalized_name: dict[str, pd.Series],
) -> list[dict[str, Any]]:
    """Build EEG channel-to-BESA coordinate validation rows.

    Args:
        file_row: File inventory row as a dictionary.
        header: Header metadata returned by ``inspect_edf_header``.
        besa_by_normalized_name: Mapping from normalized BESA labels to
            coordinate rows.

    Returns:
        One row per expected EEG channel, plus observed EEG rows when an
        unexpected scalp label appears. Read failures and missing files receive
        a single factual row with blank channel fields.

    Side effects:
        None.
    """

    if header["read_status"] != "ok":
        return [
            {
                **file_identity_fields(file_row),
                "read_status": header["read_status"],
                "raw_channel_name": "",
                "normalized_channel_name": "",
                "raw_eeg_label_present": False,
                "expected_eeg_channel": False,
                "besa_coordinate_match": False,
                "besa_channel_name": "",
                "besa_x": pd.NA,
                "besa_y": pd.NA,
                "besa_z": pd.NA,
                "montage_issue_codes": header["read_status"],
            }
        ]

    observed_by_normalized_name = {
        normalize_channel_name(channel_name): channel_name for channel_name in header["channel_names"]
    }
    expected_eeg_names = list(EXPECTED_EEG_CHANNELS)
    observed_eeg_names = [
        normalize_channel_name(channel_name)
        for channel_name in header["channel_names"]
        if expected_channel_type(normalize_channel_name(channel_name)) == "eeg"
    ]
    montage_names = sorted(set(expected_eeg_names).union(observed_eeg_names))

    rows: list[dict[str, Any]] = []
    for normalized_name in montage_names:
        raw_channel_name = observed_by_normalized_name.get(normalized_name, "")
        raw_present = bool(raw_channel_name)
        expected_eeg_channel = normalized_name in EXPECTED_EEG_CHANNELS
        besa_row = besa_by_normalized_name.get(normalized_name)
        besa_match = besa_row is not None

        # The BESA file uses mixed-case 10-10 style names (for example Fp1 and
        # Fz), while the EDF headers often use uppercase labels. Only case is
        # normalized so that a failed match remains visible as a real label or
        # coordinate mismatch, not as an implicit rename.
        issues: list[str] = []
        if expected_eeg_channel and not raw_present:
            issues.append("expected_eeg_channel_absent_from_raw_header")
        if raw_present and not besa_match:
            issues.append("raw_eeg_channel_without_besa_coordinate")
        if expected_eeg_channel and not besa_match:
            issues.append("expected_eeg_channel_without_besa_coordinate")

        rows.append(
            {
                **file_identity_fields(file_row),
                "read_status": header["read_status"],
                "raw_channel_name": raw_channel_name,
                "normalized_channel_name": normalized_name,
                "raw_eeg_label_present": raw_present,
                "expected_eeg_channel": expected_eeg_channel,
                "besa_coordinate_match": besa_match,
                "besa_channel_name": text_or_empty(besa_row["chan"]) if besa_match else "",
                "besa_x": safe_float(besa_row["x"]) if besa_match else pd.NA,
                "besa_y": safe_float(besa_row["y"]) if besa_match else pd.NA,
                "besa_z": safe_float(besa_row["z"]) if besa_match else pd.NA,
                "montage_issue_codes": ";".join(issues),
            }
        )

    return rows


def file_identity_fields(file_row: dict[str, Any]) -> dict[str, Any]:
    """Return common participant/file identity fields for output tables.

    Args:
        file_row: File inventory row as a dictionary.

    Returns:
        Dictionary with participant ID, padded ID, source filename, file role,
        and split-part number.

    Side effects:
        None.
    """

    return {
        "participant_id": file_row.get("participant_id"),
        "participant_id_padded": file_row.get("participant_id_padded", ""),
        "source_filename": file_row.get("source_filename", ""),
        "file_role": file_row.get("file_role", "other"),
        "split_part": file_row.get("split_part"),
    }


def expected_channel_validation(channel_names: list[str]) -> dict[str, Any]:
    """Validate observed channel labels against the expected DEMI labels.

    Args:
        channel_names: Channel labels returned by MNE for one EDF.

    Returns:
        Dictionary with label-presence booleans and semicolon-delimited missing
        or unexpected channel names.

    Side effects:
        None.
    """

    normalized_names = [normalize_channel_name(channel_name) for channel_name in channel_names]
    observed = set(normalized_names)
    expected = set(EXPECTED_ALL_CHANNELS)
    missing = [channel for channel in EXPECTED_ALL_CHANNELS if channel not in observed]
    unexpected = [channel for channel in normalized_names if channel not in expected]

    return {
        "expected_37_channels_present": len(channel_names) == EXPECTED_CHANNEL_COUNT,
        "expected_eeg_labels_present": all(channel in observed for channel in EXPECTED_EEG_CHANNELS),
        "expected_eog_labels_present": all(channel in observed for channel in EXPECTED_EOG_CHANNELS),
        "expected_emg_labels_present": all(channel in observed for channel in EXPECTED_EMG_CHANNELS),
        "expected_trigger_label_present": all(channel in observed for channel in EXPECTED_TRIGGER_CHANNELS),
        "unexpected_channel_names": semicolon_join(unexpected),
        "missing_expected_channel_names": semicolon_join(missing),
        "channel_order_matches_expected": tuple(normalized_names) == EXPECTED_ALL_CHANNELS,
    }


def annotation_counts_agree(file_row: dict[str, Any], mne_annotation_count: int | None) -> bool:
    """Check whether available annotation-count sources agree.

    Args:
        file_row: File inventory row as a dictionary.
        mne_annotation_count: Annotation count read from the EDF header, or
            ``None`` when the EDF could not be read.

    Returns:
        ``True`` when all available count sources are identical or only one
        source is available.

    Side effects:
        None.
    """

    values = [
        file_row.get("manifest_annotation_count"),
        file_row.get("annotation_count_from_counts_csv"),
        mne_annotation_count,
    ]
    present_values = [int(value) for value in values if value is not None and not pd.isna(value)]
    if len(present_values) <= 1:
        return True
    return len(set(present_values)) == 1


def pending_event_context_codes(participant_id: int | None, file_role: str, annotation_count: int | None) -> list[str]:
    """Return event-policy context codes that remain factual in this dry run.

    Args:
        participant_id: Parsed DEMI participant ID, or ``None``.
        file_role: Normalized file role.
        annotation_count: Annotation count from the EDF header when available.

    Returns:
        Short context codes. These are reminders that event policy is pending;
        they are not preprocessing or epoch-construction actions.

    Side effects:
        None.
    """

    codes: list[str] = []
    if file_role == "split_part" or participant_id in SPLIT_FILE_IDS:
        codes.append("event_policy_pending_split_file_continuity")
    if file_role == "concatenated" or participant_id in CONCATENATED_FILE_IDS:
        codes.append("event_policy_pending_concatenated_file_start")
    if participant_id in ZERO_NO_OFFSET_CONTEXT_IDS:
        codes.append("event_policy_pending_zero_or_no_offset_context")
    if annotation_count == 0:
        codes.append("zero_annotation_fact_only")
    elif annotation_count is not None and annotation_count < LOW_ANNOTATION_COUNT_THRESHOLD:
        codes.append("low_annotation_fact_only")
    return codes


def build_file_summary_row(file_row: dict[str, Any], header: dict[str, Any], montage_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Build one file-level dry-run summary row.

    Args:
        file_row: File inventory row as a dictionary.
        header: Header metadata returned by ``inspect_edf_header``.
        montage_rows: Montage validation rows for the same file.

    Returns:
        Dictionary for ``preprocessing_dry_run_file_summary.csv``.

    Side effects:
        None.
    """

    channel_validation = expected_channel_validation(header["channel_names"])
    mne_annotation_count = header["mne_annotation_count"]
    available_annotation_count = (
        mne_annotation_count
        if mne_annotation_count is not None
        else file_row.get("manifest_annotation_count", file_row.get("annotation_count_from_counts_csv"))
    )
    zero_annotation_fact = available_annotation_count == 0
    low_annotation_fact = (
        available_annotation_count is not None
        and available_annotation_count < LOW_ANNOTATION_COUNT_THRESHOLD
        and not zero_annotation_fact
    )

    sampling_frequency = header["sampling_frequency_hz"]
    sampling_frequency_ok = (
        sampling_frequency is not None
        and abs(float(sampling_frequency) - EXPECTED_SAMPLING_FREQUENCY_HZ) <= SAMPLING_FREQUENCY_TOLERANCE_HZ
    )
    channel_type_rows = build_channel_type_rows(file_row, header)
    unmapped_channel_names = [
        row["raw_channel_name"]
        for row in channel_type_rows
        if row.get("channel_type_assignment_status") == "unmapped_channel_label"
    ]
    montage_problem_rows = [
        row for row in montage_rows if text_or_empty(row.get("montage_issue_codes")) and row.get("read_status") == "ok"
    ]
    unmatched_besa_names = [
        row["normalized_channel_name"]
        for row in montage_problem_rows
        if "without_besa_coordinate" in text_or_empty(row.get("montage_issue_codes"))
    ]

    issue_codes: list[str] = []
    if not file_row["manifest_present"]:
        issue_codes.append("file_not_in_raw_manifest")
    if not file_row["discovered_in_raw_dir"]:
        issue_codes.append("manifest_file_not_discovered_in_raw_dir")
    if not file_row["file_exists"]:
        issue_codes.append("raw_file_missing")
    if header["read_status"] != "ok":
        issue_codes.append(header["read_status"])
    if header["read_status"] == "ok" and not sampling_frequency_ok:
        issue_codes.append("sampling_frequency_not_approximately_1000_hz")
    if header["read_status"] == "ok" and not channel_validation["expected_37_channels_present"]:
        issue_codes.append("unexpected_channel_count")
    if header["read_status"] == "ok" and channel_validation["missing_expected_channel_names"]:
        issue_codes.append("missing_expected_channel_labels")
    if header["read_status"] == "ok" and channel_validation["unexpected_channel_names"]:
        issue_codes.append("unexpected_channel_labels")
    if unmapped_channel_names:
        issue_codes.append("channel_type_assignment_unmapped_label")
    if unmatched_besa_names:
        issue_codes.append("besa_coordinate_unmatched_eeg_label")
    if not annotation_counts_agree(file_row, mne_annotation_count):
        issue_codes.append("annotation_count_sources_disagree")
    if zero_annotation_fact:
        issue_codes.append("zero_annotations")
    elif low_annotation_fact:
        issue_codes.append("low_annotation_count")
    if file_row["file_role"] == "split_part":
        issue_codes.append("split_part_file")
    elif file_row["file_role"] == "concatenated":
        issue_codes.append("concatenated_file")
    elif file_row["file_role"] == "other":
        issue_codes.append("other_file_role")

    event_codes = pending_event_context_codes(file_row["participant_id"], file_row["file_role"], available_annotation_count)

    return {
        "participant_id": file_row["participant_id"],
        "participant_id_padded": file_row["participant_id_padded"],
        "source_filename": file_row["source_filename"],
        "file_path": file_row["file_path"],
        "manifest_present": file_row["manifest_present"],
        "discovered_in_raw_dir": file_row["discovered_in_raw_dir"],
        "file_exists": file_row["file_exists"],
        "file_role": file_row["file_role"],
        "split_part": file_row["split_part"],
        "read_status": header["read_status"],
        "mne_read_error_type": header["mne_read_error_type"],
        "mne_read_error_message": header["mne_read_error_message"],
        "sampling_frequency_hz": sampling_frequency,
        "sampling_frequency_approximately_1000_hz": sampling_frequency_ok,
        "duration_seconds": header["duration_seconds"],
        "n_channels": len(header["channel_names"]) if header["read_status"] == "ok" else pd.NA,
        **channel_validation,
        "channel_type_assignment_complete": not unmapped_channel_names and header["read_status"] == "ok",
        "unmapped_channel_type_names": semicolon_join(unmapped_channel_names),
        "eeg_channels_matched_to_besa": not unmatched_besa_names and header["read_status"] == "ok",
        "besa_unmatched_eeg_channel_names": semicolon_join(unmatched_besa_names),
        "manifest_annotation_count": file_row["manifest_annotation_count"],
        "annotation_count_from_counts_csv": file_row["annotation_count_from_counts_csv"],
        "mne_annotation_count": mne_annotation_count,
        "annotation_count_sources_agree": annotation_counts_agree(file_row, mne_annotation_count),
        "zero_annotation_fact": zero_annotation_fact,
        "low_annotation_fact": low_annotation_fact,
        "file_role_context": file_row["file_role"],
        "event_context_pending_codes": ";".join(event_codes),
        "dry_run_issue_codes": ";".join(issue_codes),
        "no_preprocessing_performed": True,
    }


def split_context_for_file(source_filename: str, context: dict[str, pd.DataFrame | None]) -> dict[str, Any]:
    """Summarize split-file audit rows for one source EDF.

    Args:
        source_filename: EDF filename.
        context: Optional event-context tables.

    Returns:
        Dictionary of split-file context fields for one event-context row.

    Side effects:
        None.
    """

    split = context.get("split_file_continuity_audit")
    if split is None or "source_filename" not in split.columns:
        return {
            "split_file_context_present": False,
            "split_file_context_row_kinds": "",
            "split_file_strict_clean_timing_rows": 0,
            "split_file_raw_only_rows": 0,
            "split_file_offset_only_rows": 0,
        }

    rows = split[split["source_filename"].map(text_or_empty).eq(source_filename)]
    return {
        "split_file_context_present": not rows.empty,
        "split_file_context_row_kinds": semicolon_join(rows["row_kind"]) if "row_kind" in rows else "",
        "split_file_strict_clean_timing_rows": int(pd.to_numeric(rows.get("strict_clean_timing_rows", 0), errors="coerce").fillna(0).sum()),
        "split_file_raw_only_rows": int(pd.to_numeric(rows.get("raw_only_rows", 0), errors="coerce").fillna(0).sum()),
        "split_file_offset_only_rows": int(pd.to_numeric(rows.get("offset_only_rows", 0), errors="coerce").fillna(0).sum()),
    }


def concatenated_context_for_file(source_filename: str, context: dict[str, pd.DataFrame | None]) -> dict[str, Any]:
    """Summarize ID 5 concatenated-file audit rows for one EDF.

    Args:
        source_filename: EDF filename.
        context: Optional event-context tables.

    Returns:
        Dictionary of concatenated-file context fields for one event-context
        row.

    Side effects:
        None.
    """

    id5 = context.get("id5_concatenated_file_audit")
    if id5 is None or "source_filename" not in id5.columns:
        return {
            "concatenated_context_present": False,
            "id5_file_start_onset_seconds": pd.NA,
            "id5_segments_with_strict_clean_timing": "",
        }

    rows = id5[id5["source_filename"].map(text_or_empty).eq(source_filename)]
    if rows.empty:
        return {
            "concatenated_context_present": False,
            "id5_file_start_onset_seconds": pd.NA,
            "id5_segments_with_strict_clean_timing": "",
        }

    return {
        "concatenated_context_present": True,
        "id5_file_start_onset_seconds": safe_float(rows["file_start_onset_seconds"].dropna().iloc[0])
        if "file_start_onset_seconds" in rows and not rows["file_start_onset_seconds"].dropna().empty
        else pd.NA,
        "id5_segments_with_strict_clean_timing": semicolon_join(
            rows.loc[rows.get("matched_rows_timing_clean_under_strict_definition", "").map(text_or_empty).eq("yes"), "segment"]
        )
        if {"matched_rows_timing_clean_under_strict_definition", "segment"}.issubset(rows.columns)
        else "",
    }


def zero_no_offset_context_for_file(source_filename: str, context: dict[str, pd.DataFrame | None]) -> dict[str, Any]:
    """Summarize zero/no-offset audit context for one EDF.

    Args:
        source_filename: EDF filename.
        context: Optional event-context tables.

    Returns:
        Dictionary of zero/no-offset context fields for one event-context row.

    Side effects:
        None.
    """

    zero = context.get("zero_no_offset_case_audit")
    if zero is None or "source_filenames" not in zero.columns:
        return {
            "zero_no_offset_context_present": False,
            "zero_no_offset_factual_note": "",
            "old_compatible_offset_rows": pd.NA,
        }

    mask = zero["source_filenames"].map(
        lambda value: source_filename in [part.strip() for part in text_or_empty(value).split(";") if part.strip()]
    )
    rows = zero[mask]
    if rows.empty:
        return {
            "zero_no_offset_context_present": False,
            "zero_no_offset_factual_note": "",
            "old_compatible_offset_rows": pd.NA,
        }

    return {
        "zero_no_offset_context_present": True,
        "zero_no_offset_factual_note": semicolon_join(rows["factual_note"]) if "factual_note" in rows else "",
        "old_compatible_offset_rows": int(
            pd.to_numeric(rows.get("old_compatible_offset_rows", 0), errors="coerce").fillna(0).sum()
        ),
    }


def proposed_join_context_for_file(source_filename: str, context: dict[str, pd.DataFrame | None]) -> str:
    """Summarize proposed join-status counts for one EDF.

    Args:
        source_filename: EDF filename.
        context: Optional event-context tables.

    Returns:
        Semicolon-delimited join-status count text.

    Side effects:
        None.
    """

    join = context.get("proposed_offset_join_audit")
    if join is None or not {"source_filename", "join_status"}.issubset(join.columns):
        return ""
    rows = join[join["source_filename"].map(text_or_empty).eq(source_filename)]
    return value_counts_text(rows["join_status"]) if not rows.empty else ""


def raw_only_context_for_file(source_filename: str, context: dict[str, pd.DataFrame | None]) -> dict[str, Any]:
    """Summarize raw-only row-classification context for one EDF.

    Args:
        source_filename: EDF filename.
        context: Optional event-context tables.

    Returns:
        Dictionary with raw-only row count and category counts.

    Side effects:
        None.
    """

    raw_only = context.get("raw_only_row_classification")
    if raw_only is None or "source_filename" not in raw_only.columns:
        return {"raw_only_context_rows": 0, "raw_only_category_counts": ""}

    rows = raw_only[raw_only["source_filename"].map(text_or_empty).eq(source_filename)]
    return {
        "raw_only_context_rows": int(len(rows)),
        "raw_only_category_counts": value_counts_text(rows["raw_only_category"]) if "raw_only_category" in rows else "",
    }


def build_event_context_flags(
    file_summary: pd.DataFrame,
    context: dict[str, pd.DataFrame | None],
) -> pd.DataFrame:
    """Build one row per EDF carrying event-audit context facts.

    Args:
        file_summary: Dry-run file summary table.
        context: Optional event-context tables.

    Returns:
        Event-context flag DataFrame.

    Side effects:
        None.
    """

    source_status = optional_context_source_status(context)
    rows: list[dict[str, Any]] = []
    for _, file_row in file_summary.iterrows():
        source_filename = text_or_empty(file_row["source_filename"])
        split_context = split_context_for_file(source_filename, context)
        concatenated_context = concatenated_context_for_file(source_filename, context)
        zero_context = zero_no_offset_context_for_file(source_filename, context)
        raw_only_context = raw_only_context_for_file(source_filename, context)

        rows.append(
            {
                "participant_id": file_row["participant_id"],
                "participant_id_padded": file_row["participant_id_padded"],
                "source_filename": source_filename,
                "file_role": file_row["file_role"],
                "split_part": file_row["split_part"],
                "annotation_count": file_row["mne_annotation_count"],
                "zero_annotation_fact": file_row["zero_annotation_fact"],
                "low_annotation_fact": file_row["low_annotation_fact"],
                **split_context,
                **concatenated_context,
                **zero_context,
                "proposed_join_status_counts": proposed_join_context_for_file(source_filename, context),
                **raw_only_context,
                "event_context_pending_codes": file_row["event_context_pending_codes"],
                "event_context_source_status": source_status,
            }
        )

    return pd.DataFrame(rows, columns=EVENT_CONTEXT_COLUMNS)


def run_dry_run_validations(
    inventory: pd.DataFrame,
    besa_coordinates: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Inspect EDF headers and build file, channel, and montage tables.

    Args:
        inventory: File inventory from manifests and raw-directory discovery.
        besa_coordinates: BESA coordinate table with normalized names.

    Returns:
        ``(file_summary, channel_type_validation, montage_validation)``.

    Side effects:
        Opens each existing EDF through MNE with ``preload=False``.
    """

    besa_by_normalized_name = {
        row["normalized_channel_name"]: row for _, row in besa_coordinates.drop_duplicates("normalized_channel_name").iterrows()
    }
    file_rows: list[dict[str, Any]] = []
    channel_rows: list[dict[str, Any]] = []
    montage_rows: list[dict[str, Any]] = []

    for _, inventory_row in inventory.iterrows():
        file_row = inventory_row.to_dict()
        if not file_row["file_exists"]:
            header = {
                "read_status": "raw_file_missing",
                "mne_read_error_type": "",
                "mne_read_error_message": "",
                "sampling_frequency_hz": None,
                "duration_seconds": None,
                "channel_names": [],
                "mne_channel_types": [],
                "mne_annotation_count": None,
                "mne_annotation_description_counts_json": "{}",
            }
        else:
            header = inspect_edf_header(Path(file_row["resolved_edf_path"]))

        file_montage_rows = build_montage_coordinate_rows(file_row, header, besa_by_normalized_name)
        file_channel_rows = build_channel_type_rows(file_row, header)
        file_summary_row = build_file_summary_row(file_row, header, file_montage_rows)

        file_rows.append(file_summary_row)
        channel_rows.extend(file_channel_rows)
        montage_rows.extend(file_montage_rows)

    return (
        pd.DataFrame(file_rows, columns=FILE_SUMMARY_COLUMNS),
        pd.DataFrame(channel_rows, columns=CHANNEL_TYPE_COLUMNS),
        pd.DataFrame(montage_rows, columns=MONTAGE_COLUMNS),
    )


def count_semicolon_codes(values: pd.Series) -> Counter[str]:
    """Count semicolon-delimited issue or context codes.

    Args:
        values: Series containing semicolon-delimited code text.

    Returns:
        Counter keyed by code.

    Side effects:
        None.
    """

    counts: Counter[str] = Counter()
    for value in values.fillna(""):
        for code in text_or_empty(value).split(";"):
            if code:
                counts[code] += 1
    return counts


def format_source_filenames(frame: pd.DataFrame) -> str:
    """Format source filenames from a DataFrame for Markdown summaries.

    Args:
        frame: DataFrame with a ``source_filename`` column.

    Returns:
        Backtick-quoted filenames joined by commas, or ``none`` when empty.

    Side effects:
        None.
    """

    if frame.empty or "source_filename" not in frame.columns:
        return "none"
    return ", ".join(f"`{name}`" for name in frame["source_filename"].dropna().map(text_or_empty).sort_values())


def coerce_split_part_for_output(frame: pd.DataFrame) -> pd.DataFrame:
    """Render split-part numbers as nullable integers in output CSVs.

    Args:
        frame: DataFrame that may have a ``split_part`` column.

    Returns:
        Copy of ``frame`` with ``split_part`` coerced to pandas nullable
        integer dtype when present.

    Side effects:
        None.
    """

    out = frame.copy()
    if "split_part" in out.columns:
        out["split_part"] = pd.to_numeric(out["split_part"], errors="coerce").astype("Int64")
    return out


def write_summary_markdown(
    output_path: Path,
    file_summary: pd.DataFrame,
    channel_validation: pd.DataFrame,
    montage_validation: pd.DataFrame,
    event_context: pd.DataFrame,
    optional_context: dict[str, pd.DataFrame | None],
) -> None:
    """Write the dry-run Markdown summary.

    Args:
        output_path: Destination Markdown file.
        file_summary: File-level dry-run summary table.
        channel_validation: Channel-type validation table.
        montage_validation: Montage coordinate validation table.
        event_context: Event-context flag table.
        optional_context: Optional event-context tables.

    Returns:
        ``None``.

    Side effects:
        Writes ``output_path``.
    """

    issue_counts = count_semicolon_codes(file_summary["dry_run_issue_codes"])
    event_code_counts = count_semicolon_codes(file_summary["event_context_pending_codes"])
    role_counts = file_summary["file_role"].value_counts(dropna=False).sort_index()
    read_status_counts = file_summary["read_status"].value_counts(dropna=False).sort_index()

    read_failures = file_summary[file_summary["read_status"].ne("ok")]
    channel_mismatch_files = file_summary[
        file_summary["missing_expected_channel_names"].fillna("").ne("")
        | file_summary["unexpected_channel_names"].fillna("").ne("")
        | ~file_summary["expected_37_channels_present"].fillna(False)
    ]
    montage_problem_files = file_summary[file_summary["besa_unmatched_eeg_channel_names"].fillna("").ne("")]
    zero_annotation_files = file_summary[file_summary["zero_annotation_fact"].fillna(False)]
    low_annotation_files = file_summary[file_summary["low_annotation_fact"].fillna(False)]

    raw_only_categories = Counter()
    raw_only = optional_context.get("raw_only_row_classification")
    if raw_only is not None and "raw_only_category" in raw_only.columns:
        raw_only_categories.update(raw_only["raw_only_category"].dropna().map(text_or_empty))

    lines = [
        "# Continuous Raw Preprocessing Dry-Run Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "This was a validation/QC-manifest dry run only. No filtering, ICA, bad-channel interpolation, CSD, event repair, epoching, cleaned FIF writing, or participant-retention decision was performed.",
        "",
        "## Files inspected",
        "",
        f"- File-summary rows: {len(file_summary)}",
        f"- Manifest rows present: {int(file_summary['manifest_present'].sum())}",
        f"- EDFs discovered in raw directory: {int(file_summary['discovered_in_raw_dir'].sum())}",
        "",
        "## File roles",
        "",
    ]
    for role, count in role_counts.items():
        lines.append(f"- {role}: {int(count)}")

    lines.extend(["", "## MNE read status", ""])
    for status, count in read_status_counts.items():
        lines.append(f"- {status}: {int(count)}")
    if not read_failures.empty:
        lines.append(f"- Files needing read/file follow-up: {format_source_filenames(read_failures)}")

    lines.extend(["", "## Dry-run issue counts", ""])
    if issue_counts:
        for code, count in sorted(issue_counts.items()):
            lines.append(f"- {code}: {count}")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Channel and coordinate checks",
            "",
            f"- Files with channel-count or label mismatches: {len(channel_mismatch_files)}",
            f"- Files with unmatched BESA EEG coordinates: {len(montage_problem_files)}",
            f"- Channel-type validation rows: {len(channel_validation)}",
            f"- Montage validation rows: {len(montage_validation)}",
        ]
    )
    if not channel_mismatch_files.empty:
        lines.append(f"- Channel-label follow-up files: {format_source_filenames(channel_mismatch_files)}")
    if not montage_problem_files.empty:
        lines.append(f"- BESA-coordinate follow-up files: {format_source_filenames(montage_problem_files)}")

    lines.extend(
        [
            "",
            "## Annotation-count facts",
            "",
            f"- Zero-annotation files: {len(zero_annotation_files)} ({format_source_filenames(zero_annotation_files)})",
            f"- Low nonzero-annotation files: {len(low_annotation_files)} ({format_source_filenames(low_annotation_files)})",
            "",
            "Annotation-count flags are carried as facts only.",
            "",
            "## Event-context facts carried forward",
            "",
        ]
    )
    if event_code_counts:
        for code, count in sorted(event_code_counts.items()):
            lines.append(f"- {code}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "Optional event-context table status:", ""])
    for label, frame in optional_context.items():
        status = "present" if frame is not None else "missing"
        row_count = len(frame) if frame is not None else 0
        lines.append(f"- {label}: {status}, rows={row_count}")

    lines.extend(["", "Raw-only row categories from optional context:", ""])
    if raw_only_categories:
        for category, count in sorted(raw_only_categories.items()):
            lines.append(f"- {category}: {count}")
    else:
        lines.append("- not available")

    split_context_files = event_context[event_context["split_file_context_present"].fillna(False)]
    concatenated_context_files = event_context[event_context["concatenated_context_present"].fillna(False)]
    zero_no_offset_context_files = event_context[event_context["zero_no_offset_context_present"].fillna(False)]

    lines.extend(
        [
            "",
            "Special event-policy context remains pending for split-file continuity, ID 5 file-start handling, zero/no-offset cases, raw-only row categories, and any future use of the proposed trigger-cleanup surface.",
            f"- Split context files: {format_source_filenames(split_context_files)}",
            f"- Concatenated/file-start context files: {format_source_filenames(concatenated_context_files)}",
            f"- Zero/no-offset context files: {format_source_filenames(zero_no_offset_context_files)}",
            "",
            "## Output files",
            "",
            f"- `{FILE_SUMMARY_FILENAME}`",
            f"- `{CHANNEL_TYPE_VALIDATION_FILENAME}`",
            f"- `{MONTAGE_COORDINATE_VALIDATION_FILENAME}`",
            f"- `{EVENT_CONTEXT_FLAGS_FILENAME}`",
            "",
            "These outputs are local-only under `_Data/`.",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(
    output_dir: Path,
    file_summary: pd.DataFrame,
    channel_validation: pd.DataFrame,
    montage_validation: pd.DataFrame,
    event_context: pd.DataFrame,
    optional_context: dict[str, pd.DataFrame | None],
) -> None:
    """Write all dry-run output tables and the Markdown summary.

    Args:
        output_dir: Dry-run output directory.
        file_summary: File-level dry-run summary table.
        channel_validation: Channel-type validation table.
        montage_validation: Montage coordinate validation table.
        event_context: Event-context flag table.
        optional_context: Optional event-context tables.

    Returns:
        ``None``.

    Side effects:
        Creates ``output_dir`` when needed and writes local-only outputs below
        it.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    file_summary = coerce_split_part_for_output(file_summary)
    channel_validation = coerce_split_part_for_output(channel_validation)
    montage_validation = coerce_split_part_for_output(montage_validation)
    event_context = coerce_split_part_for_output(event_context)

    file_summary.to_csv(output_dir / FILE_SUMMARY_FILENAME, index=False)
    channel_validation.to_csv(output_dir / CHANNEL_TYPE_VALIDATION_FILENAME, index=False)
    montage_validation.to_csv(output_dir / MONTAGE_COORDINATE_VALIDATION_FILENAME, index=False)
    event_context.to_csv(output_dir / EVENT_CONTEXT_FLAGS_FILENAME, index=False)
    write_summary_markdown(
        output_dir / DRY_RUN_SUMMARY_FILENAME,
        file_summary,
        channel_validation,
        montage_validation,
        event_context,
        optional_context,
    )


def main() -> None:
    """Run the continuous raw preprocessing dry run.

    Args:
        None. All paths are resolved relative to the repository root.

    Returns:
        ``None``.

    Side effects:
        Reads required local manifests, the BESA coordinate CSV, optional
        event-audit context tables, and EDF headers with MNE ``preload=False``.
        Writes local-only QC outputs under
        ``_Data/eeg/mne_preprocessing/dry_run/``.
    """

    repo_root = repo_root_from_script()
    raw_dir = repo_root / RAW_EEG_DIR
    output_dir = repo_root / DRY_RUN_OUTPUT_DIR

    manifest = read_csv_checked(
        repo_root / RAW_MANIFEST_PATH,
        {"participant_id", "source_filename", "file_path", "file_role", "annotation_count"},
        "raw EEG file manifest",
    )
    annotation_counts = read_csv_checked(
        repo_root / ANNOTATION_COUNTS_PATH,
        {"source_filename", "annotation_description", "annotation_count"},
        "raw EEG annotation-count CSV",
    )
    besa_coordinates = load_besa_coordinates(repo_root / BESA_COORDINATES_PATH)
    optional_context = read_optional_event_context(repo_root)

    edf_files = discover_edf_files(raw_dir)
    annotation_totals = aggregate_annotation_counts(annotation_counts)
    inventory = build_file_inventory(manifest, annotation_totals, edf_files, raw_dir, repo_root)
    file_summary, channel_validation, montage_validation = run_dry_run_validations(inventory, besa_coordinates)
    event_context = build_event_context_flags(file_summary, optional_context)

    write_outputs(output_dir, file_summary, channel_validation, montage_validation, event_context, optional_context)

    print(f"Dry-run inspected {len(file_summary)} EDF file row(s).")
    read_status_counts = {
        str(status): int(count)
        for status, count in file_summary["read_status"].value_counts(dropna=False).sort_index().items()
    }
    print(f"MNE read status counts: {read_status_counts}")
    print(f"Wrote dry-run QC outputs to {output_dir}")


if __name__ == "__main__":
    main()
