"""Config-driven raw QC for future continuous MNE preprocessing.

This script is the first continuous-raw preprocessing entry point for the new
DEMI EEG reanalysis path, but it deliberately does not preprocess signals. It
now reads a private raw-QC configuration, checks that every signal-mutating or
event-dependent stage is still disabled, opens each raw EDF with MNE using
``preload=False`` for header-level validation, applies approved channel-type
overrides in memory, validates EEG channel names against the tracked BESA-81
coordinate table after case normalization, and writes local raw-QC summaries.

The first raw-QC stage comes before actual continuous preprocessing because the
active MNE workflow is restarting from raw EDF recordings, while final
preprocessing, event, and epoch policies are still pending. Header, channel,
coordinate, file-role, annotation-count, sampled raw-amplitude, sampled raw-PSD,
and event-context checks can be made now without filtering data, notch
filtering, re-referencing, repairing triggers, constructing epochs, or making
participant-level analytic decisions. This script therefore creates a review
surface for Tony before any later script writes cleaned continuous FIF files.

The script reads these required local inputs:

- ``_Data/eeg/manifest/raw_eeg_file_manifest.csv``;
- ``_Data/eeg/manifest/raw_eeg_annotation_counts.csv``;
- ``_Data/eeg/BESA-81.csv``;
- raw EDF files discovered directly below ``_Data/eeg/raw/``;
- by default, ``_Private/config/preprocessing_raw_qc_config.yaml``.

When present, it also reads these local event-audit context tables:

- ``_Data/eeg/event_sequence_audit/proposed_offset_join_audit.csv``;
- ``_Data/eeg/event_special_case_audit/split_file_continuity_audit.csv``;
- ``_Data/eeg/event_special_case_audit/id5_concatenated_file_audit.csv``;
- ``_Data/eeg/event_special_case_audit/zero_no_offset_case_audit.csv``;
- ``_Data/eeg/event_special_case_audit/raw_only_row_classification.csv``.

Generated outputs are local-only raw-QC files under
``_Data/eeg/mne_preprocessing/raw_qc/`` by default:

- ``raw_qc_run_manifest.json``;
- ``preprocessing_raw_qc_file_summary.csv``;
- ``channel_type_validation.csv``;
- ``montage_coordinate_validation.csv``;
- ``event_context_flags.csv``;
- ``amplitude_range_summary_by_channel_type.csv``;
- ``raw_psd_summary.csv``;
- ``raw_time_series_snapshot_figures.csv``;
- ``preprocessing_raw_qc_summary.md``;
- raw PSD and time-series snapshot figures below ``figures/`` when sampled
  signal reads succeed.

What this script explicitly does not do:

- no EEG filtering;
- no notch filtering;
- no re-referencing;
- no bad-channel detection;
- no ICA;
- no bad-channel interpolation;
- no CSD computation;
- no EMG movement-summary computation;
- no EMG-derived trial handling;
- no epoch construction;
- no cleaned FIF or EDF derivative writing;
- no event repair;
- no participant-level decision.

Safety boundaries:

- EDFs are opened with ``mne.io.read_raw_edf(..., preload=False)`` for
  header-level validation.
- Signal reads for raw QC are sampled, one EDF at a time, and are not saved as
  modified MNE objects.
- Raw EDF files are not modified.
- One missing or problematic EDF records a factual row and does not stop the
  raw-QC run.
- Optional event-audit tables are treated as context only. Their facts are not
  converted into event policy.
- All generated files are below ignored ``_Data/`` paths.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/04_preprocess_continuous_raw.py

To use a different private config path:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/04_preprocess_continuous_raw.py --config _Private/config/preprocessing_raw_qc_config.yaml

Paths are resolved relative to this script, so the command also works when
launched from another current working directory.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import yaml

# MNE may import matplotlib depending on the local environment. Keep any cache
# writes inside a temp directory so this raw-QC script does not depend on
# the user's global matplotlib configuration path.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "demi_matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import pandas as pd
from scipy.signal import welch


# Keep the small unit contract importable both when this numbered script is
# executed directly and when tests load it by file path.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from channel_qc import classify_unit_status  # noqa: E402


RAW_EEG_DIR = Path("_Data") / "eeg" / "raw"
RAW_MANIFEST_PATH = Path("_Data") / "eeg" / "manifest" / "raw_eeg_file_manifest.csv"
ANNOTATION_COUNTS_PATH = Path("_Data") / "eeg" / "manifest" / "raw_eeg_annotation_counts.csv"
BESA_COORDINATES_PATH = Path("_Data") / "eeg" / "BESA-81.csv"
DEFAULT_CONFIG_PATH = Path("_Private") / "config" / "preprocessing_raw_qc_config.yaml"
RAW_QC_OUTPUT_DIR = Path("_Data") / "eeg" / "mne_preprocessing" / "raw_qc"

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

RUN_MANIFEST_FILENAME = "raw_qc_run_manifest.json"
FILE_SUMMARY_FILENAME = "preprocessing_raw_qc_file_summary.csv"
CHANNEL_TYPE_VALIDATION_FILENAME = "channel_type_validation.csv"
MONTAGE_COORDINATE_VALIDATION_FILENAME = "montage_coordinate_validation.csv"
EVENT_CONTEXT_FLAGS_FILENAME = "event_context_flags.csv"
AMPLITUDE_RANGE_FILENAME = "amplitude_range_summary_by_channel_type.csv"
PSD_SUMMARY_FILENAME = "raw_psd_summary.csv"
FIGURE_MANIFEST_FILENAME = "raw_time_series_snapshot_figures.csv"
RAW_QC_SUMMARY_FILENAME = "preprocessing_raw_qc_summary.md"
FIGURE_DIRNAME = "figures"
TIME_SERIES_FIGURE_DIRNAME = "time_series_snapshots"
PSD_FIGURE_DIRNAME = "raw_psd"

EXPECTED_SAMPLING_FREQUENCY_HZ = 1000.0
SAMPLING_FREQUENCY_TOLERANCE_HZ = 0.01
EXPECTED_CHANNEL_COUNT = 37
LOW_ANNOTATION_COUNT_THRESHOLD = 600
RAW_QC_SAMPLE_WINDOW_SECONDS = 10.0
RAW_QC_SAMPLE_FRACTIONS = (0.10, 0.50, 0.90)
RAW_QC_PSD_MAX_FREQUENCY_HZ = 120.0
RAW_QC_PSD_NPERSEG_SECONDS = 2.0
RAW_QC_FIGURE_DPI = 120

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
    "channel_type_override_status",
    "channel_type_override_error_type",
    "channel_type_override_error_message",
    "eeg_channels_matched_to_besa",
    "besa_unmatched_eeg_channel_names",
    "standard_1005_expected_eeg_names_matched",
    "standard_1005_unmatched_eeg_channel_names",
    "manifest_annotation_count",
    "annotation_count_from_counts_csv",
    "mne_annotation_count",
    "annotation_count_sources_agree",
    "zero_annotation_fact",
    "low_annotation_fact",
    "file_role_context",
    "event_context_pending_codes",
    "raw_qc_issue_codes",
    "raw_qc_only_no_signal_mutation",
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
    "mne_reported_type_before_override",
    "channel_type_after_overrides",
    "expected_raw_qc_type",
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
    "standard_1005_name_match",
    "standard_1005_channel_name",
    "standard_1005_validation_only",
    "standard_1005_issue_codes",
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

AMPLITUDE_RANGE_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_role",
    "split_part",
    "raw_qc_status",
    "raw_qc_error_type",
    "raw_qc_error_message",
    "sample_strategy",
    "sample_window_count",
    "sampled_seconds_per_channel",
    "channel_type",
    "channel_group",
    "n_channels",
    "n_samples_per_channel",
    "mne_data_units",
    "amplitude_min",
    "amplitude_max",
    "amplitude_peak_to_peak",
    "amplitude_mean",
    "amplitude_median",
    "amplitude_std",
    "amplitude_rms",
    "amplitude_mean_abs",
]

PSD_SUMMARY_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_role",
    "split_part",
    "raw_qc_status",
    "raw_qc_error_type",
    "raw_qc_error_message",
    "sample_strategy",
    "sample_window_count",
    "sampled_seconds_per_channel",
    "channel_type",
    "channel_group",
    "n_channels",
    "welch_nperseg",
    "frequency_min_hz",
    "frequency_max_hz",
    "median_psd",
    "mean_psd",
    "median_psd_1_4_hz",
    "median_psd_4_8_hz",
    "median_psd_8_13_hz",
    "median_psd_13_30_hz",
    "median_psd_30_50_hz",
    "median_psd_55_65_hz",
    "psd_note",
]

FIGURE_MANIFEST_COLUMNS = [
    "participant_id",
    "participant_id_padded",
    "source_filename",
    "file_role",
    "split_part",
    "figure_type",
    "figure_path",
    "raw_qc_status",
    "raw_qc_error_type",
    "raw_qc_error_message",
    "figure_note",
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


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the raw-QC run.

    Args:
        None.

    Returns:
        Parsed arguments. The only current option is ``--config``, which
        defaults to the private raw-QC config approved for this stage.

    Side effects:
        Reads command-line arguments from ``sys.argv`` through ``argparse``.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Run DEMI raw EDF QC using the approved private config. "
            "The script fails closed if any signal-mutating stage is enabled."
        )
    )
    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG_PATH.as_posix(),
        help=(
            "Path to the private raw-QC YAML config. Relative paths are resolved "
            "from the repository root."
        ),
    )
    return parser.parse_args()


def resolve_repo_path(repo_root: Path, value: Any) -> Path:
    """Resolve a config path relative to the repository root.

    Args:
        repo_root: Absolute repository root.
        value: Path-like value from the config or command line.

    Returns:
        Absolute ``Path``.

    Side effects:
        None.
    """

    path = Path(text_or_empty(value))
    return path if path.is_absolute() else repo_root / path


def load_raw_qc_config(config_path: Path) -> dict[str, Any]:
    """Read the private YAML raw-QC config.

    Args:
        config_path: Absolute path to a YAML config.

    Returns:
        Parsed config dictionary.

    Side effects:
        Reads ``config_path``. Raises ``RuntimeError`` for missing, empty, or
        malformed configs. No config values are applied here.
    """

    if not config_path.exists():
        raise RuntimeError(f"raw-QC config is missing: {config_path.as_posix()}")

    try:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        raise RuntimeError(f"failed to parse raw-QC YAML config: {error}") from error

    if not isinstance(config, dict):
        raise RuntimeError("raw-QC config must parse to a mapping")
    return config


def config_value(config: dict[str, Any], path: tuple[str, ...], default: Any = None) -> Any:
    """Return a nested config value with a default.

    Args:
        config: Parsed raw-QC config.
        path: Sequence of nested mapping keys.
        default: Value returned when any key is absent.

    Returns:
        Config value or ``default``.

    Side effects:
        None.
    """

    current: Any = config
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def require_config_value(config: dict[str, Any], path: tuple[str, ...]) -> Any:
    """Return a nested config value, failing if it is absent.

    Args:
        config: Parsed raw-QC config.
        path: Sequence of nested mapping keys.

    Returns:
        Config value at ``path``.

    Side effects:
        Raises ``RuntimeError`` when the path is missing.
    """

    marker = object()
    value = config_value(config, path, marker)
    if value is marker:
        raise RuntimeError(f"raw-QC config is missing required key: {'.'.join(path)}")
    return value


def configured_channel_type_by_normalized_name(config: dict[str, Any]) -> dict[str, str]:
    """Build the approved channel-type map from the config.

    Args:
        config: Parsed raw-QC config.

    Returns:
        Mapping from normalized channel names to MNE channel types.

    Side effects:
        Raises ``RuntimeError`` if the configured mapping is not a mapping of
        channel-type labels to channel-name lists.
    """

    type_mapping = require_config_value(config, ("channels", "type_mapping"))
    if not isinstance(type_mapping, dict):
        raise RuntimeError("channels.type_mapping must be a mapping")

    out: dict[str, str] = {}
    for channel_type, channel_names in type_mapping.items():
        if not isinstance(channel_names, list):
            raise RuntimeError(f"channels.type_mapping.{channel_type} must be a list")
        for channel_name in channel_names:
            normalized_name = normalize_channel_name(channel_name)
            if normalized_name in out and out[normalized_name] != channel_type:
                raise RuntimeError(f"channel {normalized_name} has conflicting configured types")
            out[normalized_name] = text_or_empty(channel_type).lower()
    return out


def configured_overrides_by_normalized_name(config: dict[str, Any]) -> dict[str, str]:
    """Return the in-memory channel-type overrides approved by the config.

    Args:
        config: Parsed raw-QC config.

    Returns:
        Mapping from normalized raw channel names to MNE channel types.

    Side effects:
        Raises ``RuntimeError`` if the override section is malformed.
    """

    overrides = require_config_value(config, ("channels", "future_mne_type_overrides"))
    if not isinstance(overrides, dict):
        raise RuntimeError("channels.future_mne_type_overrides must be a mapping")
    return {normalize_channel_name(channel_name): text_or_empty(channel_type).lower() for channel_name, channel_type in overrides.items()}


def validate_raw_qc_config(config: dict[str, Any], config_path: Path, repo_root: Path) -> None:
    """Fail closed unless the config matches the approved raw-QC boundary.

    Args:
        config: Parsed raw-QC config.
        config_path: Path the config was read from, used in error messages.
        repo_root: Absolute repository root.

    Returns:
        ``None``.

    Side effects:
        Raises ``RuntimeError`` before any EDF is read if a disabled stage is
        enabled, if cleaned outputs are requested, or if the channel mapping no
        longer matches the approved raw-QC scope.
    """

    disabled_sections = (
        "filtering",
        "notch_filtering",
        "reference",
        "bad_channel_detection",
        "interpolation",
        "ica",
        "csd",
        "cleaned_fif_writing",
        "event_repair",
        "epoching",
    )
    for section in disabled_sections:
        enabled = require_config_value(config, ("disabled_signal_mutating_sections", section, "enabled"))
        if bool(enabled):
            raise RuntimeError(f"raw-QC config enables disabled stage: {section}")

    fail_closed_false_paths = (
        ("outputs", "write_cleaned_signal_derivatives"),
        ("outputs", "write_cleaned_fif"),
        ("outputs", "write_preprocessed_edf"),
        ("outputs", "write_epochs"),
        ("raw_qc_outputs", "signal_mutation_allowed"),
        ("raw_qc_outputs", "event_context", "make_event_policy_decisions"),
        ("channels", "trigger_policy", "use_for_event_repair"),
        ("channels", "trigger_policy", "use_for_epoching"),
        ("montage_validation", "apply_active_montage_to_signal"),
    )
    for path in fail_closed_false_paths:
        if bool(require_config_value(config, path)):
            raise RuntimeError(f"raw-QC config value must stay false: {'.'.join(path)}")

    if not bool(require_config_value(config, ("raw_qc_outputs", "enabled"))):
        raise RuntimeError("raw-QC config must enable raw_qc_outputs.enabled")
    if not bool(require_config_value(config, ("outputs", "local_only"))):
        raise RuntimeError("raw-QC config must keep outputs.local_only=true")

    output_root = resolve_repo_path(repo_root, require_config_value(config, ("outputs", "output_root")))
    try:
        output_relative = output_root.resolve().relative_to(repo_root)
    except ValueError as error:
        raise RuntimeError("raw-QC output_root must stay inside the repository") from error
    if not output_relative.parts or output_relative.parts[0] != "_Data":
        raise RuntimeError("raw-QC output_root must stay under _Data/")

    expected_mapping = EXPECTED_CHANNEL_TYPE_BY_NAME
    configured_mapping = configured_channel_type_by_normalized_name(config)
    if configured_mapping != expected_mapping:
        raise RuntimeError(
            "raw-QC config channel mapping differs from the approved DEMI channel set; "
            f"config={config_path.as_posix()}"
        )

    expected_overrides = {
        "HEO": "eog",
        "VEO": "eog",
        "EMG-L": "emg",
        "EMG-A": "emg",
        "TRIGGER": "stim",
    }
    configured_overrides = configured_overrides_by_normalized_name(config)
    for channel_name, expected_type in expected_overrides.items():
        if configured_overrides.get(channel_name) != expected_type:
            raise RuntimeError(
                "raw-QC config must keep the approved in-memory override "
                f"{channel_name}={expected_type}"
            )


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


def optional_event_context_paths(config: dict[str, Any]) -> dict[str, Path]:
    """Return optional event-context paths from the config.

    Args:
        config: Parsed raw-QC config.

    Returns:
        Mapping from context-table label to repository-relative path.

    Side effects:
        None.
    """

    configured = config_value(config, ("inputs", "optional_event_context_tables"), {})
    if not isinstance(configured, dict):
        raise RuntimeError("inputs.optional_event_context_tables must be a mapping when present")
    return {label: Path(text_or_empty(path)) for label, path in configured.items()}


def read_optional_event_context(repo_root: Path, context_paths: dict[str, Path]) -> dict[str, pd.DataFrame | None]:
    """Read optional event-audit context tables when they are present.

    Args:
        repo_root: Absolute repository root path.
        context_paths: Mapping from context-table labels to paths from the
            raw-QC config.

    Returns:
        Dictionary keyed by context-table label. Values are DataFrames for
        present tables and ``None`` for absent tables.

    Side effects:
        Reads local CSV files below ``_Data/`` when they exist.
    """

    context: dict[str, pd.DataFrame | None] = {}
    for label, relative_path in context_paths.items():
        path = resolve_repo_path(repo_root, relative_path)
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


def standard_1005_name_map_if_requested(config: dict[str, Any]) -> tuple[dict[str, str] | None, str]:
    """Build a validation-only MNE ``standard_1005`` name map when requested.

    Args:
        config: Parsed raw-QC config.

    Returns:
        ``(mapping, status)``. ``mapping`` is ``None`` when the check is not
        requested or cannot be built. ``status`` records what happened for the
        run manifest and Markdown summary.

    Side effects:
        Calls ``mne.channels.make_standard_montage`` when the config requests
        the validation-only check. It does not apply the montage to any Raw
        object.
    """

    standard_config = config_value(config, ("montage_validation", "standard_1005_check"), {})
    if not isinstance(standard_config, dict) or not bool(standard_config.get("enabled_if_straightforward", False)):
        return None, "not_requested"
    if bool(standard_config.get("active_montage", False)):
        raise RuntimeError("standard_1005_check.active_montage must remain false")
    if not bool(standard_config.get("validation_only", False)):
        raise RuntimeError("standard_1005_check.validation_only must remain true")

    try:
        montage = mne.channels.make_standard_montage("standard_1005")
    except Exception as error:  # noqa: BLE001 - validation-only check should not stop raw QC.
        return None, f"not_available:{type(error).__name__}:{error}"

    return {normalize_channel_name(channel_name): channel_name for channel_name in montage.ch_names}, "available"


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
    """Return the raw-QC channel type expected for a normalized channel name.

    Args:
        normalized_channel_name: Uppercase channel label.

    Returns:
        ``eeg``, ``eog``, ``emg``, ``stim``, or ``unmapped``.

    Side effects:
        None.
    """

    return EXPECTED_CHANNEL_TYPE_BY_NAME.get(normalized_channel_name, "unmapped")


def channel_group_for_type(channel_type: str) -> str:
    """Return a broad channel group label for raw-QC output.

    Args:
        channel_type: Expected raw-QC channel type.

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


def apply_configured_channel_type_overrides(
    raw: mne.io.BaseRaw,
    overrides_by_normalized_name: dict[str, str],
) -> dict[str, str]:
    """Apply approved channel-type overrides to an in-memory MNE Raw object.

    Args:
        raw: Raw EDF object opened by MNE. The object may use
            ``preload=False``.
        overrides_by_normalized_name: Mapping from normalized channel names to
            MNE channel types.

    Returns:
        Mapping from raw channel labels present in ``raw`` to the channel type
        applied.

    Side effects:
        Calls ``raw.set_channel_types`` on the in-memory Raw object. This does
        not edit the source EDF and does not save a derivative.
    """

    present_overrides = {
        channel_name: overrides_by_normalized_name[normalize_channel_name(channel_name)]
        for channel_name in raw.ch_names
        if normalize_channel_name(channel_name) in overrides_by_normalized_name
    }
    if not present_overrides:
        return {}

    try:
        raw.set_channel_types(present_overrides, on_unit_change="ignore")
    except TypeError:
        # Older MNE versions did not expose ``on_unit_change``. The fallback is
        # still an in-memory metadata change only.
        raw.set_channel_types(present_overrides)
    return present_overrides


def inspect_edf_header(edf_path: Path, overrides_by_normalized_name: dict[str, str]) -> dict[str, Any]:
    """Read EDF header metadata through MNE without loading signal data.

    Args:
        edf_path: Path to one raw EDF file.
        overrides_by_normalized_name: Approved in-memory channel-type override
            mapping from the raw-QC config.

    Returns:
        Dictionary with read status, sampling frequency, duration, channel
        labels, original MNE-reported channel types, post-override channel
        types, and annotation count. Error rows keep the exception type and
        message.

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
        mne_channel_types_before = list(raw.get_channel_types())
        override_status = "not_needed"
        override_error_type = ""
        override_error_message = ""
        applied_overrides: dict[str, str] = {}
        try:
            applied_overrides = apply_configured_channel_type_overrides(raw, overrides_by_normalized_name)
            override_status = "applied" if applied_overrides else "not_needed"
        except Exception as error:  # noqa: BLE001 - record the problem and continue the file audit.
            override_status = "failed"
            override_error_type = type(error).__name__
            override_error_message = str(error)
        channel_types_after = list(raw.get_channel_types())
        annotation_count = int(len(raw.annotations))
        annotation_description_counts = Counter(str(description) for description in raw.annotations.description)

        return {
            "read_status": "ok",
            "mne_read_error_type": "",
            "mne_read_error_message": "",
            "sampling_frequency_hz": sfreq,
            "duration_seconds": duration_seconds,
            "channel_names": channel_names,
            "mne_channel_types_before_overrides": mne_channel_types_before,
            "channel_types_after_overrides": channel_types_after,
            "applied_channel_type_overrides_json": json.dumps(applied_overrides, sort_keys=True),
            "channel_type_override_status": override_status,
            "channel_type_override_error_type": override_error_type,
            "channel_type_override_error_message": override_error_message,
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
            "mne_channel_types_before_overrides": [],
            "channel_types_after_overrides": [],
            "applied_channel_type_overrides_json": "{}",
            "channel_type_override_status": "not_attempted",
            "channel_type_override_error_type": "",
            "channel_type_override_error_message": "",
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
                "mne_reported_type_before_override": "",
                "channel_type_after_overrides": "",
                "expected_raw_qc_type": "",
                "channel_type_assignment_status": "file_not_read",
                "channel_group": "",
                "channel_type_issue_codes": header["read_status"],
            }
        ]

    rows: list[dict[str, Any]] = []
    for index, (channel_name, mne_type_before, channel_type_after) in enumerate(
        zip(
            header["channel_names"],
            header["mne_channel_types_before_overrides"],
            header["channel_types_after_overrides"],
        ),
        start=1,
    ):
        normalized_name = normalize_channel_name(channel_name)
        raw_qc_type = expected_channel_type(normalized_name)

        # Channel typing is intentionally name-based here. EDF channel metadata
        # can be sparse, while the DEMI acquisition labels have stable names for
        # scalp EEG, eye, muscle, and trigger channels. The approved raw-QC
        # config applies those labels in memory before any signal summaries.
        if raw_qc_type == "unmapped":
            status = "unmapped_channel_label"
            issue_codes = "unmapped_channel_label"
        elif raw_qc_type == channel_type_after:
            status = "matches_configured_type_after_override"
            issue_codes = ""
        else:
            status = "configured_type_mismatch_after_override"
            issue_codes = "channel_type_after_override_differs_from_expected"

        rows.append(
            {
                **file_identity_fields(file_row),
                "read_status": header["read_status"],
                "channel_index": index,
                "raw_channel_name": channel_name,
                "normalized_channel_name": normalized_name,
                "mne_reported_type_before_override": mne_type_before,
                "channel_type_after_overrides": channel_type_after,
                "expected_raw_qc_type": raw_qc_type,
                "channel_type_assignment_status": status,
                "channel_group": channel_group_for_type(raw_qc_type),
                "channel_type_issue_codes": issue_codes,
            }
        )

    return rows


def build_montage_coordinate_rows(
    file_row: dict[str, Any],
    header: dict[str, Any],
    besa_by_normalized_name: dict[str, pd.Series],
    standard_1005_by_normalized_name: dict[str, str] | None,
) -> list[dict[str, Any]]:
    """Build EEG channel-to-BESA coordinate validation rows.

    Args:
        file_row: File inventory row as a dictionary.
        header: Header metadata returned by ``inspect_edf_header``.
        besa_by_normalized_name: Mapping from normalized BESA labels to
            coordinate rows.
        standard_1005_by_normalized_name: Optional validation-only mapping from
            normalized MNE ``standard_1005`` labels to original montage labels.

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
                "standard_1005_name_match": False,
                "standard_1005_channel_name": "",
                "standard_1005_validation_only": standard_1005_by_normalized_name is not None,
                "standard_1005_issue_codes": header["read_status"],
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
        standard_1005_name = (
            standard_1005_by_normalized_name.get(normalized_name, "")
            if standard_1005_by_normalized_name is not None
            else ""
        )
        standard_1005_match = bool(standard_1005_name)

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

        standard_issues: list[str] = []
        if standard_1005_by_normalized_name is not None and expected_eeg_channel and not standard_1005_match:
            standard_issues.append("expected_eeg_channel_without_standard_1005_name_match")

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
                "standard_1005_name_match": standard_1005_match,
                "standard_1005_channel_name": standard_1005_name,
                "standard_1005_validation_only": standard_1005_by_normalized_name is not None,
                "standard_1005_issue_codes": ";".join(standard_issues),
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
    """Build one file-level raw-QC summary row.

    Args:
        file_row: File inventory row as a dictionary.
        header: Header metadata returned by ``inspect_edf_header``.
        montage_rows: Montage validation rows for the same file.

    Returns:
        Dictionary for ``preprocessing_raw_qc_file_summary.csv``.

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
    channel_type_mismatch_names = [
        row["raw_channel_name"]
        for row in channel_type_rows
        if row.get("channel_type_assignment_status") == "configured_type_mismatch_after_override"
    ]
    montage_problem_rows = [
        row for row in montage_rows if text_or_empty(row.get("montage_issue_codes")) and row.get("read_status") == "ok"
    ]
    unmatched_besa_names = [
        row["normalized_channel_name"]
        for row in montage_problem_rows
        if "without_besa_coordinate" in text_or_empty(row.get("montage_issue_codes"))
    ]
    standard_1005_unmatched_names = [
        row["normalized_channel_name"]
        for row in montage_rows
        if "without_standard_1005_name_match" in text_or_empty(row.get("standard_1005_issue_codes"))
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
    if channel_type_mismatch_names:
        issue_codes.append("channel_type_after_override_mismatch")
    if header.get("channel_type_override_status") == "failed":
        issue_codes.append("channel_type_override_failed")
    if unmatched_besa_names:
        issue_codes.append("besa_coordinate_unmatched_eeg_label")
    if standard_1005_unmatched_names:
        issue_codes.append("standard_1005_name_unmatched_eeg_label")
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
        "channel_type_assignment_complete": not unmapped_channel_names
        and not channel_type_mismatch_names
        and header["read_status"] == "ok"
        and header.get("channel_type_override_status") != "failed",
        "unmapped_channel_type_names": semicolon_join(unmapped_channel_names),
        "channel_type_override_status": header.get("channel_type_override_status", ""),
        "channel_type_override_error_type": header.get("channel_type_override_error_type", ""),
        "channel_type_override_error_message": header.get("channel_type_override_error_message", ""),
        "eeg_channels_matched_to_besa": not unmatched_besa_names and header["read_status"] == "ok",
        "besa_unmatched_eeg_channel_names": semicolon_join(unmatched_besa_names),
        "standard_1005_expected_eeg_names_matched": not standard_1005_unmatched_names and header["read_status"] == "ok",
        "standard_1005_unmatched_eeg_channel_names": semicolon_join(standard_1005_unmatched_names),
        "manifest_annotation_count": file_row["manifest_annotation_count"],
        "annotation_count_from_counts_csv": file_row["annotation_count_from_counts_csv"],
        "mne_annotation_count": mne_annotation_count,
        "annotation_count_sources_agree": annotation_counts_agree(file_row, mne_annotation_count),
        "zero_annotation_fact": zero_annotation_fact,
        "low_annotation_fact": low_annotation_fact,
        "file_role_context": file_row["file_role"],
        "event_context_pending_codes": ";".join(event_codes),
        "raw_qc_issue_codes": ";".join(issue_codes),
        "raw_qc_only_no_signal_mutation": True,
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
        file_summary: Raw-QC file summary table.
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


def run_header_validations(
    inventory: pd.DataFrame,
    besa_coordinates: pd.DataFrame,
    overrides_by_normalized_name: dict[str, str],
    standard_1005_by_normalized_name: dict[str, str] | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Inspect EDF headers and build file, channel, and montage tables.

    Args:
        inventory: File inventory from manifests and raw-directory discovery.
        besa_coordinates: BESA coordinate table with normalized names.
        overrides_by_normalized_name: Approved in-memory channel-type
            overrides from the config.
        standard_1005_by_normalized_name: Optional validation-only
            ``standard_1005`` channel-name mapping.

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
                "mne_channel_types_before_overrides": [],
                "channel_types_after_overrides": [],
                "applied_channel_type_overrides_json": "{}",
                "channel_type_override_status": "not_attempted",
                "channel_type_override_error_type": "",
                "channel_type_override_error_message": "",
                "mne_annotation_count": None,
                "mne_annotation_description_counts_json": "{}",
            }
        else:
            header = inspect_edf_header(Path(file_row["resolved_edf_path"]), overrides_by_normalized_name)

        file_montage_rows = build_montage_coordinate_rows(
            file_row,
            header,
            besa_by_normalized_name,
            standard_1005_by_normalized_name,
        )
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


def sample_window_bounds(n_times: int, sfreq: float) -> list[tuple[int, int]]:
    """Choose short raw-data windows for sampled signal QC.

    Args:
        n_times: Number of samples in the EDF.
        sfreq: Sampling frequency in Hz.

    Returns:
        List of ``(start, stop)`` sample bounds. Bounds are unique, sorted, and
        limited to the available recording length.

    Side effects:
        None.
    """

    if n_times <= 0 or sfreq <= 0:
        return []
    window_samples = max(1, min(int(round(RAW_QC_SAMPLE_WINDOW_SECONDS * sfreq)), n_times))
    bounds: dict[tuple[int, int], None] = {}
    for fraction in RAW_QC_SAMPLE_FRACTIONS:
        center = int(round(n_times * fraction))
        start = max(0, min(center - window_samples // 2, n_times - window_samples))
        stop = min(n_times, start + window_samples)
        bounds[(start, stop)] = None
    return sorted(bounds.keys())


def sampled_seconds(bounds: list[tuple[int, int]], sfreq: float) -> float:
    """Return sampled seconds represented by window bounds.

    Args:
        bounds: Sample bounds from ``sample_window_bounds``.
        sfreq: Sampling frequency in Hz.

    Returns:
        Total sampled duration in seconds.

    Side effects:
        None.
    """

    if sfreq <= 0:
        return 0.0
    return float(sum(stop - start for start, stop in bounds) / sfreq)


def sample_strategy_text(bounds: list[tuple[int, int]]) -> str:
    """Describe the raw-QC signal sampling strategy.

    Args:
        bounds: Sample bounds chosen for one file.

    Returns:
        Short text stored in raw-QC tables.

    Side effects:
        None.
    """

    return (
        f"{len(bounds)} sampled window(s); "
        f"{RAW_QC_SAMPLE_WINDOW_SECONDS:g}s requested at fractions "
        f"{','.join(str(value) for value in RAW_QC_SAMPLE_FRACTIONS)}"
    )


def read_sampled_raw_data(raw: mne.io.BaseRaw, bounds: list[tuple[int, int]]) -> np.ndarray:
    """Read sampled raw signal windows without preloading a full EDF.

    Args:
        raw: MNE Raw object opened with ``preload=False``.
        bounds: Sample bounds to read.

    Returns:
        Array shaped ``n_channels x n_sampled_times``.

    Side effects:
        Reads sampled signal data from the EDF. The Raw object is not modified.
    """

    segments = [raw.get_data(start=start, stop=stop) for start, stop in bounds if stop > start]
    if not segments:
        return np.empty((len(raw.ch_names), 0), dtype=float)
    return np.concatenate(segments, axis=1)


def channel_type_indices(raw: mne.io.BaseRaw) -> dict[str, list[int]]:
    """Group Raw channel indices by post-override channel type.

    Args:
        raw: MNE Raw object after approved in-memory channel-type overrides.

    Returns:
        Mapping from channel type to channel indices.

    Side effects:
        None.
    """

    groups: dict[str, list[int]] = {}
    for index, channel_type in enumerate(raw.get_channel_types()):
        groups.setdefault(channel_type, []).append(index)
    return groups


def aggregated_unit_label(channel_type: str, original_units: list[object]) -> str:
    """Return an honest unit label for an aggregated channel-type row.

    Args:
        channel_type: Post-override MNE channel type.
        original_units: Original EDF dimension values for channels in the
            aggregate.

    Returns:
        Concise unit/semantics text. Mixed classifications are surfaced rather
        than collapsed to a blanket SI claim.

    Side effects:
        None.
    """

    statuses = [classify_unit_status(channel_type, value) for value in original_units]
    status_names = sorted({status.value_status for status in statuses})
    if status_names == ["mne_scaled_calibrated_voltage"]:
        original = sorted({status.edf_original_unit_normalized for status in statuses})
        return f"V from MNE scaling of calibrated EDF {','.join(original)}"
    if status_names == ["arbitrary_unknown_acquisition_units"]:
        return "unknown acquisition units; not calibrated SI"
    if status_names == ["digital_stim_values"]:
        return "digital trigger states; continuous amplitude not interpreted"
    if status_names == ["calibrated_physical_units"]:
        units = sorted({status.in_memory_value_unit for status in statuses})
        return f"declared EDF physical units: {','.join(units)}"
    return f"mixed unit statuses: {','.join(status_names)}"


def base_raw_qc_identity(file_row: dict[str, Any]) -> dict[str, Any]:
    """Return common identity fields for sampled raw-QC outputs.

    Args:
        file_row: File inventory row as a dictionary.

    Returns:
        Participant/file identity fields.

    Side effects:
        None.
    """

    return file_identity_fields(file_row)


def error_amplitude_rows(file_row: dict[str, Any], error_type: str, error_message: str) -> list[dict[str, Any]]:
    """Build an amplitude-summary error row for one file.

    Args:
        file_row: File inventory row as a dictionary.
        error_type: Error class or short status.
        error_message: Human-readable error text.

    Returns:
        One-row list for the amplitude summary table.

    Side effects:
        None.
    """

    return [
        {
            **base_raw_qc_identity(file_row),
            "raw_qc_status": "signal_qc_error",
            "raw_qc_error_type": error_type,
            "raw_qc_error_message": error_message,
            "sample_strategy": "",
            "sample_window_count": 0,
            "sampled_seconds_per_channel": 0.0,
            "channel_type": "",
            "channel_group": "",
            "n_channels": 0,
            "n_samples_per_channel": 0,
            "mne_data_units": "",
            "amplitude_min": pd.NA,
            "amplitude_max": pd.NA,
            "amplitude_peak_to_peak": pd.NA,
            "amplitude_mean": pd.NA,
            "amplitude_median": pd.NA,
            "amplitude_std": pd.NA,
            "amplitude_rms": pd.NA,
            "amplitude_mean_abs": pd.NA,
        }
    ]


def error_psd_rows(file_row: dict[str, Any], error_type: str, error_message: str) -> list[dict[str, Any]]:
    """Build a PSD-summary error row for one file.

    Args:
        file_row: File inventory row as a dictionary.
        error_type: Error class or short status.
        error_message: Human-readable error text.

    Returns:
        One-row list for the PSD summary table.

    Side effects:
        None.
    """

    return [
        {
            **base_raw_qc_identity(file_row),
            "raw_qc_status": "signal_qc_error",
            "raw_qc_error_type": error_type,
            "raw_qc_error_message": error_message,
            "sample_strategy": "",
            "sample_window_count": 0,
            "sampled_seconds_per_channel": 0.0,
            "channel_type": "",
            "channel_group": "",
            "n_channels": 0,
            "welch_nperseg": pd.NA,
            "frequency_min_hz": pd.NA,
            "frequency_max_hz": pd.NA,
            "median_psd": pd.NA,
            "mean_psd": pd.NA,
            "median_psd_1_4_hz": pd.NA,
            "median_psd_4_8_hz": pd.NA,
            "median_psd_8_13_hz": pd.NA,
            "median_psd_13_30_hz": pd.NA,
            "median_psd_30_50_hz": pd.NA,
            "median_psd_55_65_hz": pd.NA,
            "psd_note": "",
        }
    ]


def build_amplitude_rows(
    file_row: dict[str, Any],
    raw: mne.io.BaseRaw,
    sampled_data: np.ndarray,
    bounds: list[tuple[int, int]],
) -> list[dict[str, Any]]:
    """Summarize sampled raw amplitude/range by channel type.

    Args:
        file_row: File inventory row as a dictionary.
        raw: MNE Raw object after approved channel-type overrides.
        sampled_data: Sampled signal array shaped ``n_channels x n_times``.
        bounds: Sample bounds read from the EDF.

    Returns:
        List of amplitude summary rows, one per channel type present.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    groups = channel_type_indices(raw)
    sfreq = float(raw.info["sfreq"])
    seconds = sampled_seconds(bounds, sfreq)
    strategy = sample_strategy_text(bounds)

    for channel_type, indices in sorted(groups.items()):
        group_data = sampled_data[indices, :]
        original_units = [getattr(raw, "_orig_units", {}).get(raw.ch_names[index], "") for index in indices]
        unit_label = aggregated_unit_label(channel_type, original_units)
        finite = group_data[np.isfinite(group_data)]
        if channel_type == "stim" or finite.size == 0:
            stats = {name: pd.NA for name in ("min", "max", "ptp", "mean", "median", "std", "rms", "mean_abs")}
        else:
            stats = {
                "min": float(np.min(finite)),
                "max": float(np.max(finite)),
                "ptp": float(np.ptp(finite)),
                "mean": float(np.mean(finite)),
                "median": float(np.median(finite)),
                "std": float(np.std(finite)),
                "rms": float(np.sqrt(np.mean(np.square(finite)))),
                "mean_abs": float(np.mean(np.abs(finite))),
            }

        rows.append(
            {
                **base_raw_qc_identity(file_row),
                "raw_qc_status": "ok",
                "raw_qc_error_type": "",
                "raw_qc_error_message": "",
                "sample_strategy": strategy,
                "sample_window_count": len(bounds),
                "sampled_seconds_per_channel": seconds,
                "channel_type": channel_type,
                "channel_group": channel_group_for_type(channel_type),
                "n_channels": len(indices),
                "n_samples_per_channel": int(group_data.shape[1]),
                "mne_data_units": unit_label,
                "amplitude_min": stats["min"],
                "amplitude_max": stats["max"],
                "amplitude_peak_to_peak": stats["ptp"],
                "amplitude_mean": stats["mean"],
                "amplitude_median": stats["median"],
                "amplitude_std": stats["std"],
                "amplitude_rms": stats["rms"],
                "amplitude_mean_abs": stats["mean_abs"],
            }
        )

    return rows


def median_psd_in_band(frequencies: np.ndarray, psd_values: np.ndarray, low_hz: float, high_hz: float) -> float | Any:
    """Return median PSD within a frequency band.

    Args:
        frequencies: Welch frequency bins.
        psd_values: PSD values shaped ``n_channels x n_frequencies``.
        low_hz: Inclusive lower frequency bound.
        high_hz: Inclusive upper frequency bound.

    Returns:
        Median PSD as a float, or ``pd.NA`` when the band has no bins.

    Side effects:
        None.
    """

    mask = (frequencies >= low_hz) & (frequencies <= high_hz)
    if not np.any(mask):
        return pd.NA
    values = psd_values[:, mask]
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return pd.NA
    return float(np.median(finite))


def build_psd_rows(
    file_row: dict[str, Any],
    raw: mne.io.BaseRaw,
    sampled_data: np.ndarray,
    bounds: list[tuple[int, int]],
) -> list[dict[str, Any]]:
    """Summarize sampled raw PSD by channel type.

    Args:
        file_row: File inventory row as a dictionary.
        raw: MNE Raw object after approved channel-type overrides.
        sampled_data: Sampled signal array shaped ``n_channels x n_times``.
        bounds: Sample bounds read from the EDF.

    Returns:
        List of PSD summary rows. Stim channels receive a factual skipped row
        because trigger pulses are not continuous EEG/auxiliary signals.

    Side effects:
        Computes Welch PSD summaries from sampled raw data. It does not modify
        the Raw object and does not filter the signal.
    """

    rows: list[dict[str, Any]] = []
    groups = channel_type_indices(raw)
    sfreq = float(raw.info["sfreq"])
    seconds = sampled_seconds(bounds, sfreq)
    strategy = sample_strategy_text(bounds)
    nperseg = max(1, min(int(round(RAW_QC_PSD_NPERSEG_SECONDS * sfreq)), sampled_data.shape[1]))

    for channel_type, indices in sorted(groups.items()):
        base = {
            **base_raw_qc_identity(file_row),
            "raw_qc_status": "ok",
            "raw_qc_error_type": "",
            "raw_qc_error_message": "",
            "sample_strategy": strategy,
            "sample_window_count": len(bounds),
            "sampled_seconds_per_channel": seconds,
            "channel_type": channel_type,
            "channel_group": channel_group_for_type(channel_type),
            "n_channels": len(indices),
            "welch_nperseg": nperseg,
        }
        if channel_type == "stim":
            rows.append(
                {
                    **base,
                    "frequency_min_hz": pd.NA,
                    "frequency_max_hz": pd.NA,
                    "median_psd": pd.NA,
                    "mean_psd": pd.NA,
                    "median_psd_1_4_hz": pd.NA,
                    "median_psd_4_8_hz": pd.NA,
                    "median_psd_8_13_hz": pd.NA,
                    "median_psd_13_30_hz": pd.NA,
                    "median_psd_30_50_hz": pd.NA,
                    "median_psd_55_65_hz": pd.NA,
                    "psd_note": "stim_channel_group_not_summarized_as_continuous_psd",
                }
            )
            continue

        group_data = sampled_data[indices, :]
        frequencies, psd_values = welch(group_data, fs=sfreq, nperseg=nperseg, axis=1)
        frequency_mask = frequencies <= RAW_QC_PSD_MAX_FREQUENCY_HZ
        frequencies = frequencies[frequency_mask]
        psd_values = psd_values[:, frequency_mask]
        finite = psd_values[np.isfinite(psd_values)]

        rows.append(
            {
                **base,
                "frequency_min_hz": float(frequencies.min()) if frequencies.size else pd.NA,
                "frequency_max_hz": float(frequencies.max()) if frequencies.size else pd.NA,
                "median_psd": float(np.median(finite)) if finite.size else pd.NA,
                "mean_psd": float(np.mean(finite)) if finite.size else pd.NA,
                "median_psd_1_4_hz": median_psd_in_band(frequencies, psd_values, 1.0, 4.0),
                "median_psd_4_8_hz": median_psd_in_band(frequencies, psd_values, 4.0, 8.0),
                "median_psd_8_13_hz": median_psd_in_band(frequencies, psd_values, 8.0, 13.0),
                "median_psd_13_30_hz": median_psd_in_band(frequencies, psd_values, 13.0, 30.0),
                "median_psd_30_50_hz": median_psd_in_band(frequencies, psd_values, 30.0, 50.0),
                "median_psd_55_65_hz": median_psd_in_band(frequencies, psd_values, 55.0, 65.0),
                "psd_note": "sampled_raw_welch_psd_no_filtering",
            }
        )

    return rows


def safe_figure_stem(file_row: dict[str, Any]) -> str:
    """Return a filesystem-safe stem for raw-QC figure filenames.

    Args:
        file_row: File inventory row as a dictionary.

    Returns:
        Safe filename stem based on source EDF name.

    Side effects:
        None.
    """

    source_name = Path(text_or_empty(file_row.get("source_filename"))).stem
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", source_name).strip("_") or "unknown_file"


def write_time_series_snapshot_figure(
    figure_path: Path,
    raw: mne.io.BaseRaw,
    sampled_data: np.ndarray,
    bounds: list[tuple[int, int]],
) -> None:
    """Write a compact sampled raw time-series snapshot figure.

    Args:
        figure_path: Destination PNG path.
        raw: MNE Raw object after channel-type overrides.
        sampled_data: Sampled signal array.
        bounds: Sample bounds read from the EDF.

    Returns:
        ``None``.

    Side effects:
        Writes a PNG figure. Data are standardized per channel for display only
        so EEG/EOG/EMG/stim traces can be seen together without changing any
        saved signal object.
    """

    figure_path.parent.mkdir(parents=True, exist_ok=True)
    preferred_channels = ["FZ", "CZ", "PZ", "OZ", "HEO", "VEO", "EMG-L", "EMG-A", "Trigger"]
    channel_indices = [raw.ch_names.index(name) for name in preferred_channels if name in raw.ch_names]
    if not channel_indices:
        channel_indices = list(range(min(8, len(raw.ch_names))))

    sfreq = float(raw.info["sfreq"])
    max_samples = min(sampled_data.shape[1], int(round(8 * sfreq)))
    plot_data = sampled_data[channel_indices, :max_samples].copy()
    times = np.arange(max_samples) / sfreq

    # Standardize only the plotted copy. This avoids visually losing small EEG
    # channels next to trigger pulses or auxiliary channels with different
    # units, while leaving raw QC tables in MNE-returned units.
    centered = plot_data - np.nanmedian(plot_data, axis=1, keepdims=True)
    scale = np.nanstd(centered, axis=1, keepdims=True)
    scale[~np.isfinite(scale) | (scale == 0)] = 1.0
    display_data = centered / scale

    fig, ax = plt.subplots(figsize=(11, 6))
    offsets = np.arange(len(channel_indices))[::-1] * 4.0
    for row_index, channel_index in enumerate(channel_indices):
        ax.plot(times, display_data[row_index] + offsets[row_index], linewidth=0.8)
    ax.set_yticks(offsets)
    ax.set_yticklabels([raw.ch_names[index] for index in channel_indices])
    ax.set_xlabel("Seconds from first sampled QC window")
    ax.set_title("Raw sampled time-series snapshot (standardized for display)")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=RAW_QC_FIGURE_DPI)
    plt.close(fig)


def write_psd_figure(
    figure_path: Path,
    raw: mne.io.BaseRaw,
    sampled_data: np.ndarray,
) -> None:
    """Write a compact sampled raw PSD figure by channel type.

    Args:
        figure_path: Destination PNG path.
        raw: MNE Raw object after channel-type overrides.
        sampled_data: Sampled signal array.

    Returns:
        ``None``.

    Side effects:
        Writes a PNG figure from sampled raw data. The signal is not filtered
        and no Raw object is saved.
    """

    figure_path.parent.mkdir(parents=True, exist_ok=True)
    sfreq = float(raw.info["sfreq"])
    nperseg = max(1, min(int(round(RAW_QC_PSD_NPERSEG_SECONDS * sfreq)), sampled_data.shape[1]))
    groups = channel_type_indices(raw)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    plotted = False
    for channel_type, indices in sorted(groups.items()):
        if channel_type == "stim":
            continue
        frequencies, psd_values = welch(sampled_data[indices, :], fs=sfreq, nperseg=nperseg, axis=1)
        mask = (frequencies > 0) & (frequencies <= RAW_QC_PSD_MAX_FREQUENCY_HZ)
        if not np.any(mask):
            continue
        median_psd = np.nanmedian(psd_values[:, mask], axis=0)
        ax.plot(frequencies[mask], 10 * np.log10(np.maximum(median_psd, np.finfo(float).tiny)), label=channel_type)
        plotted = True

    if not plotted:
        ax.text(0.5, 0.5, "No continuous channel groups available for PSD", ha="center", va="center")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Median PSD (dB, sampled raw data)")
    ax.set_title("Raw sampled PSD by channel type")
    ax.grid(True, alpha=0.25)
    if plotted:
        ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path, dpi=RAW_QC_FIGURE_DPI)
    plt.close(fig)


def figure_manifest_error_rows(file_row: dict[str, Any], error_type: str, error_message: str) -> list[dict[str, Any]]:
    """Build figure-manifest error rows for one file.

    Args:
        file_row: File inventory row as a dictionary.
        error_type: Error class or short status.
        error_message: Human-readable error text.

    Returns:
        Error rows for the figure manifest.

    Side effects:
        None.
    """

    return [
        {
            **base_raw_qc_identity(file_row),
            "figure_type": "raw_signal_qc",
            "figure_path": "",
            "raw_qc_status": "signal_qc_error",
            "raw_qc_error_type": error_type,
            "raw_qc_error_message": error_message,
            "figure_note": "",
        }
    ]


def run_sampled_signal_qc(
    inventory: pd.DataFrame,
    overrides_by_normalized_name: dict[str, str],
    output_dir: Path,
    repo_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run sampled raw-signal QC one EDF at a time.

    Args:
        inventory: File inventory from manifests and raw-directory discovery.
        overrides_by_normalized_name: Approved in-memory channel-type override
            mapping.
        output_dir: Raw-QC output directory.
        repo_root: Repository root for rendering relative figure paths.

    Returns:
        ``(amplitude_summary, psd_summary, figure_manifest)``.

    Side effects:
        Opens each EDF with ``preload=False``, reads sampled raw data windows,
        writes PNG figures below ``output_dir / figures``, and closes each Raw
        object before moving to the next file. No signal derivatives are saved.
    """

    amplitude_rows: list[dict[str, Any]] = []
    psd_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []

    time_series_dir = output_dir / FIGURE_DIRNAME / TIME_SERIES_FIGURE_DIRNAME
    psd_dir = output_dir / FIGURE_DIRNAME / PSD_FIGURE_DIRNAME

    for _, inventory_row in inventory.iterrows():
        file_row = inventory_row.to_dict()
        if not file_row["file_exists"]:
            message = "raw file missing"
            amplitude_rows.extend(error_amplitude_rows(file_row, "raw_file_missing", message))
            psd_rows.extend(error_psd_rows(file_row, "raw_file_missing", message))
            figure_rows.extend(figure_manifest_error_rows(file_row, "raw_file_missing", message))
            continue

        raw: mne.io.BaseRaw | None = None
        try:
            raw = mne.io.read_raw_edf(Path(file_row["resolved_edf_path"]), preload=False, verbose="ERROR")
            apply_configured_channel_type_overrides(raw, overrides_by_normalized_name)
            sfreq = float(raw.info["sfreq"])
            bounds = sample_window_bounds(raw.n_times, sfreq)
            sampled_data = read_sampled_raw_data(raw, bounds)

            amplitude_rows.extend(build_amplitude_rows(file_row, raw, sampled_data, bounds))
            psd_rows.extend(build_psd_rows(file_row, raw, sampled_data, bounds))

            stem = safe_figure_stem(file_row)
            time_series_path = time_series_dir / f"{stem}_raw_timeseries_snapshot.png"
            psd_path = psd_dir / f"{stem}_raw_psd.png"
            write_time_series_snapshot_figure(time_series_path, raw, sampled_data, bounds)
            write_psd_figure(psd_path, raw, sampled_data)

            figure_rows.extend(
                [
                    {
                        **base_raw_qc_identity(file_row),
                        "figure_type": "raw_time_series_snapshot",
                        "figure_path": relative_to_repo(time_series_path, repo_root),
                        "raw_qc_status": "ok",
                        "raw_qc_error_type": "",
                        "raw_qc_error_message": "",
                        "figure_note": "sampled_raw_data_standardized_for_display",
                    },
                    {
                        **base_raw_qc_identity(file_row),
                        "figure_type": "raw_psd",
                        "figure_path": relative_to_repo(psd_path, repo_root),
                        "raw_qc_status": "ok",
                        "raw_qc_error_type": "",
                        "raw_qc_error_message": "",
                        "figure_note": "sampled_raw_welch_psd_no_filtering",
                    },
                ]
            )
        except Exception as error:  # noqa: BLE001 - one problematic EDF should not stop raw QC.
            error_type = type(error).__name__
            error_message = str(error)
            amplitude_rows.extend(error_amplitude_rows(file_row, error_type, error_message))
            psd_rows.extend(error_psd_rows(file_row, error_type, error_message))
            figure_rows.extend(figure_manifest_error_rows(file_row, error_type, error_message))
        finally:
            if raw is not None and hasattr(raw, "close"):
                raw.close()

    return (
        pd.DataFrame(amplitude_rows, columns=AMPLITUDE_RANGE_COLUMNS),
        pd.DataFrame(psd_rows, columns=PSD_SUMMARY_COLUMNS),
        pd.DataFrame(figure_rows, columns=FIGURE_MANIFEST_COLUMNS),
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


def filter_inventory_for_config(inventory: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Apply config file-role switches to the raw-QC inventory.

    Args:
        inventory: Inventory built from the raw manifest and raw directory.
        config: Parsed raw-QC config.

    Returns:
        Inventory rows whose file role is enabled for raw QC.

    Side effects:
        None. This is a config-scoped processing selection, not a participant
        decision.
    """

    role_switches = require_config_value(config, ("file_classes", "process_for_raw_qc"))
    if not isinstance(role_switches, dict):
        raise RuntimeError("file_classes.process_for_raw_qc must be a mapping")
    enabled_roles = {role for role, enabled in role_switches.items() if bool(enabled)}
    return inventory[inventory["file_role"].isin(enabled_roles)].reset_index(drop=True)


def build_run_manifest(
    config: dict[str, Any],
    config_path: Path,
    output_dir: Path,
    file_summary: pd.DataFrame,
    amplitude_summary: pd.DataFrame,
    standard_1005_status: str,
    repo_root: Path,
) -> dict[str, Any]:
    """Build a machine-readable raw-QC run manifest.

    Args:
        config: Parsed raw-QC config.
        config_path: Config path used for this run.
        output_dir: Raw-QC output directory.
        file_summary: File-level raw-QC summary table.
        amplitude_summary: Sampled amplitude/range table.
        standard_1005_status: Validation-only standard montage check status.
        repo_root: Repository root.

    Returns:
        JSON-serializable run manifest dictionary.

    Side effects:
        None.
    """

    return {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "script": relative_to_repo(Path(__file__), repo_root),
        "config_path": relative_to_repo(config_path, repo_root),
        "run_label": config_value(config, ("run", "label"), ""),
        "output_root": relative_to_repo(output_dir, repo_root),
        "file_summary_rows": int(len(file_summary)),
        "read_status_counts": {
            str(status): int(count)
            for status, count in file_summary["read_status"].value_counts(dropna=False).sort_index().items()
        },
        "raw_signal_qc_status_counts": {
            str(status): int(count)
            for status, count in amplitude_summary["raw_qc_status"].value_counts(dropna=False).sort_index().items()
        }
        if not amplitude_summary.empty
        else {},
        "standard_1005_validation_status": standard_1005_status,
        "raw_qc_only": True,
        "disabled_stages_confirmed_off": [
            "filtering",
            "notch_filtering",
            "reference",
            "bad_channel_detection",
            "interpolation",
            "ica",
            "csd",
            "emg_movement_summary",
            "emg_derived_trial_handling",
            "event_repair",
            "epoching",
            "cleaned_fif_writing",
        ],
    }


def write_summary_markdown(
    output_path: Path,
    config_path: Path,
    run_manifest: dict[str, Any],
    file_summary: pd.DataFrame,
    channel_validation: pd.DataFrame,
    montage_validation: pd.DataFrame,
    event_context: pd.DataFrame,
    amplitude_summary: pd.DataFrame,
    psd_summary: pd.DataFrame,
    figure_manifest: pd.DataFrame,
    optional_context: dict[str, pd.DataFrame | None],
) -> None:
    """Write the raw-QC Markdown summary.

    Args:
        output_path: Destination Markdown file.
        config_path: Config path used for the run.
        run_manifest: Run manifest dictionary written as JSON.
        file_summary: File-level raw-QC summary table.
        channel_validation: Channel-type validation table.
        montage_validation: Montage coordinate validation table.
        event_context: Event-context flag table.
        amplitude_summary: Sampled amplitude/range summary table.
        psd_summary: Sampled PSD summary table.
        figure_manifest: Figure manifest table.
        optional_context: Optional event-context tables.

    Returns:
        ``None``.

    Side effects:
        Writes ``output_path``.
    """

    issue_counts = count_semicolon_codes(file_summary["raw_qc_issue_codes"])
    event_code_counts = count_semicolon_codes(file_summary["event_context_pending_codes"])
    role_counts = file_summary["file_role"].value_counts(dropna=False).sort_index()
    read_status_counts = file_summary["read_status"].value_counts(dropna=False).sort_index()
    override_status_counts = file_summary["channel_type_override_status"].value_counts(dropna=False).sort_index()

    read_failures = file_summary[file_summary["read_status"].ne("ok")]
    channel_mismatch_files = file_summary[
        file_summary["missing_expected_channel_names"].fillna("").ne("")
        | file_summary["unexpected_channel_names"].fillna("").ne("")
        | ~file_summary["expected_37_channels_present"].fillna(False)
    ]
    montage_problem_files = file_summary[file_summary["besa_unmatched_eeg_channel_names"].fillna("").ne("")]
    standard_problem_files = file_summary[file_summary["standard_1005_unmatched_eeg_channel_names"].fillna("").ne("")]
    zero_annotation_files = file_summary[file_summary["zero_annotation_fact"].fillna(False)]
    low_annotation_files = file_summary[file_summary["low_annotation_fact"].fillna(False)]
    signal_qc_errors = amplitude_summary[amplitude_summary["raw_qc_status"].ne("ok")] if not amplitude_summary.empty else amplitude_summary

    raw_only_categories = Counter()
    raw_only = optional_context.get("raw_only_row_classification")
    if raw_only is not None and "raw_only_category" in raw_only.columns:
        raw_only_categories.update(raw_only["raw_only_category"].dropna().map(text_or_empty))

    lines = [
        "# Continuous Raw EEG QC Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Config: `{relative_to_repo(config_path, repo_root_from_script())}`",
        f"Run label: `{run_manifest.get('run_label', '')}`",
        "",
        "This was a raw-QC/config-validation run only. No filtering, notch filtering, re-reference, bad-channel detection, interpolation, ICA, CSD, EMG movement-summary computation, EMG-derived trial handling, event repair, epoching, cleaned FIF writing, or participant-level decision was performed.",
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

    lines.extend(["", "## Channel type overrides", ""])
    for status, count in override_status_counts.items():
        lines.append(f"- {status}: {int(count)}")

    lines.extend(["", "## Raw-QC issue counts", ""])
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
            f"- Files with unmatched validation-only standard_1005 names: {len(standard_problem_files)}",
            f"- Channel-type validation rows: {len(channel_validation)}",
            f"- Montage validation rows: {len(montage_validation)}",
        ]
    )
    if not channel_mismatch_files.empty:
        lines.append(f"- Channel-label follow-up files: {format_source_filenames(channel_mismatch_files)}")
    if not montage_problem_files.empty:
        lines.append(f"- BESA-coordinate follow-up files: {format_source_filenames(montage_problem_files)}")
    if not standard_problem_files.empty:
        lines.append(f"- standard_1005 name-check follow-up files: {format_source_filenames(standard_problem_files)}")

    lines.extend(
        [
            "",
            "## Sampled raw signal QC",
            "",
            f"- Amplitude/range summary rows: {len(amplitude_summary)}",
            f"- PSD summary rows: {len(psd_summary)}",
            f"- Figure manifest rows: {len(figure_manifest)}",
            f"- Signal-QC error files: {format_source_filenames(signal_qc_errors)}",
            "",
            "Signal QC used short sampled windows from each EDF, one file at a time. These summaries are raw QC only and are not cleaned derivatives.",
        ]
    )

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
            f"- `{AMPLITUDE_RANGE_FILENAME}`",
            f"- `{PSD_SUMMARY_FILENAME}`",
            f"- `{FIGURE_MANIFEST_FILENAME}`",
            f"- `{RUN_MANIFEST_FILENAME}`",
            "",
            "These outputs are local-only under `_Data/`.",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(
    output_dir: Path,
    config_path: Path,
    run_manifest: dict[str, Any],
    file_summary: pd.DataFrame,
    channel_validation: pd.DataFrame,
    montage_validation: pd.DataFrame,
    event_context: pd.DataFrame,
    amplitude_summary: pd.DataFrame,
    psd_summary: pd.DataFrame,
    figure_manifest: pd.DataFrame,
    optional_context: dict[str, pd.DataFrame | None],
) -> None:
    """Write all raw-QC output tables, figures manifest, and summary.

    Args:
        output_dir: Raw-QC output directory.
        config_path: Config path used for the run.
        run_manifest: Run manifest dictionary.
        file_summary: File-level raw-QC summary table.
        channel_validation: Channel-type validation table.
        montage_validation: Montage coordinate validation table.
        event_context: Event-context flag table.
        amplitude_summary: Sampled amplitude/range summary table.
        psd_summary: Sampled PSD summary table.
        figure_manifest: Figure manifest table.
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
    amplitude_summary = coerce_split_part_for_output(amplitude_summary)
    psd_summary = coerce_split_part_for_output(psd_summary)
    figure_manifest = coerce_split_part_for_output(figure_manifest)

    (output_dir / RUN_MANIFEST_FILENAME).write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    file_summary.to_csv(output_dir / FILE_SUMMARY_FILENAME, index=False)
    channel_validation.to_csv(output_dir / CHANNEL_TYPE_VALIDATION_FILENAME, index=False)
    montage_validation.to_csv(output_dir / MONTAGE_COORDINATE_VALIDATION_FILENAME, index=False)
    event_context.to_csv(output_dir / EVENT_CONTEXT_FLAGS_FILENAME, index=False)
    amplitude_summary.to_csv(output_dir / AMPLITUDE_RANGE_FILENAME, index=False)
    psd_summary.to_csv(output_dir / PSD_SUMMARY_FILENAME, index=False)
    figure_manifest.to_csv(output_dir / FIGURE_MANIFEST_FILENAME, index=False)
    write_summary_markdown(
        output_dir / RAW_QC_SUMMARY_FILENAME,
        config_path,
        run_manifest,
        file_summary,
        channel_validation,
        montage_validation,
        event_context,
        amplitude_summary,
        psd_summary,
        figure_manifest,
        optional_context,
    )


def main() -> None:
    """Run the config-driven continuous raw EEG QC stage.

    Args:
        None. Command-line arguments are parsed by ``parse_args``.

    Returns:
        ``None``.

    Side effects:
        Reads the private raw-QC config, required local manifests, BESA
        coordinates, optional event-audit context tables, EDF headers with MNE
        ``preload=False``, and sampled raw signal windows for QC. Writes
        local-only QC outputs under the configured ``_Data/`` output root.
    """

    args = parse_args()
    repo_root = repo_root_from_script()
    config_path = resolve_repo_path(repo_root, args.config)
    config = load_raw_qc_config(config_path)
    validate_raw_qc_config(config, config_path, repo_root)

    raw_dir = resolve_repo_path(repo_root, config_value(config, ("inputs", "raw_edf_root"), RAW_EEG_DIR))
    output_dir = resolve_repo_path(repo_root, config_value(config, ("outputs", "output_root"), RAW_QC_OUTPUT_DIR))
    overrides_by_normalized_name = configured_overrides_by_normalized_name(config)
    standard_1005_by_normalized_name, standard_1005_status = standard_1005_name_map_if_requested(config)

    manifest = read_csv_checked(
        resolve_repo_path(repo_root, config_value(config, ("inputs", "raw_manifest"), RAW_MANIFEST_PATH)),
        {"participant_id", "source_filename", "file_path", "file_role", "annotation_count"},
        "raw EEG file manifest",
    )
    annotation_counts = read_csv_checked(
        resolve_repo_path(repo_root, config_value(config, ("inputs", "raw_annotation_counts"), ANNOTATION_COUNTS_PATH)),
        {"source_filename", "annotation_description", "annotation_count"},
        "raw EEG annotation-count CSV",
    )
    besa_coordinates = load_besa_coordinates(
        resolve_repo_path(repo_root, config_value(config, ("inputs", "besa_81_coordinates"), BESA_COORDINATES_PATH))
    )
    optional_context = read_optional_event_context(repo_root, optional_event_context_paths(config))

    edf_files = discover_edf_files(raw_dir)
    annotation_totals = aggregate_annotation_counts(annotation_counts)
    inventory = build_file_inventory(manifest, annotation_totals, edf_files, raw_dir, repo_root)
    inventory = filter_inventory_for_config(inventory, config)

    expected_file_count = safe_int(config_value(config, ("expected_raw_metadata", "expected_file_count")))
    if expected_file_count is not None and len(inventory) != expected_file_count:
        raise RuntimeError(
            "raw-QC inventory row count does not match expected_raw_metadata.expected_file_count: "
            f"observed={len(inventory)}, expected={expected_file_count}"
        )

    file_summary, channel_validation, montage_validation = run_header_validations(
        inventory,
        besa_coordinates,
        overrides_by_normalized_name,
        standard_1005_by_normalized_name,
    )
    event_context = build_event_context_flags(file_summary, optional_context)
    amplitude_summary, psd_summary, figure_manifest = run_sampled_signal_qc(
        inventory,
        overrides_by_normalized_name,
        output_dir,
        repo_root,
    )
    run_manifest = build_run_manifest(
        config,
        config_path,
        output_dir,
        file_summary,
        amplitude_summary,
        standard_1005_status,
        repo_root,
    )

    write_outputs(
        output_dir,
        config_path,
        run_manifest,
        file_summary,
        channel_validation,
        montage_validation,
        event_context,
        amplitude_summary,
        psd_summary,
        figure_manifest,
        optional_context,
    )

    print(f"Raw QC inspected {len(file_summary)} EDF file row(s).")
    read_status_counts = {
        str(status): int(count)
        for status, count in file_summary["read_status"].value_counts(dropna=False).sort_index().items()
    }
    print(f"MNE read status counts: {read_status_counts}")
    signal_qc_counts = {
        str(status): int(count)
        for status, count in amplitude_summary["raw_qc_status"].value_counts(dropna=False).sort_index().items()
    }
    print(f"Sampled raw signal QC status counts: {signal_qc_counts}")
    print(f"Wrote raw-QC outputs to {output_dir}")


if __name__ == "__main__":
    main()
