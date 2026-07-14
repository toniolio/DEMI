"""Versioned configuration and channel contracts for continuous preprocessing.

This module centralizes the accepted DEMI channel surface and fail-closed
validation of the tracked production configuration. It exists so scientific
parameters, safety boundaries, and derivative retention choices are explicit
and testable rather than scattered through a command-line script.

Inputs are the tracked YAML configuration and ordinary paths. Returned values
are immutable Python mappings and tuples used by the production stages.
Loading the contract reads one YAML file and writes nothing.

This module does not read signal samples, mutate EDF files, construct epochs,
run AutoReject, compute CSD, or make participant inclusion decisions.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import yaml


CONFIG_SCHEMA_VERSION = 1
PIPELINE_VERSION = "continuous_preprocessing_v1"

SCALP_SOURCE_CHANNELS = (
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
    "TP7",
    "CP3",
    "CPZ",
    "CP4",
    "TP8",
    "P7",
    "P3",
    "PZ",
    "P4",
    "P8",
    "O1",
    "OZ",
    "O2",
)
MASTOID_CHANNELS = ("M1", "M2")

# This order matches the physical EDF channel order and remains the saved EEG
# target order. M1/M2 are intentionally retained where they were recorded.
EEG_TARGET_CHANNELS = (
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
EOG_CHANNELS = ("HEO", "VEO")
EMG_CHANNELS = ("EMG-L", "EMG-A")
STIM_CHANNELS = ("Trigger",)
ALL_RETAINED_CHANNELS = EEG_TARGET_CHANNELS + EOG_CHANNELS + EMG_CHANNELS + STIM_CHANNELS

CHANNEL_TYPE_BY_NAME = MappingProxyType(
    {
        **{channel: "eeg" for channel in EEG_TARGET_CHANNELS},
        **{channel: "eog" for channel in EOG_CHANNELS},
        **{channel: "emg" for channel in EMG_CHANNELS},
        **{channel: "stim" for channel in STIM_CHANNELS},
    }
)

PRIMARY_PYPREP_CRITERIA = (
    "bad_by_nan",
    "bad_by_flat",
    "bad_by_deviation",
    "bad_by_hf_noise",
    "bad_by_correlation",
    "bad_by_SNR",
    "bad_by_dropout",
    "bad_by_ransac",
)
REPORT_ONLY_PYPREP_CRITERIA = ("bad_by_psd",)

STAGE_NAMES = (
    "source_validation_and_read",
    "channel_typing",
    "montage_application",
    "optional_line_noise_removal",
    "analysis_filtering",
    "detector_input_preparation",
    "pyprep_criterion_detection",
    "accepted_global_bad_union",
    "reference_estimation",
    "reference_application",
    "interpolation",
    "post_interpolation_validation",
    "ica_fitting_copy_preparation",
    "rank_estimation",
    "ica_fit",
    "component_scoring",
    "automatic_component_proposal",
    "exception_or_continuation_decision",
    "derivative_write",
    "manifest_finalization",
)


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize JSON-compatible data deterministically.

    Args:
        value: JSON-compatible object.

    Returns:
        UTF-8 bytes with sorted keys and no insignificant whitespace.

    Side effects:
        None.
    """

    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    """Return the hexadecimal SHA-256 digest of bytes without side effects."""

    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    """Hash one file through bounded reads.

    Args:
        path: File to hash.
        chunk_bytes: Maximum bytes read for each update.

    Returns:
        Hexadecimal SHA-256 digest.

    Side effects:
        Opens ``path`` read-only.
    """

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def _require_equal(actual: Any, expected: Any, label: str) -> None:
    """Raise a configuration error when a fixed accepted value changed."""

    if actual != expected:
        raise ValueError(f"{label} must be {expected!r}; found {actual!r}.")


def _require_channel_list(config: Mapping[str, Any], key: str, expected: Sequence[str]) -> None:
    """Validate one ordered channel list from the configuration."""

    actual = config.get("channels", {}).get(key)
    _require_equal(actual, list(expected), f"channels.{key}")


def validate_config(config: Mapping[str, Any]) -> None:
    """Fail closed unless the configuration matches the accepted v1 contract.

    Args:
        config: Parsed YAML mapping.

    Returns:
        None.

    Raises:
        ValueError: If a fixed parameter, channel surface, storage boundary, or
            prohibited operation differs from the accepted contract.

    Side effects:
        None.
    """

    _require_equal(config.get("schema_version"), CONFIG_SCHEMA_VERSION, "schema_version")
    _require_equal(config.get("pipeline_version"), PIPELINE_VERSION, "pipeline_version")

    safety = config.get("safety", {})
    for key in (
        "mutate_source_edf",
        "write_edf",
        "construct_epochs",
        "run_autoreject",
        "compute_csd",
        "change_event_eligibility",
        "participant_inclusion_decisions",
        "parallel_recordings",
    ):
        _require_equal(safety.get(key), False, f"safety.{key}")
    _require_equal(safety.get("process_one_recording_at_a_time"), True, "safety.process_one_recording_at_a_time")

    _require_channel_list(config, "detector_source_channels", SCALP_SOURCE_CHANNELS)
    _require_channel_list(config, "reference_source_channels", SCALP_SOURCE_CHANNELS)
    _require_channel_list(config, "eeg_target_channels", EEG_TARGET_CHANNELS)
    _require_channel_list(config, "mastoid_source_excluded_target_retained", MASTOID_CHANNELS)
    _require_channel_list(config, "eog_channels", EOG_CHANNELS)
    _require_channel_list(config, "emg_channels", EMG_CHANNELS)
    _require_channel_list(config, "stim_channels", STIM_CHANNELS)
    _require_channel_list(config, "all_retained_channels", ALL_RETAINED_CHANNELS)

    _require_equal(config.get("montage", {}).get("name"), "standard_1005", "montage.name")
    _require_equal(config.get("montage", {}).get("historical_besa_is_active"), False, "montage.historical_besa_is_active")

    line = config.get("line_noise", {})
    if not isinstance(line.get("enabled"), bool):
        raise ValueError("line_noise.enabled must be one fixed run-level boolean.")
    _require_equal(float(line.get("frequency_hz", 0.0)), 60.0, "line_noise.frequency_hz")
    _require_equal(line.get("method"), "spectrum_fit", "line_noise.method")
    _require_equal(line.get("filter_length"), "10s", "line_noise.filter_length")
    _require_equal(float(line.get("multitaper_bandwidth_hz", 0.0)), 2.0, "line_noise.multitaper_bandwidth_hz")
    _require_equal(float(line.get("p_value", 0.0)), 0.01, "line_noise.p_value")

    analysis_filter = config.get("analysis_filter", {})
    expected_analysis = {
        "high_pass_hz": 0.5,
        "low_pass_hz": 45.0,
        "method": "fir",
        "phase": "zero",
        "fir_design": "firwin",
        "pad": "reflect_limited",
        "skip_by_annotation": ["edge", "bad_acq_skip"],
    }
    for key, expected in expected_analysis.items():
        _require_equal(analysis_filter.get(key), expected, f"analysis_filter.{key}")

    detector = config.get("detector", {})
    _require_equal(detector.get("input_branch"), "native_unfiltered_unnotched_unreferenced", "detector.input_branch")
    _require_equal(detector.get("resample_hz"), None, "detector.resample_hz")
    _require_equal(detector.get("random_seed"), 20260712, "detector.random_seed")
    for key in ("do_detrend", "correlation", "ransac"):
        _require_equal(detector.get(key), True, f"detector.{key}")
    _require_equal(detector.get("determinism_repeats"), 2, "detector.determinism_repeats")
    _require_equal(detector.get("primary_criteria"), list(PRIMARY_PYPREP_CRITERIA), "detector.primary_criteria")
    _require_equal(detector.get("report_only_criteria"), list(REPORT_ONLY_PYPREP_CRITERIA), "detector.report_only_criteria")

    reference = config.get("reference", {})
    _require_equal(reference.get("kind"), "explicit_scalp_average", "reference.kind")
    _require_equal(reference.get("minimum_usable_scalp_channels"), 23, "reference.minimum_usable_scalp_channels")
    _require_equal(reference.get("apply_to_all_32_eeg_targets"), True, "reference.apply_to_all_32_eeg_targets")

    interpolation = config.get("interpolation", {})
    _require_equal(interpolation.get("method"), "spherical_spline", "interpolation.method")
    _require_equal(interpolation.get("mastoid_candidates"), False, "interpolation.mastoid_candidates")
    _require_equal(float(interpolation.get("stop_proportion", 0.0)), 0.25, "interpolation.stop_proportion")
    _require_equal(interpolation.get("proportion_denominator"), 30, "interpolation.proportion_denominator")
    _require_equal(interpolation.get("reset_bads_after_success"), True, "interpolation.reset_bads_after_success")

    ica = config.get("ica", {})
    expected_ica = {
        "high_pass_hz": 1.0,
        "low_pass_hz": 45.0,
        "method": "infomax",
        "extended": True,
        "random_seed": 20260712,
        "n_components": "estimated_eeg_rank",
        "rank_tolerance": 1.0e-6,
        "rank_tolerance_kind": "relative",
        "eog_channels": list(EOG_CHANNELS),
        "eog_measure": "zscore",
        "eog_threshold": 3.0,
        "eog_score_low_hz": 1.0,
        "eog_score_high_hz": 10.0,
        "component_selection_rule": "historical_all_find_bads_eog_components",
        "ordinary_zero_component_action": "stop_as_historical_pipeline_error",
        "component_review_exception_ids": [86],
    }
    for key, expected in expected_ica.items():
        _require_equal(ica.get(key), expected, f"ica.{key}")
    fit_decim = int(ica.get("fit_decim", 0))
    if fit_decim < 1:
        raise ValueError("ica.fit_decim must be a positive fixed integer.")

    derivatives = config.get("derivatives", {})
    _require_equal(derivatives.get("fif_format"), "single", "derivatives.fif_format")
    _require_equal(derivatives.get("ordinary_signal_derivative"), "post_ica_only", "derivatives.ordinary_signal_derivative")
    _require_equal(derivatives.get("retain_ica_object"), True, "derivatives.retain_ica_object")
    _require_equal(derivatives.get("retain_pre_ica_for_review_stops"), True, "derivatives.retain_pre_ica_for_review_stops")
    _require_equal(derivatives.get("atomic_directory_publish"), True, "derivatives.atomic_directory_publish")

    output_root = str(config.get("paths", {}).get("output_root", ""))
    _require_equal(
        output_root,
        "_Data/eeg/mne_preprocessing/continuous_validation_v1",
        "paths.output_root",
    )
    production_output_root = str(
        config.get("paths", {}).get("production_output_root", "")
    )
    _require_equal(
        production_output_root,
        "_Data/eeg/mne_preprocessing/continuous_v1",
        "paths.production_output_root",
    )
    if production_output_root == output_root:
        raise ValueError("Validation and production output roots must be distinct.")

    production = config.get("production_surface", {})
    _require_equal(
        production.get("raw_manifest"),
        "_Data/eeg/manifest/raw_eeg_file_manifest.csv",
        "production_surface.raw_manifest",
    )
    _require_equal(
        production.get("expected_readable_edf_count"),
        95,
        "production_surface.expected_readable_edf_count",
    )
    _require_equal(
        production.get("deterministic_order"),
        "parsed_recording_id_then_source_filename",
        "production_surface.deterministic_order",
    )
    _require_equal(
        production.get("minimum_free_bytes"),
        40 * 1024**3,
        "production_surface.minimum_free_bytes",
    )


def load_config(path: Path) -> dict[str, Any]:
    """Read, validate, and return the tracked continuous-preprocessing config.

    Args:
        path: YAML file path.

    Returns:
        Parsed configuration mapping with a ``_provenance`` hash entry.

    Side effects:
        Reads one file. It writes nothing.
    """

    raw_bytes = path.read_bytes()
    parsed = yaml.safe_load(raw_bytes)
    if not isinstance(parsed, dict):
        raise ValueError("Continuous-preprocessing config must be a YAML mapping.")
    validate_config(parsed)
    config = dict(parsed)
    config["_provenance"] = {
        "path": path.as_posix(),
        "sha256": sha256_bytes(raw_bytes),
    }
    return config
