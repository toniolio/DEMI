"""Contracts and tested helpers for DEMI trial-level time-frequency power.

This module implements the accepted stage-16 signal-processing contract while
keeping orchestration, file discovery, and aggregate run publication in the
numbered driver.  It validates the exact channel/frequency/cycle/time surface,
resamples in-memory Epochs copies with MNE's anti-aliased polyphase route,
computes Morlet power, derives trial-matched ``red_on`` normalization, and
publishes reopenable NumPy arrays atomically.

Inputs:
    Already accepted ``epochs_v1`` response-onset, response-end, and ``red_on``
    Epochs objects plus their row-aligned metadata.

Outputs:
    Float32 raw/dB task arrays, float32 baseline-mean arrays, deterministic
    log-power views, diagnostic metadata, hashes, and validation evidence.

This module explicitly does not modify accepted FIFs, interpolate or reject an
epoch, run AutoReject, compute complex/phase/connectivity output, apply CSD,
average bands, define ROIs, exclude participants, or fit statistical models.
"""

from __future__ import annotations

from dataclasses import dataclass
import gc
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import uuid

import mne
from mne.time_frequency import morlet, tfr_array_morlet
import numpy as np
import pandas as pd
import yaml

from epoch_construction import (
    atomic_write_csv,
    atomic_write_json,
    atomic_write_parquet,
    atomic_write_text,
    content_fingerprint,
    ordered_key_hash,
    sha256_file,
)
from event_epoch_eligibility import stable_frame_hash


TFR_NAMESPACE = "tfr_v1"
TFR_SCHEMA_VERSION = 1
EXPECTED_RECORDING_COUNT = 81
EXPECTED_TRIAL_COUNT = 8_798
EXPECTED_STRICT_CLEAN_COUNT = 8_789
EXPECTED_DURATION_WARNING_COUNT = 9
EXPECTED_FILE49_COUNT = 117
INPUT_SFREQ = 1_000.0
TARGET_SFREQ = 100.0
RESAMPLE_METHOD = "polyphase"
RESAMPLE_WINDOW = ("kaiser", 5.0)
RESAMPLE_PAD = "reflect"
TASK_TMIN = -0.5
TASK_TMAX = 1.5
BASELINE_TMIN = -0.5
BASELINE_TMAX = -0.2
EXPECTED_TASK_SAMPLE_COUNT = 201
TFR_DECIM = 1
CALCULATION_DTYPE = np.dtype("float64")
PERSISTED_DTYPE = np.dtype("float32")
TASK_FAMILIES = ("response_onset", "response_end")
ALL_FAMILIES = (*TASK_FAMILIES, "red_on")

PRIMARY_SCALP_CHANNELS: tuple[str, ...] = (
    "FP1", "FP2", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ",
    "FC4", "FT8", "T7", "C3", "CZ", "C4", "T8", "TP7", "CP3", "CPZ",
    "CP4", "TP8", "P7", "P3", "PZ", "P4", "P8", "O1", "OZ", "O2",
)
EOG_CHANNELS = ("HEO", "VEO")
EMG_CHANNELS = ("EMG-L", "EMG-A")
EXCLUDED_PRIMARY_CHANNELS = ("M1", "M2", *EOG_CHANNELS, *EMG_CHANNELS, "Trigger")

FREQUENCIES_HZ = np.arange(4.0, 41.0, 1.0, dtype=np.float64)
N_CYCLES = 4.0 + (FREQUENCIES_HZ - 4.0) / 6.0
EXPECTED_TASK_TIMES = np.arange(-50, 151, dtype=np.float64) / TARGET_SFREQ

CROSS_FAMILY_IDENTITY_COLUMNS: tuple[str, ...] = (
    "canonical_order_index",
    "canonical_event_key",
    "ledger_version",
    "event_policy_version",
    "behavioural_id",
    "eeg_source_id",
    "source_recording_filename",
    "source_file_role",
    "task_filename",
    "offset_session",
    "offset_block",
    "offset_trial",
    "audit_trial_count",
    "raw_trial_sequence",
    "join_trial_count",
    "offset_trial_count",
    "physical",
    "condition_semantics",
    "event_source_type",
    "selected_event_source",
    "trial_event_provenance_json",
    "primary_eligibility",
    "primary_eligibility_reason_code",
    "strict_clean_eligibility",
    "strict_clean_eligibility_reason_code",
    "duration_warning_flag",
    "event_policy_reason_code",
    "accepted_policy_ids",
    "accepted_specific_policy_id",
    "accepted_policy_action",
    "continuous_qc_warning",
    "interpolation_count",
    "interpolation_proportion",
    "interpolation_denominator",
    "continuous_v2_run_id",
    "continuous_v2_manifest_path",
    "continuous_v2_manifest_id",
    "post_ica_derivative_kind",
    "post_ica_derivative_path",
    "continuous_source_sha256",
    "ica_terminal_route",
    "post_ica_availability_reason_code",
    "future_readiness_reason_code",
)


@dataclass(frozen=True)
class ArrayContract:
    """Expected shape, dtype, axes, and numerical constraints for one array."""

    shape: tuple[int, ...]
    dtype: np.dtype
    axis_order: tuple[str, ...]
    finite: bool = True
    positive: bool = False


def load_and_validate_config(path: Path) -> tuple[dict[str, Any], str]:
    """Load the tracked YAML configuration and enforce accepted invariants.

    Args:
        path: Tracked ``tfr_config_v1.yaml`` path.

    Returns:
        Parsed configuration and SHA-256 file digest.

    Raises:
        RuntimeError: If any accepted scientific or safety value differs.

    Side effects:
        Reads and hashes the configuration.
    """

    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    expected = {
        "pipeline_version": TFR_NAMESPACE,
        "input_sfreq": INPUT_SFREQ,
        "target_sfreq": TARGET_SFREQ,
        "frequencies": FREQUENCIES_HZ.tolist(),
        "cycles": N_CYCLES.tolist(),
        "channels": list(PRIMARY_SCALP_CHANNELS),
        "task_times": [TASK_TMIN, TASK_TMAX, EXPECTED_TASK_SAMPLE_COUNT],
        "baseline": [BASELINE_TMIN, BASELINE_TMAX],
    }
    observed_frequencies = np.arange(
        float(config["tfr"]["frequencies_hz"]["start"]),
        float(config["tfr"]["frequencies_hz"]["stop_inclusive"]) + 0.5,
        float(config["tfr"]["frequencies_hz"]["step"]),
    )
    observed = {
        "pipeline_version": config.get("pipeline_version"),
        "input_sfreq": float(config["accepted_surface"]["input_sampling_frequency_hz"]),
        "target_sfreq": float(config["sampling"]["target_frequency_hz"]),
        "frequencies": observed_frequencies.tolist(),
        "cycles": (4.0 + (observed_frequencies - 4.0) / 6.0).tolist(),
        "channels": list(config["channels"]["primary_scalp"]),
        "task_times": [
            float(config["time"]["task_tmin_seconds"]),
            float(config["time"]["task_tmax_seconds"]),
            int(config["time"]["task_sample_count"]),
        ],
        "baseline": [
            float(config["time"]["baseline_tmin_seconds"]),
            float(config["time"]["baseline_tmax_seconds"]),
        ],
    }
    if observed != expected:
        raise RuntimeError(f"TFR configuration differs from accepted contract: {observed}")
    required_exact = {
        ("sampling", "method"): "polyphase",
        ("sampling", "pad"): "reflect",
        ("tfr", "method"): "morlet",
        ("tfr", "output"): "power",
        ("tfr", "decim"): 1,
        ("tfr", "zero_mean"): True,
        ("tfr", "use_fft"): True,
        ("precision", "calculation_dtype"): "float64",
        ("precision", "persisted_power_dtype"): "float32",
        ("precision", "persisted_db_dtype"): "float32",
        ("precision", "persisted_baseline_dtype"): "float32",
    }
    for keys, expected_value in required_exact.items():
        value: Any = config
        for key in keys:
            value = value[key]
        if value != expected_value:
            raise RuntimeError(f"configuration {'.'.join(keys)} differs: {value}")
    if list(config["sampling"]["window"]) != ["kaiser", 5.0]:
        raise RuntimeError("resampling window must be [kaiser, 5.0]")
    forbidden_true = (
        "mutate_accepted_epochs",
        "epoch_interpolation",
        "epoch_rejection",
        "fixed_amplitude_rejection",
        "run_autoreject",
        "compute_csd",
        "compute_complex_output",
        "compute_phase_output",
        "compute_connectivity",
        "compute_band_power",
        "define_roi",
        "participant_inclusion_decisions",
    )
    if any(bool(config["safety"][key]) for key in forbidden_true):
        raise RuntimeError("TFR safety configuration enables an unauthorized operation")
    return config, sha256_file(path)


def expected_wavelet_support() -> dict[str, Any]:
    """Return exact installed-MNE kernel support for the accepted schedule."""

    wavelets = morlet(
        sfreq=TARGET_SFREQ,
        freqs=FREQUENCIES_HZ,
        n_cycles=N_CYCLES,
        zero_mean=True,
    )
    samples = np.asarray([len(item) for item in wavelets], dtype=int)
    half_support = (samples - 1) / (2.0 * TARGET_SFREQ)
    maximum = float(half_support.max())
    if maximum > 0.8:
        raise RuntimeError(f"wavelet support exceeds accepted edge allowance: {maximum}")
    if TASK_TMIN < -1.5 + maximum or TASK_TMAX > 2.5 - maximum:
        raise RuntimeError("accepted task retained interval is not edge-safe")
    if BASELINE_TMIN < -1.5 + maximum or BASELINE_TMAX > 0.8 - maximum:
        raise RuntimeError("accepted red_on baseline is not edge-safe")
    return {
        "wavelet_samples": samples.tolist(),
        "half_support_seconds": half_support.tolist(),
        "maximum_half_support_seconds": maximum,
        "task_edge_safe_interval_seconds": [-1.5 + maximum, 2.5 - maximum],
        "red_on_edge_safe_interval_seconds": [-1.5 + maximum, 0.8 - maximum],
    }


def validate_channel_surface(epochs: mne.BaseEpochs) -> None:
    """Require the complete accepted input and exact primary scalp types/order."""

    missing = set((*PRIMARY_SCALP_CHANNELS, *EXCLUDED_PRIMARY_CHANNELS)) - set(
        epochs.ch_names
    )
    if missing:
        raise RuntimeError(f"accepted epoch shard lacks required channels: {sorted(missing)}")
    types = dict(zip(epochs.ch_names, epochs.get_channel_types(), strict=True))
    if any(types[channel] != "eeg" for channel in PRIMARY_SCALP_CHANNELS):
        raise RuntimeError("one or more primary scalp channels are not typed EEG")
    expected_input_order = [
        channel for channel in epochs.ch_names if channel in PRIMARY_SCALP_CHANNELS
    ]
    if expected_input_order != list(PRIMARY_SCALP_CHANNELS):
        raise RuntimeError("accepted primary scalp channel order differs")


def validate_family_identity(metadata: Mapping[str, pd.DataFrame]) -> str:
    """Require one-to-one, row-aligned identity and provenance across families."""

    if set(metadata) != set(ALL_FAMILIES):
        raise RuntimeError("all three accepted family metadata frames are required")
    reference = metadata["response_onset"].reset_index(drop=True)
    if reference["canonical_event_key"].duplicated().any():
        raise RuntimeError("baseline_match_duplicate_canonical_event_key")
    missing_columns = set(CROSS_FAMILY_IDENTITY_COLUMNS) - set(reference.columns)
    if missing_columns:
        raise RuntimeError(f"accepted metadata lacks identity columns: {sorted(missing_columns)}")
    reference_identity = reference.loc[:, CROSS_FAMILY_IDENTITY_COLUMNS]
    for family in ALL_FAMILIES:
        frame = metadata[family].reset_index(drop=True)
        if len(frame) != len(reference):
            raise RuntimeError(f"baseline_match_missing_rows:{family}")
        if frame["canonical_event_key"].duplicated().any():
            raise RuntimeError(f"baseline_match_duplicate_canonical_event_key:{family}")
        try:
            pd.testing.assert_frame_equal(
                frame.loc[:, CROSS_FAMILY_IDENTITY_COLUMNS],
                reference_identity,
                check_dtype=True,
                check_exact=True,
            )
        except AssertionError as error:
            raise RuntimeError(f"baseline_match_identity_or_provenance_mismatch:{family}") from error
    return ordered_key_hash(reference["canonical_event_key"].astype(str))


def resample_for_tfr(epochs: mne.BaseEpochs) -> mne.Epochs:
    """Copy, select the 30 scalp channels, and anti-alias resample to 100 Hz."""

    if not epochs.preload:
        raise RuntimeError("accepted Epochs must be preloaded before in-memory resampling")
    if not np.isclose(epochs.info["sfreq"], INPUT_SFREQ, atol=1e-9, rtol=0.0):
        raise RuntimeError(f"incorrect input sampling rate: {epochs.info['sfreq']}")
    validate_channel_surface(epochs)
    # The accepted input order has already been checked above. MNE 1.12's
    # BaseEpochs.pick preserves that input order and does not expose an
    # ``ordered`` keyword.
    work = epochs.copy().pick(list(PRIMARY_SCALP_CHANNELS))
    work.resample(
        TARGET_SFREQ,
        method=RESAMPLE_METHOD,
        window=RESAMPLE_WINDOW,
        pad=RESAMPLE_PAD,
        n_jobs=1,
        verbose="ERROR",
    )
    if not np.isclose(work.info["sfreq"], TARGET_SFREQ, atol=1e-12, rtol=0.0):
        raise RuntimeError(f"incorrect resampled rate: {work.info['sfreq']}")
    if work.ch_names != list(PRIMARY_SCALP_CHANNELS):
        raise RuntimeError("resampled primary channel order differs")
    # MNE resamples the finite array length, so 4,001 native samples become
    # 400 samples at one tenth the rate (ending at +2.49 s), and 2,301 become
    # 230 (ending at +0.79 s). The required -0.5/+1.5 and -0.5/-0.2 endpoints
    # remain exactly represented and are checked separately below.
    expected_count = int(round(len(epochs.times) * TARGET_SFREQ / INPUT_SFREQ))
    expected_times = epochs.tmin + np.arange(expected_count) / TARGET_SFREQ
    if len(work.times) != expected_count or not np.allclose(
        work.times, expected_times, atol=1e-12, rtol=0.0
    ):
        raise RuntimeError("resampled time vector differs from explicit 100-Hz grid")
    return work


def compute_morlet_power(epochs: mne.BaseEpochs) -> tuple[np.ndarray, np.ndarray]:
    """Compute float64 trial-level Morlet power with no TFR decimation."""

    data = epochs.get_data(copy=False)
    if data.dtype != CALCULATION_DTYPE:
        data = np.asarray(data, dtype=CALCULATION_DTYPE)
    if not np.isfinite(data).all():
        raise RuntimeError("accepted_signal_nonfinite")
    power = tfr_array_morlet(
        data,
        sfreq=TARGET_SFREQ,
        freqs=FREQUENCIES_HZ,
        n_cycles=N_CYCLES,
        zero_mean=True,
        use_fft=True,
        decim=TFR_DECIM,
        output="power",
        n_jobs=1,
        verbose="ERROR",
    )
    if power.dtype != CALCULATION_DTYPE:
        power = np.asarray(power, dtype=CALCULATION_DTYPE)
    expected_shape = (len(epochs), len(PRIMARY_SCALP_CHANNELS), len(FREQUENCIES_HZ), len(epochs.times))
    if power.shape != expected_shape:
        raise RuntimeError(f"unexpected Morlet output shape: {power.shape} != {expected_shape}")
    if not np.isfinite(power).all():
        raise RuntimeError("raw_power_nonfinite")
    if np.any(power <= 0):
        raise RuntimeError("raw_power_nonpositive_log_view_unavailable")
    return power, epochs.times.copy()


def inclusive_time_indices(times: np.ndarray, tmin: float, tmax: float) -> np.ndarray:
    """Return inclusive exact-grid indices and fail on missing endpoints."""

    indices = np.flatnonzero((times >= tmin - 1e-12) & (times <= tmax + 1e-12))
    if not len(indices):
        raise RuntimeError(f"time interval absent: {tmin} to {tmax}")
    if not np.isclose(times[indices[0]], tmin, atol=1e-12, rtol=0.0):
        raise RuntimeError(f"time interval start absent: {tmin}")
    if not np.isclose(times[indices[-1]], tmax, atol=1e-12, rtol=0.0):
        raise RuntimeError(f"time interval end absent: {tmax}")
    return indices


def derive_baseline_mean(red_power: np.ndarray, red_times: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    """Average red_on power over time first, separately by trial/channel/frequency."""

    indices = inclusive_time_indices(red_times, BASELINE_TMIN, BASELINE_TMAX)
    baseline = red_power[..., indices].mean(axis=-1, dtype=np.float64)
    if baseline.dtype != CALCULATION_DTYPE:
        raise RuntimeError("baseline calculation did not remain float64")
    invalid_nonfinite = ~np.isfinite(baseline)
    invalid_nonpositive = baseline <= 0
    if invalid_nonfinite.any():
        raise RuntimeError(f"baseline_nonfinite:{int(invalid_nonfinite.sum())}")
    if invalid_nonpositive.any():
        raise RuntimeError(f"baseline_nonpositive:{int(invalid_nonpositive.sum())}")
    return baseline, {
        "sample_count": int(len(indices)),
        "first_time_seconds": float(red_times[indices[0]]),
        "last_time_seconds": float(red_times[indices[-1]]),
        "nonfinite_count": 0,
        "nonpositive_count": 0,
        "minimum": float(baseline.min()),
        "maximum": float(baseline.max()),
    }


def crop_task_power(power: np.ndarray, times: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Copy the exact inclusive -0.5 to +1.5-s task power surface."""

    indices = inclusive_time_indices(times, TASK_TMIN, TASK_TMAX)
    retained_times = times[indices].copy()
    if len(indices) != EXPECTED_TASK_SAMPLE_COUNT or not np.allclose(
        retained_times, EXPECTED_TASK_TIMES, atol=1e-12, rtol=0.0
    ):
        raise RuntimeError("retained task time vector differs from accepted 201-sample grid")
    retained = np.asarray(power[..., indices], dtype=np.float64)
    return retained, retained_times


def normalize_db(task_power: np.ndarray, baseline_mean: np.ndarray) -> np.ndarray:
    """Apply the exact float64 trial-matched power-ratio dB transform."""

    if task_power.dtype != CALCULATION_DTYPE or baseline_mean.dtype != CALCULATION_DTYPE:
        raise RuntimeError("dB inputs must remain float64")
    if task_power.shape[:3] != baseline_mean.shape:
        raise RuntimeError(
            f"task/baseline shape mismatch: {task_power.shape[:3]} != {baseline_mean.shape}"
        )
    if not np.isfinite(baseline_mean).all():
        raise RuntimeError("baseline_nonfinite")
    if np.any(baseline_mean <= 0):
        raise RuntimeError("baseline_nonpositive")
    db = 10.0 * np.log10(task_power / baseline_mean[..., np.newaxis])
    if db.dtype != CALCULATION_DTYPE or not np.isfinite(db).all():
        raise RuntimeError("db_nonfinite")
    return db


def derive_log_power(raw_power: np.ndarray) -> np.ndarray:
    """Return deterministic unnormalized ``10 * log10(raw_power)`` in float64."""

    values = np.asarray(raw_power, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("raw power contains non-finite values")
    if np.any(values <= 0):
        raise ValueError("raw power must be strictly positive for log power")
    result = 10.0 * np.log10(values)
    if not np.isfinite(result).all():
        raise ValueError("derived log power is non-finite")
    return result


def robust_upper_z(values: np.ndarray) -> np.ndarray:
    """Return accepted recording-relative robust z scores on log10 p2p."""

    logged = np.log10(np.maximum(values, np.finfo(float).tiny))
    median = np.median(logged, axis=0, keepdims=True)
    mad = np.median(np.abs(logged - median), axis=0, keepdims=True)
    scale = np.where(mad > 0, mad / 0.6744897501960817, np.inf)
    return (logged - median) / scale


def artifact_diagnostics(epochs: mne.BaseEpochs, family: str) -> pd.DataFrame:
    """Calculate non-decisional trial diagnostics without changing eligibility."""

    if family not in ALL_FAMILIES or epochs.metadata is None:
        raise RuntimeError("artifact diagnostics require accepted family metadata")
    validate_channel_surface(epochs)
    scalp = epochs.get_data(picks=list(PRIMARY_SCALP_CHANNELS), copy=True) * 1e6
    eog = epochs.get_data(picks=list(EOG_CHANNELS), copy=True) * 1e6
    emg = epochs.get_data(picks=list(EMG_CHANNELS), copy=True) * 1e6
    if not (np.isfinite(scalp).all() and np.isfinite(eog).all() and np.isfinite(emg).all()):
        raise RuntimeError("artifact_diagnostic_input_nonfinite")
    scalp_p2p = np.ptp(scalp, axis=-1)
    scalp_jump = np.max(np.abs(np.diff(scalp, axis=-1)), axis=-1)
    robust_z = robust_upper_z(scalp_p2p)
    metadata = epochs.metadata.reset_index(drop=True)
    frame = metadata[
        [
            "canonical_order_index",
            "canonical_event_key",
            "behavioural_id",
            "eeg_source_id",
            "source_recording_filename",
            "condition_semantics",
            "physical",
            "strict_clean_eligibility",
            "duration_warning_flag",
            "continuous_qc_warning",
        ]
    ].copy()
    frame.insert(0, "family", family)
    frame["scalp_p2p_max_uv"] = scalp_p2p.max(axis=1)
    frame["scalp_p2p_median_uv"] = np.median(scalp_p2p, axis=1)
    frame["n_scalp_channels_robust_logp2p_z_gt_6"] = (robust_z > 6.0).sum(axis=1)
    frame["diagnostic_robust_logp2p_z_gt_6"] = (
        frame["n_scalp_channels_robust_logp2p_z_gt_6"] > 0
    )
    frame["scalp_jump_max_uv_per_native_sample"] = scalp_jump.max(axis=1)
    frame["diagnostic_jump_gt_50uv_per_native_sample"] = np.any(
        scalp_jump > 50.0, axis=1
    )
    frame["scalp_p2p_min_uv"] = scalp_p2p.min(axis=1)
    frame["diagnostic_flat_p2p_lt_1uv"] = np.any(scalp_p2p < 1.0, axis=1)
    frame["eog_p2p_max_uv"] = np.ptp(eog, axis=-1).max(axis=1)
    frame["emg_p2p_max_uv"] = np.ptp(emg, axis=-1).max(axis=1)
    frame["diagnostic_only"] = True
    frame["changes_primary_eligibility"] = False
    del scalp, eog, emg, scalp_p2p, scalp_jump, robust_z
    gc.collect()
    return frame


def atomic_write_npy(path: Path, array: np.ndarray) -> dict[str, Any]:
    """Persist one C-contiguous NumPy array through an atomic sibling rename."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}.npy")
    with temporary.open("xb") as handle:
        np.save(handle, np.ascontiguousarray(array), allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    return {
        "format": "numpy_npy",
        "path": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "c_contiguous": bool(array.flags.c_contiguous),
    }


def _chunked_numeric_scan(array: np.ndarray, *, positive: bool) -> dict[str, Any]:
    """Scan a memory-mapped array by trial without allocating a full mask."""

    nonfinite = 0
    nonpositive = 0
    minimum = math.inf
    maximum = -math.inf
    for index in range(array.shape[0]):
        chunk = np.asarray(array[index])
        nonfinite += int((~np.isfinite(chunk)).sum())
        if positive:
            nonpositive += int((chunk <= 0).sum())
        minimum = min(minimum, float(np.min(chunk)))
        maximum = max(maximum, float(np.max(chunk)))
    return {
        "nonfinite_count": nonfinite,
        "nonpositive_count": nonpositive,
        "minimum": minimum,
        "maximum": maximum,
    }


def reopen_validate_npy(
    path: Path,
    contract: ArrayContract,
    *,
    expected_sha256: str | None = None,
) -> dict[str, Any]:
    """Reopen, hash, and bounded-scan one persisted NumPy array."""

    array = np.load(path, mmap_mode="r", allow_pickle=False)
    if array.shape != contract.shape:
        raise RuntimeError(f"array shape mismatch for {path}: {array.shape}")
    if array.dtype != contract.dtype:
        raise RuntimeError(f"array dtype mismatch for {path}: {array.dtype}")
    if not array.flags.c_contiguous:
        raise RuntimeError(f"array is not C-contiguous: {path}")
    observed_sha = sha256_file(path)
    if expected_sha256 is not None and observed_sha != expected_sha256:
        raise RuntimeError(f"array hash mismatch: {path}")
    scan = _chunked_numeric_scan(array, positive=contract.positive)
    if contract.finite and scan["nonfinite_count"]:
        raise RuntimeError(f"array contains non-finite values: {path}")
    if contract.positive and scan["nonpositive_count"]:
        raise RuntimeError(f"array contains non-positive values: {path}")
    return {
        "valid": True,
        "path": path.as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": observed_sha,
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "axis_order": list(contract.axis_order),
        **scan,
    }


def validate_stored_formula(
    raw_path: Path,
    db_path: Path,
    baseline_path: Path,
) -> dict[str, Any]:
    """Reproduce dB and log-power formulas on deterministic stored samples."""

    raw = np.load(raw_path, mmap_mode="r", allow_pickle=False)
    db = np.load(db_path, mmap_mode="r", allow_pickle=False)
    baseline = np.load(baseline_path, mmap_mode="r", allow_pickle=False)
    trial_indices = sorted({0, raw.shape[0] // 2, raw.shape[0] - 1})
    channel_indices = sorted({0, raw.shape[1] // 2, raw.shape[1] - 1})
    frequency_indices = sorted({0, raw.shape[2] // 3, (2 * raw.shape[2]) // 3, raw.shape[2] - 1})
    time_indices = sorted({0, raw.shape[3] // 2, raw.shape[3] - 1})
    db_differences = []
    log_differences = []
    for i in trial_indices:
        for c in channel_indices:
            for f in frequency_indices:
                for t in time_indices:
                    raw_value = float(raw[i, c, f, t])
                    expected_db = 10.0 * math.log10(raw_value / float(baseline[i, c, f]))
                    db_differences.append(abs(float(db[i, c, f, t]) - expected_db))
                    expected_log = 10.0 * math.log10(raw_value)
                    derived_log = float(derive_log_power(np.asarray([raw_value]))[0])
                    log_differences.append(abs(expected_log - derived_log))
    maximum_db_difference = max(db_differences)
    if maximum_db_difference > 2e-4:
        raise RuntimeError(f"stored dB formula differs by {maximum_db_difference}")
    if max(log_differences) > 1e-12:
        raise RuntimeError("deterministic log-power derivation differs")
    return {
        "sample_count": len(db_differences),
        "maximum_absolute_db_difference": maximum_db_difference,
        "maximum_absolute_log_derivation_difference": max(log_differences),
        "db_formula": "10 * log10(raw_power / baseline_mean_power)",
        "log_power_formula": "10 * log10(raw_power)",
    }


def source_snapshot(paths: Sequence[Path]) -> list[dict[str, Any]]:
    """Capture accepted source path, size, and nanosecond modification time."""

    return [
        {
            "path": path.resolve().as_posix(),
            "size_bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for path in sorted((item.resolve() for item in paths), key=str)
    ]


def compare_source_snapshots(
    before: Sequence[Mapping[str, Any]], after: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """Fail closed if an accepted epoch source size or mtime changed."""

    if list(before) != list(after):
        raise RuntimeError("accepted_epoch_source_size_or_mtime_changed")
    return {
        "unchanged": True,
        "file_count": len(before),
        "total_bytes": sum(int(row["size_bytes"]) for row in before),
    }


def recording_fingerprint(value: Mapping[str, Any]) -> str:
    """Return a deterministic complete cache fingerprint for one recording."""

    return content_fingerprint(dict(value))


def assess_cached_recording(
    recording_dir: Path,
    expected_fingerprint: str,
    expected_contracts: Mapping[str, ArrayContract],
) -> dict[str, Any]:
    """Validate whether a complete recording result can be safely reused."""

    manifest_path = recording_dir / "manifest.json"
    if not recording_dir.exists():
        return {"action": "build", "reason": "no_previous_result"}
    if not recording_dir.is_dir() or not manifest_path.is_file():
        return {"action": "rebuild", "reason": "incomplete_previous_result"}
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"action": "rebuild", "reason": "unreadable_previous_manifest"}
    if manifest.get("status") != "complete":
        return {"action": "rebuild", "reason": "previous_status_not_complete"}
    if manifest.get("fingerprint") != expected_fingerprint:
        return {"action": "rebuild", "reason": "provenance_fingerprint_changed"}
    arrays = manifest.get("arrays", {})
    try:
        for label, contract in expected_contracts.items():
            descriptor = arrays[label]
            path = recording_dir / descriptor["relative_path"]
            if path.stat().st_size != int(descriptor["size_bytes"]):
                raise RuntimeError("size")
            reopen_validate_npy(path, contract, expected_sha256=descriptor["sha256"])
        for task_family in TASK_FAMILIES:
            validate_stored_formula(
                recording_dir / arrays[f"{task_family}_raw"]["relative_path"],
                recording_dir / arrays[f"{task_family}_db"]["relative_path"],
                recording_dir / arrays["red_on_baseline_mean"]["relative_path"],
            )
    except (KeyError, OSError, RuntimeError, ValueError):
        return {"action": "rebuild", "reason": "cached_artifact_validation_failed"}
    return {"action": "reuse", "reason": "complete_hash_reopen_valid", "manifest": manifest}


def require_ignored_output_root(repo_root: Path, output_root: Path) -> None:
    """Require exactly the ignored versioned `_Data/eeg/tfr_v1` namespace."""

    expected = (repo_root / "_Data/eeg/tfr_v1").resolve()
    if output_root.resolve() != expected:
        raise RuntimeError(f"TFR output root must be exactly {expected}")
    result = os.spawnlp(
        os.P_WAIT,
        "git",
        "git",
        "-C",
        str(repo_root),
        "check-ignore",
        "--quiet",
        "_Data/eeg/tfr_v1",
    )
    if result != 0:
        raise RuntimeError("_Data/eeg/tfr_v1 is not ignored by git")


def forbidden_output_scan(output_root: Path) -> list[str]:
    """Return output paths that imply an unauthorized analysis expansion."""

    forbidden_terms = (
        "autoreject",
        "interpolated_epoch",
        "rejected_epoch",
        "complex_coeff",
        "phase_lock",
        "connectivity",
        "csd",
        "band_power",
        "roi",
        "model",
        "id86",
        "54_1",
    )
    matches = []
    for path in output_root.rglob("*"):
        if path.is_file() and any(term in path.name.lower() for term in forbidden_terms):
            matches.append(path.relative_to(output_root).as_posix())
    return sorted(matches)


__all__ = [
    "ALL_FAMILIES",
    "ArrayContract",
    "BASELINE_TMAX",
    "BASELINE_TMIN",
    "CALCULATION_DTYPE",
    "CROSS_FAMILY_IDENTITY_COLUMNS",
    "EXPECTED_DURATION_WARNING_COUNT",
    "EXPECTED_FILE49_COUNT",
    "EXPECTED_RECORDING_COUNT",
    "EXPECTED_STRICT_CLEAN_COUNT",
    "EXPECTED_TASK_SAMPLE_COUNT",
    "EXPECTED_TASK_TIMES",
    "EXPECTED_TRIAL_COUNT",
    "EXCLUDED_PRIMARY_CHANNELS",
    "FREQUENCIES_HZ",
    "INPUT_SFREQ",
    "N_CYCLES",
    "PERSISTED_DTYPE",
    "PRIMARY_SCALP_CHANNELS",
    "RESAMPLE_METHOD",
    "RESAMPLE_PAD",
    "RESAMPLE_WINDOW",
    "TARGET_SFREQ",
    "TASK_FAMILIES",
    "TASK_TMAX",
    "TASK_TMIN",
    "TFR_DECIM",
    "TFR_NAMESPACE",
    "TFR_SCHEMA_VERSION",
    "artifact_diagnostics",
    "assess_cached_recording",
    "atomic_write_csv",
    "atomic_write_json",
    "atomic_write_npy",
    "atomic_write_parquet",
    "atomic_write_text",
    "compare_source_snapshots",
    "compute_morlet_power",
    "content_fingerprint",
    "crop_task_power",
    "derive_baseline_mean",
    "derive_log_power",
    "expected_wavelet_support",
    "forbidden_output_scan",
    "inclusive_time_indices",
    "load_and_validate_config",
    "normalize_db",
    "ordered_key_hash",
    "recording_fingerprint",
    "reopen_validate_npy",
    "require_ignored_output_root",
    "sha256_file",
    "source_snapshot",
    "stable_frame_hash",
    "validate_family_identity",
    "validate_stored_formula",
]
