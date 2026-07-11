"""Run deterministic, unit-aware, per-channel QC on raw DEMI EDF files.

This script turns the earlier sampled, channel-type-level raw-QC surface into
one reproducible row per EDF channel. It exists to provide evidence for later
manual bad-channel, interpolation, line-noise, and ICA planning without making
any of those decisions.

Inputs:

* raw EDF files below ``_Data/eeg/raw``;
* original per-signal EDF header metadata read directly by ``channel_qc``;
* MNE metadata and read-only values returned by ``Raw.get_data()``.

Local outputs below ``_Data/eeg/mne_preprocessing/channel_qc_v1``:

* ``channel_qc_metrics.csv``: one row per file/channel, including original EDF
  unit/range/calibration, MNE unit metadata, full-duration descriptive metrics,
  scale-free spectral ratios, coverage, and descriptive ranks;
* ``channel_qc_file_summary.csv``: compact per-file counts;
* ``channel_qc_candidate_review.csv``: a small descriptive review queue based
  on hard exceptions and extreme ranks, never final bad-channel labels;
* ``stim_transition_summary.csv``: Trigger transitions/digital states kept
  separate from continuous physical-amplitude interpretation;
* ``channel_qc_summary.md`` and ``channel_qc_run_manifest.json``.

Coverage and memory behavior:

* each EDF is opened and closed before the next file;
* every signal sample is scanned once in deterministic contiguous chunks;
* exact full-duration counts/ranges/moments, repeated samples, run lengths,
  rails, transitions, and maximum steps are computed;
* robust quantiles and step distributions use a deterministic, duration-wide
  decimated sample whose stride is recorded;
* Welch summaries cover all complete analysis chunks and record their coverage.

Unit boundary:

* scalp EEG declares ``uV`` in the EDF and is returned by MNE in volts;
* HEO, VEO, EMG-L, and EMG-A have unknown EDF physical dimensions, so their
  absolute amplitudes/PSD are labelled unknown and are not ranked across files;
* their morphology, frequency ratios, flatness, rail, and within-file metrics
  remain useful descriptive evidence;
* Trigger is summarized only as digital states/transitions, not continuous
  amplitude or PSD.

Safety boundaries:

* raw EDFs are read-only and never rewritten;
* no filter or notch filter is applied;
* no reference is changed;
* no bad channel is marked;
* no channel is interpolated;
* no ICA or CSD is run;
* no cleaned FIF/EDF or epochs are written;
* no participant, file, channel, or trial inclusion decision is made.

Run from the repository root:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/07_run_unit_aware_channel_qc.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.signal import welch

# MNE can initialize matplotlib even though this script creates no figures.
# Keep any cache writes out of the repository and user configuration.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "demi_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "demi_cache"))

import mne  # noqa: E402

from channel_qc import (  # noqa: E402
    EdfSignalHeader,
    classify_unit_status,
    deterministic_chunk_bounds,
    deterministic_sample_stride,
    longest_true_run,
    mne_scaling_for_unit_status,
    read_edf_signal_headers,
)


SCHEMA_VERSION = "channel_qc_v1"
RAW_EEG_DIR = Path("_Data") / "eeg" / "raw"
OUTPUT_DIR = Path("_Data") / "eeg" / "mne_preprocessing" / SCHEMA_VERSION

CHANNEL_METRICS_FILENAME = "channel_qc_metrics.csv"
FILE_SUMMARY_FILENAME = "channel_qc_file_summary.csv"
CANDIDATE_FILENAME = "channel_qc_candidate_review.csv"
STIM_SUMMARY_FILENAME = "stim_transition_summary.csv"
SUMMARY_FILENAME = "channel_qc_summary.md"
RUN_MANIFEST_FILENAME = "channel_qc_run_manifest.json"

FULL_SCAN_CHUNK_SECONDS = 120.0
MAXIMUM_ROBUST_SAMPLES_PER_CHANNEL = 200_000
WELCH_NPERSEG_SECONDS = 4.0
WELCH_OVERLAP_PROPORTION = 0.5
CANDIDATE_EXTREME_PERCENTILE = 0.995

EXPECTED_EEG_CHANNELS = (
    "FP1", "FP2", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ",
    "FC4", "FT8", "T7", "C3", "CZ", "C4", "T8", "M1", "TP7", "CP3",
    "CPZ", "CP4", "TP8", "M2", "P7", "P3", "PZ", "P4", "P8", "O1",
    "OZ", "O2",
)
EXPECTED_EOG_CHANNELS = ("HEO", "VEO")
EXPECTED_EMG_CHANNELS = ("EMG-L", "EMG-A")
EXPECTED_STIM_CHANNELS = ("TRIGGER",)
EXPECTED_CHANNEL_TYPE = {
    **{name: "eeg" for name in EXPECTED_EEG_CHANNELS},
    **{name: "eog" for name in EXPECTED_EOG_CHANNELS},
    **{name: "emg" for name in EXPECTED_EMG_CHANNELS},
    **{name: "stim" for name in EXPECTED_STIM_CHANNELS},
}
EXPECTED_CHANNEL_ORDER = tuple(EXPECTED_CHANNEL_TYPE)

EDF_FILENAME_RE = re.compile(
    r"^demi_(?P<eeg_recording_id>\d{1,3})(?:_(?P<split_part>\d+))?"
    r"(?P<label>.*)\.edf$",
    flags=re.IGNORECASE,
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script path."""

    return Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    """Parse optional local-development limits.

    ``--max-files`` is useful only for a focused smoke run. The requested
    evidence run uses the default and therefore processes every EDF.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Process only the first N sorted EDFs for a local smoke test.",
    )
    return parser.parse_args()


def normalize_channel_name(value: object) -> str:
    """Normalize channel labels by stripping whitespace and uppercasing."""

    return str(value or "").strip().upper()


def parse_edf_filename(path: Path) -> dict[str, Any]:
    """Parse recording identity and file role from a DEMI EDF filename."""

    match = EDF_FILENAME_RE.match(path.name)
    if match is None:
        raise ValueError(f"unexpected DEMI EDF filename: {path.name}")
    recording_id = int(match.group("eeg_recording_id"))
    split_part_text = match.group("split_part")
    split_part = int(split_part_text) if split_part_text else None
    label = match.group("label").lower()
    if split_part is not None:
        role = "split_part"
    elif "concatenated" in label:
        role = "concatenated"
    else:
        role = "single"
    return {
        "eeg_recording_id": recording_id,
        "eeg_recording_id_padded": f"{recording_id:03d}",
        "source_filename": path.name,
        "file_role": role,
        "split_part": split_part,
    }


def discover_edf_files(raw_dir: Path, maximum: int | None = None) -> list[Path]:
    """Return sorted EDF paths, optionally limited for a smoke run."""

    if not raw_dir.is_dir():
        raise RuntimeError(f"missing raw EDF directory: {raw_dir}")
    paths = sorted(
        path for path in raw_dir.iterdir()
        if path.is_file() and path.suffix.lower() == ".edf"
    )
    if not paths:
        raise RuntimeError(f"no EDF files found below {raw_dir}")
    if maximum is not None:
        if maximum <= 0:
            raise ValueError("--max-files must be positive")
        paths = paths[:maximum]
    return paths


def expected_type_for_channel(channel_name: str) -> str:
    """Return the fail-closed configured type for a DEMI channel label."""

    normalized = normalize_channel_name(channel_name)
    if normalized not in EXPECTED_CHANNEL_TYPE:
        raise RuntimeError(f"unmapped raw channel label: {channel_name}")
    return EXPECTED_CHANNEL_TYPE[normalized]


def apply_channel_types(raw: mne.io.BaseRaw) -> None:
    """Apply reviewed in-memory auxiliary/stim types without touching disk."""

    observed = tuple(normalize_channel_name(name) for name in raw.ch_names)
    if observed != EXPECTED_CHANNEL_ORDER:
        raise RuntimeError(
            "raw channel order differs from the reviewed 37-channel DEMI contract"
        )
    # Only change types that are wrong in the EDF/MNE header. Re-applying an
    # already-correct stim type can replace FIFF_UNIT_NONE with a voltage unit
    # in some MNE versions, which would obscure the Trigger unit contract.
    current_types = raw.get_channel_types()
    mapping = {
        name: expected_type_for_channel(name)
        for name, current_type in zip(raw.ch_names, current_types)
        if current_type != expected_type_for_channel(name)
    }
    if not mapping:
        return
    try:
        raw.set_channel_types(mapping, on_unit_change="ignore")
    except TypeError:
        raw.set_channel_types(mapping)


def edf_headers_by_channel(path: Path) -> dict[str, EdfSignalHeader]:
    """Read and validate original EDF signal metadata for observed channels."""

    rows = read_edf_signal_headers(path)
    by_name: dict[str, EdfSignalHeader] = {}
    for row in rows:
        normalized = normalize_channel_name(row.label)
        if normalized == "EDF ANNOTATIONS":
            continue
        if normalized in by_name:
            raise RuntimeError(f"duplicate EDF signal label {row.label!r} in {path.name}")
        by_name[normalized] = row
    if tuple(by_name) != EXPECTED_CHANNEL_ORDER:
        raise RuntimeError(
            f"direct EDF signal-header order differs from contract in {path.name}"
        )
    return by_name


def update_boolean_run(
    flags: np.ndarray,
    carry: int,
    longest: int,
) -> tuple[int, int]:
    """Update a run length across chunk boundaries.

    Args:
        flags: One chunk's ordered difference-level flatness flags.
        carry: True-run length continuing from the previous chunk.
        longest: Longest true-run observed previously.

    Returns:
        Updated ``(carry, longest)`` difference counts.
    """

    values = np.asarray(flags, dtype=bool)
    if values.size == 0:
        return carry, longest
    if np.all(values):
        carry += int(values.size)
        return carry, max(longest, carry)
    first_false = int(np.flatnonzero(~values)[0])
    merged_prefix = carry + first_false
    suffix = int(values.size - 1 - np.flatnonzero(~values)[-1])
    chunk_longest = longest_true_run(values)
    return suffix, max(longest, merged_prefix, chunk_longest)


def mean_psd_in_band(
    frequencies: np.ndarray,
    psd: np.ndarray,
    low_hz: float,
    high_hz: float,
) -> float | None:
    """Return the finite mean PSD within an inclusive band."""

    mask = (frequencies >= low_hz) & (frequencies <= high_hz)
    values = psd[mask]
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if values.size else None


def finalize_psd_metrics(
    frequencies: np.ndarray | None,
    psd: np.ndarray | None,
) -> dict[str, float | None]:
    """Reduce one channel's accumulated Welch PSD to output summaries."""

    empty = {
        "median_psd_1_100_hz": None,
        "mean_psd_1_4_hz": None,
        "mean_psd_30_45_hz": None,
        "mean_psd_59_61_hz": None,
        "mean_psd_60hz_sidebands": None,
        "line_noise_ratio_60hz": None,
        "low_high_frequency_ratio": None,
    }
    if frequencies is None or psd is None:
        return empty
    broadband = psd[(frequencies >= 1.0) & (frequencies <= 100.0)]
    low = mean_psd_in_band(frequencies, psd, 1.0, 4.0)
    high = mean_psd_in_band(frequencies, psd, 30.0, 45.0)
    line = mean_psd_in_band(frequencies, psd, 59.0, 61.0)
    side_mask = (
        ((frequencies >= 55.0) & (frequencies <= 58.0))
        | ((frequencies >= 62.0) & (frequencies <= 65.0))
    )
    side = float(np.mean(psd[side_mask])) if np.any(side_mask) else None
    return {
        "median_psd_1_100_hz": (
            float(np.median(broadband[np.isfinite(broadband)]))
            if np.any(np.isfinite(broadband)) else None
        ),
        "mean_psd_1_4_hz": low,
        "mean_psd_30_45_hz": high,
        "mean_psd_59_61_hz": line,
        "mean_psd_60hz_sidebands": side,
        "line_noise_ratio_60hz": (
            float(line / side) if line is not None and side is not None and side > 0 else None
        ),
        "low_high_frequency_ratio": (
            float(low / high) if low is not None and high is not None and high > 0 else None
        ),
    }


def scan_raw_file(path: Path, repo_root: Path) -> list[dict[str, Any]]:
    """Scan one EDF fully and return one descriptive row per channel.

    Args:
        path: Raw EDF path.
        repo_root: Repository root used for stable relative paths.

    Returns:
        Per-channel metadata and metric dictionaries.

    Side effects:
        Opens the EDF read-only through MNE and reads every signal sample in
        bounded-memory chunks. The Raw object is closed before return.
    """

    identity = parse_edf_filename(path)
    header_by_name = edf_headers_by_channel(path)
    raw: mne.io.BaseRaw | None = None
    try:
        raw = mne.io.read_raw_edf(path, preload=False, verbose="ERROR")
        apply_channel_types(raw)
        sfreq = float(raw.info["sfreq"])
        n_times = int(raw.n_times)
        n_channels = len(raw.ch_names)
        channel_types = raw.get_channel_types()
        if n_channels != len(EXPECTED_CHANNEL_ORDER):
            raise RuntimeError(f"unexpected MNE channel count in {path.name}")

        unit_statuses = []
        header_rows = []
        mne_scales = np.zeros(n_channels, dtype=float)
        digital_steps = np.zeros(n_channels, dtype=float)
        rail_minima = np.zeros(n_channels, dtype=float)
        rail_maxima = np.zeros(n_channels, dtype=float)
        tolerances = np.zeros(n_channels, dtype=float)
        for index, channel_name in enumerate(raw.ch_names):
            header = header_by_name[normalize_channel_name(channel_name)]
            status = classify_unit_status(channel_types[index], header.physical_dimension)
            scale = mne_scaling_for_unit_status(status)
            step = abs(header.calibration_slope() * scale)
            unit_statuses.append(status)
            header_rows.append(header)
            mne_scales[index] = scale
            digital_steps[index] = step
            rail_minima[index] = header.physical_minimum * scale
            rail_maxima[index] = header.physical_maximum * scale
            tolerances[index] = max(step * 0.5, np.finfo(float).eps)

        chunk_samples = max(1, int(round(FULL_SCAN_CHUNK_SECONDS * sfreq)))
        bounds = deterministic_chunk_bounds(n_times, chunk_samples)
        robust_stride = deterministic_sample_stride(
            n_times, MAXIMUM_ROBUST_SAMPLES_PER_CHANNEL
        )

        finite_counts = np.zeros(n_channels, dtype=np.int64)
        nonfinite_counts = np.zeros(n_channels, dtype=np.int64)
        sums = np.zeros(n_channels, dtype=np.float64)
        sums_of_squares = np.zeros(n_channels, dtype=np.float64)
        minima = np.full(n_channels, np.inf)
        maxima = np.full(n_channels, -np.inf)
        valid_difference_counts = np.zeros(n_channels, dtype=np.int64)
        unchanged_counts = np.zeros(n_channels, dtype=np.int64)
        transition_counts = np.zeros(n_channels, dtype=np.int64)
        maximum_steps = np.zeros(n_channels, dtype=np.float64)
        rail_counts = np.zeros(n_channels, dtype=np.int64)
        flat_carries = np.zeros(n_channels, dtype=np.int64)
        longest_flat_differences = np.zeros(n_channels, dtype=np.int64)
        previous_values = np.full(n_channels, np.nan)
        robust_samples: list[list[np.ndarray]] = [[] for _ in range(n_channels)]
        step_samples: list[list[np.ndarray]] = [[] for _ in range(n_channels)]
        stim_unique_values: list[set[float]] = [set() for _ in range(n_channels)]

        nonstim_indices = [index for index, kind in enumerate(channel_types) if kind != "stim"]
        nperseg = max(2, int(round(WELCH_NPERSEG_SECONDS * sfreq)))
        noverlap = int(round(nperseg * WELCH_OVERLAP_PROPORTION))
        psd_frequencies: np.ndarray | None = None
        psd_weighted_sum: np.ndarray | None = None
        psd_segment_counts = np.zeros(n_channels, dtype=np.int64)
        psd_coverage_samples = np.zeros(n_channels, dtype=np.int64)

        for start, stop in bounds:
            data = raw.get_data(start=start, stop=stop)
            if data.shape != (n_channels, stop - start):
                raise RuntimeError(f"unexpected MNE data shape in {path.name}")
            finite_mask = np.isfinite(data)
            finite_counts += finite_mask.sum(axis=1)
            nonfinite_counts += (~finite_mask).sum(axis=1)
            safe_data = np.where(finite_mask, data, 0.0)
            sums += safe_data.sum(axis=1)
            sums_of_squares += np.square(safe_data).sum(axis=1)
            for index in range(n_channels):
                finite_values = data[index, finite_mask[index]]
                if finite_values.size:
                    minima[index] = min(minima[index], float(np.min(finite_values)))
                    maxima[index] = max(maxima[index], float(np.max(finite_values)))
                    rail_counts[index] += np.count_nonzero(
                        (np.abs(finite_values - rail_minima[index]) <= tolerances[index])
                        | (np.abs(finite_values - rail_maxima[index]) <= tolerances[index])
                    )
                if channel_types[index] == "stim":
                    stim_unique_values[index].update(float(value) for value in finite_values)

                row = data[index]
                if start == 0:
                    difference_values = np.diff(row)
                else:
                    difference_values = np.diff(
                        np.concatenate(([previous_values[index]], row))
                    )
                if row.size:
                    previous_values[index] = row[-1]
                difference_finite = np.isfinite(difference_values)
                valid_differences = difference_values[difference_finite]
                absolute_steps = np.abs(valid_differences)
                valid_difference_counts[index] += valid_differences.size
                unchanged_counts[index] += np.count_nonzero(valid_differences == 0.0)
                transition_counts[index] += np.count_nonzero(valid_differences != 0.0)
                if absolute_steps.size:
                    maximum_steps[index] = max(
                        maximum_steps[index], float(np.max(absolute_steps))
                    )
                    step_samples[index].append(
                        absolute_steps[::robust_stride].astype(float, copy=False)
                    )
                near_flat = np.zeros(difference_values.size, dtype=bool)
                near_flat[difference_finite] = (
                    absolute_steps <= tolerances[index]
                )
                flat_carries[index], longest_flat_differences[index] = update_boolean_run(
                    near_flat,
                    int(flat_carries[index]),
                    int(longest_flat_differences[index]),
                )

            robust_local_indices = np.flatnonzero(
                np.arange(start, stop, dtype=np.int64) % robust_stride == 0
            )
            if robust_local_indices.size:
                sampled = data[:, robust_local_indices]
                for index in range(n_channels):
                    robust_samples[index].append(sampled[index].astype(float, copy=False))

            # Welch runs on complete fixed-length segments only. EDFs are finite
            # in practice; if a channel ever contains nonfinite data, its PSD is
            # omitted while dropout metrics remain available.
            if stop - start >= nperseg and nonstim_indices:
                for index in nonstim_indices:
                    row = data[index]
                    if not np.all(np.isfinite(row)):
                        continue
                    frequencies, psd = welch(
                        row,
                        fs=sfreq,
                        nperseg=nperseg,
                        noverlap=noverlap,
                        detrend="constant",
                        scaling="density",
                    )
                    segments = 1 + (row.size - nperseg) // (nperseg - noverlap)
                    if psd_frequencies is None:
                        psd_frequencies = frequencies
                        psd_weighted_sum = np.zeros((n_channels, psd.size), dtype=float)
                    if not np.array_equal(frequencies, psd_frequencies):
                        raise RuntimeError("inconsistent Welch frequency grid")
                    assert psd_weighted_sum is not None
                    psd_weighted_sum[index] += psd * segments
                    psd_segment_counts[index] += segments
                    psd_coverage_samples[index] += row.size

        rows: list[dict[str, Any]] = []
        for index, channel_name in enumerate(raw.ch_names):
            header = header_rows[index]
            status = unit_statuses[index]
            finite_count = int(finite_counts[index])
            samples = (
                np.concatenate(robust_samples[index])
                if robust_samples[index]
                else np.array([], dtype=float)
            )
            samples = samples[np.isfinite(samples)]
            steps = (
                np.concatenate(step_samples[index])
                if step_samples[index]
                else np.array([], dtype=float)
            )
            if samples.size:
                median = float(np.median(samples))
                mad = float(np.median(np.abs(samples - median)))
                q01, q99 = np.quantile(samples, [0.01, 0.99])
            else:
                median = mad = q01 = q99 = np.nan
            mean = sums[index] / finite_count if finite_count else np.nan
            variance = (
                max(0.0, sums_of_squares[index] / finite_count - mean * mean)
                if finite_count else np.nan
            )
            median_step = float(np.median(steps)) if steps.size else np.nan
            p999_step = float(np.quantile(steps, 0.999)) if steps.size else np.nan
            step_denominator = max(
                median_step if np.isfinite(median_step) else 0.0,
                digital_steps[index],
                np.finfo(float).eps,
            )
            continuous_applicable = channel_types[index] != "stim"
            averaged_psd = None
            if (
                continuous_applicable
                and psd_weighted_sum is not None
                and psd_segment_counts[index] > 0
            ):
                averaged_psd = psd_weighted_sum[index] / psd_segment_counts[index]
            spectral = finalize_psd_metrics(psd_frequencies, averaged_psd)
            mne_channel = raw.info["chs"][index]
            mne_unit_value = mne_channel.get("unit")
            mne_unit_code = int(mne_unit_value) if mne_unit_value is not None else None
            mne_unit_label = str(mne_unit_value) if mne_unit_value is not None else ""
            longest_samples = (
                int(longest_flat_differences[index] + 1) if n_times else 0
            )

            amplitude_values: dict[str, Any]
            if continuous_applicable:
                amplitude_values = {
                    "robust_median": median if np.isfinite(median) else None,
                    "robust_mad": mad if np.isfinite(mad) else None,
                    "robust_scale_mad": 1.4826 * mad if np.isfinite(mad) else None,
                    "robust_q01": q01 if np.isfinite(q01) else None,
                    "robust_q99": q99 if np.isfinite(q99) else None,
                    "robust_peak_to_peak_q01_q99": (
                        q99 - q01 if np.isfinite(q01) and np.isfinite(q99) else None
                    ),
                    "full_minimum": minima[index] if finite_count else None,
                    "full_maximum": maxima[index] if finite_count else None,
                    "full_peak_to_peak": (
                        maxima[index] - minima[index] if finite_count else None
                    ),
                    "full_mean": mean if np.isfinite(mean) else None,
                    "full_variance": variance if np.isfinite(variance) else None,
                }
            else:
                amplitude_values = {
                    name: None for name in (
                        "robust_median", "robust_mad", "robust_scale_mad",
                        "robust_q01", "robust_q99",
                        "robust_peak_to_peak_q01_q99", "full_minimum",
                        "full_maximum", "full_peak_to_peak", "full_mean",
                        "full_variance",
                    )
                }

            rows.append(
                {
                    **identity,
                    "file_path": path.relative_to(repo_root).as_posix(),
                    "channel_index": index,
                    "raw_channel_name": channel_name,
                    "normalized_channel_name": normalize_channel_name(channel_name),
                    "configured_channel_type": channel_types[index],
                    "channel_semantic_group": (
                        "scalp_eeg" if channel_types[index] == "eeg"
                        else "auxiliary" if channel_types[index] in {"eog", "emg"}
                        else "trigger"
                    ),
                    "edf_original_physical_dimension": header.physical_dimension,
                    "edf_original_unit_display": header.physical_dimension or "blank/unknown",
                    "edf_transducer_type": header.transducer_type,
                    "edf_prefiltering": header.prefiltering,
                    "edf_physical_minimum": header.physical_minimum,
                    "edf_physical_maximum": header.physical_maximum,
                    "edf_digital_minimum": header.digital_minimum,
                    "edf_digital_maximum": header.digital_maximum,
                    "edf_samples_per_record": header.samples_per_record,
                    "edf_calibration_slope_physical_units_per_digital_code": header.calibration_slope(),
                    "edf_calibration_offset_physical_units": header.calibration_offset(),
                    "edf_to_mne_memory_scale": mne_scales[index],
                    "digital_step_in_memory_units": digital_steps[index],
                    "rail_minimum_in_memory_units": rail_minima[index],
                    "rail_maximum_in_memory_units": rail_maxima[index],
                    "mne_unit_code": mne_unit_code,
                    "mne_unit_label": mne_unit_label,
                    "mne_info_cal": float(mne_channel.get("cal", np.nan)),
                    "mne_info_range": float(mne_channel.get("range", np.nan)),
                    **status.as_dict(),
                    "continuous_amplitude_metrics_applicable": continuous_applicable,
                    "sampling_frequency_hz": sfreq,
                    "recording_n_times": n_times,
                    "recording_duration_seconds": n_times / sfreq,
                    "coverage_strategy": "full_duration_chunked_scan",
                    "full_scan_chunk_seconds": FULL_SCAN_CHUNK_SECONDS,
                    "full_scan_chunk_count": len(bounds),
                    "full_scan_samples": n_times,
                    "full_scan_coverage_fraction": 1.0 if n_times else 0.0,
                    "robust_sample_stride": robust_stride,
                    "robust_sample_count": int(samples.size),
                    "robust_sample_coverage_first_sample": 0 if samples.size else None,
                    "robust_sample_coverage_last_sample": (
                        ((n_times - 1) // robust_stride) * robust_stride
                        if samples.size else None
                    ),
                    "nonfinite_sample_count": int(nonfinite_counts[index]),
                    "finite_sample_proportion": (
                        finite_count / n_times if n_times else None
                    ),
                    **amplitude_values,
                    "unchanged_difference_count": int(unchanged_counts[index]),
                    "unchanged_difference_proportion": (
                        unchanged_counts[index] / valid_difference_counts[index]
                        if valid_difference_counts[index] else None
                    ),
                    "near_flat_tolerance_in_memory_units": tolerances[index],
                    "longest_near_flat_run_samples": longest_samples,
                    "longest_near_flat_run_seconds": longest_samples / sfreq,
                    "rail_metadata_available": True,
                    "rail_sample_count": int(rail_counts[index]),
                    "rail_sample_proportion": (
                        rail_counts[index] / finite_count if finite_count else None
                    ),
                    "maximum_absolute_step": (
                        maximum_steps[index] if valid_difference_counts[index] else None
                    ),
                    "sampled_median_absolute_step": (
                        median_step if np.isfinite(median_step) else None
                    ),
                    "sampled_p999_absolute_step": (
                        p999_step if np.isfinite(p999_step) else None
                    ),
                    "extreme_step_ratio": (
                        maximum_steps[index] / step_denominator
                        if valid_difference_counts[index] else None
                    ),
                    "welch_nperseg": nperseg if continuous_applicable else None,
                    "welch_noverlap": noverlap if continuous_applicable else None,
                    "welch_segment_count": int(psd_segment_counts[index]),
                    "psd_coverage_samples": int(psd_coverage_samples[index]),
                    "psd_coverage_fraction": (
                        psd_coverage_samples[index] / n_times if n_times else 0.0
                    ),
                    **spectral,
                    "stim_transition_count": (
                        int(transition_counts[index]) if channel_types[index] == "stim" else None
                    ),
                    "stim_unique_digital_state_count": (
                        len(stim_unique_values[index]) if channel_types[index] == "stim" else None
                    ),
                    "stim_unique_digital_states": (
                        ";".join(f"{value:g}" for value in sorted(stim_unique_values[index]))
                        if channel_types[index] == "stim" else ""
                    ),
                    "final_bad_channel_label": "",
                    "automatic_bad_channel_assignment_performed": False,
                }
            )
        return rows
    finally:
        if raw is not None and hasattr(raw, "close"):
            raw.close()


def add_descriptive_ranks(metrics: pd.DataFrame) -> pd.DataFrame:
    """Add within-type ranks without comparing unknown auxiliary amplitudes."""

    ranked = metrics.copy()
    scale_free_metrics = (
        "line_noise_ratio_60hz",
        "longest_near_flat_run_seconds",
        "rail_sample_proportion",
        "extreme_step_ratio",
    )
    for metric in scale_free_metrics:
        ranked[f"rank_{metric}_descending_within_type"] = ranked.groupby(
            "configured_channel_type", dropna=False
        )[metric].rank(method="min", ascending=False)
        ranked[f"percentile_{metric}_within_type"] = ranked.groupby(
            "configured_channel_type", dropna=False
        )[metric].rank(method="average", pct=True)

    # Absolute scale is ranked only where the unit contract permits physical
    # cross-file magnitude comparison. Unknown EOG/EMG values remain blank.
    comparable = ranked["absolute_cross_file_magnitude_interpretable"].astype(bool)
    ranked["rank_robust_scale_mad_descending_calibrated_type"] = np.nan
    ranked["percentile_robust_scale_mad_calibrated_type"] = np.nan
    subset = ranked.loc[comparable]
    ranks = subset.groupby("configured_channel_type")["robust_scale_mad"].rank(
        method="min", ascending=False
    )
    percentiles = subset.groupby("configured_channel_type")["robust_scale_mad"].rank(
        method="average", pct=True
    )
    ranked.loc[comparable, "rank_robust_scale_mad_descending_calibrated_type"] = ranks
    ranked.loc[comparable, "percentile_robust_scale_mad_calibrated_type"] = percentiles
    return ranked


def candidate_reasons_for_row(row: pd.Series) -> list[str]:
    """Return descriptive review reasons without assigning a bad channel."""

    reasons: list[str] = []
    channel_type = str(row["configured_channel_type"])
    if int(row["nonfinite_sample_count"]) > 0:
        reasons.append("nonfinite_samples_present")
    if channel_type == "stim":
        if int(row["stim_transition_count"]) == 0:
            reasons.append("stim_zero_transitions_factual_exception")
        return reasons
    robust_scale = row.get("robust_scale_mad")
    if pd.notna(robust_scale) and float(robust_scale) == 0.0:
        reasons.append("zero_robust_scale")
    if float(row["longest_near_flat_run_seconds"]) >= 5.0:
        reasons.append("near_flat_run_at_least_5_seconds")
    if float(row["unchanged_difference_proportion"]) >= 0.05:
        reasons.append("unchanged_sample_proportion_at_least_5_percent")
    if float(row["rail_sample_proportion"]) >= 0.001:
        reasons.append("edf_rail_samples_at_least_0.1_percent")

    for metric in (
        "line_noise_ratio_60hz",
        "longest_near_flat_run_seconds",
        "rail_sample_proportion",
        "extreme_step_ratio",
    ):
        value = row.get(metric)
        percentile = row.get(f"percentile_{metric}_within_type")
        if pd.notna(value) and float(value) > 0 and pd.notna(percentile):
            if float(percentile) >= CANDIDATE_EXTREME_PERCENTILE:
                reasons.append(f"extreme_within_type_rank_{metric}")

    # Only calibrated channels may receive absolute-scale rank reasons.
    scale_percentile = row.get("percentile_robust_scale_mad_calibrated_type")
    if bool(row["absolute_cross_file_magnitude_interpretable"]) and pd.notna(scale_percentile):
        if float(scale_percentile) >= CANDIDATE_EXTREME_PERCENTILE:
            reasons.append("extreme_high_calibrated_robust_scale_rank")
        if float(scale_percentile) <= 1.0 - CANDIDATE_EXTREME_PERCENTILE:
            reasons.append("extreme_low_calibrated_robust_scale_rank")
    return reasons


def build_candidate_table(metrics: pd.DataFrame) -> pd.DataFrame:
    """Build a compact descriptive review queue from per-channel evidence."""

    candidates = metrics.copy()
    candidates["candidate_reasons"] = candidates.apply(
        lambda row: ";".join(candidate_reasons_for_row(row)), axis=1
    )
    candidates = candidates[candidates["candidate_reasons"].ne("")].copy()
    candidates["candidate_reason_count"] = candidates["candidate_reasons"].str.count(";") + 1
    candidates["candidate_priority"] = np.where(
        candidates["candidate_reasons"].str.contains(
            "nonfinite|zero_robust|at_least|stim_zero", regex=True
        ),
        "hard_descriptive_exception",
        "extreme_rank_review",
    )
    candidates["final_bad_channel_label"] = ""
    candidates["candidate_is_not_bad_channel_assignment"] = True
    return candidates.sort_values(
        ["candidate_priority", "candidate_reason_count", "eeg_recording_id", "channel_index"],
        ascending=[True, False, True, True],
        kind="mergesort",
    )


def build_file_summary(metrics: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    """Reduce per-channel metrics and candidate rows to one row per EDF."""

    identity_columns = [
        "eeg_recording_id", "eeg_recording_id_padded", "source_filename",
        "file_path", "file_role", "split_part",
    ]
    summaries: list[dict[str, Any]] = []
    candidate_counts = candidates.groupby("source_filename").size().to_dict()
    hard_counts = (
        candidates[candidates["candidate_priority"].eq("hard_descriptive_exception")]
        .groupby("source_filename").size().to_dict()
    )
    for source_filename, group in metrics.groupby("source_filename", sort=False):
        row = {column: group.iloc[0][column] for column in identity_columns}
        row.update(
            {
                "channel_row_count": len(group),
                "eeg_channel_count": int(group["configured_channel_type"].eq("eeg").sum()),
                "eog_channel_count": int(group["configured_channel_type"].eq("eog").sum()),
                "emg_channel_count": int(group["configured_channel_type"].eq("emg").sum()),
                "stim_channel_count": int(group["configured_channel_type"].eq("stim").sum()),
                "unknown_unit_auxiliary_channel_count": int(
                    group["value_status"].eq("arbitrary_unknown_acquisition_units").sum()
                ),
                "nonfinite_sample_count_all_channels": int(group["nonfinite_sample_count"].sum()),
                "candidate_channel_count": int(candidate_counts.get(source_filename, 0)),
                "hard_descriptive_exception_channel_count": int(hard_counts.get(source_filename, 0)),
                "final_bad_channel_count": 0,
                "automatic_bad_channel_assignment_performed": False,
                "full_duration_scan_complete": bool(
                    np.allclose(group["full_scan_coverage_fraction"], 1.0)
                ),
            }
        )
        summaries.append(row)
    return pd.DataFrame(summaries)


def sha256_file(path: Path) -> str:
    """Return SHA-256 for a tracked code/input file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_text(repo_root: Path, arguments: list[str]) -> str:
    """Return compact Git provenance text without changing repository state."""

    result = subprocess.run(
        ["git", *arguments],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def input_inventory_digest(paths: list[Path], repo_root: Path) -> str:
    """Hash stable filenames, sizes, and mtimes without hashing 32 GB of raw data."""

    rows = [
        {
            "path": path.relative_to(repo_root).as_posix(),
            "size_bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for path in paths
    ]
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def package_versions() -> dict[str, str]:
    """Return relevant runtime versions for the run manifest."""

    import scipy

    return {
        "python": platform.python_version(),
        "mne": mne.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def write_summary(
    path: Path,
    metrics: pd.DataFrame,
    files: pd.DataFrame,
    candidates: pd.DataFrame,
) -> None:
    """Write a concise Markdown inventory of unit and candidate findings."""

    status_counts = metrics["value_status"].value_counts().to_dict()
    type_counts = metrics["configured_channel_type"].value_counts().to_dict()
    candidate_type_counts = candidates["configured_channel_type"].value_counts().to_dict()
    reason_counts = Counter(
        reason
        for value in candidates["candidate_reasons"]
        for reason in str(value).split(";")
        if reason
    )
    lines = [
        "# Unit-aware per-channel raw QC v1",
        "",
        "This is descriptive raw evidence only. No signal was filtered, re-referenced, marked bad, interpolated, ICA-cleaned, converted to CSD, epoched, or written as a cleaned derivative.",
        "",
        "## Output shape",
        "",
        f"- EDF files: {len(files)}",
        f"- channel rows: {len(metrics)}",
        f"- candidate review rows: {len(candidates)}",
        f"- hard descriptive exception rows: {int(candidates['candidate_priority'].eq('hard_descriptive_exception').sum()) if len(candidates) else 0}",
        f"- channel types: `{json.dumps(type_counts, sort_keys=True)}`",
        f"- value statuses: `{json.dumps(status_counts, sort_keys=True)}`",
        f"- candidate channel types: `{json.dumps(candidate_type_counts, sort_keys=True)}`",
        "",
        "## Unit contract",
        "",
        "Scalp EEG has calibrated EDF microvolt metadata and MNE returns volts. HEO, VEO, EMG-L, and EMG-A remain unknown acquisition units because their original EDF dimensions are blank. Trigger is digital and is not summarized as a continuous amplitude or PSD channel. Unknown-unit auxiliary absolute amplitudes are not ranked across files.",
        "",
        "## Candidate reasons",
        "",
    ]
    if reason_counts:
        lines.extend(f"- {reason}: {count}" for reason, count in reason_counts.most_common())
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "Candidate rows are a compact manual-review queue, not bad-channel labels. Extreme-rank reasons use scale-free morphology/frequency metrics; calibrated amplitude ranks are restricted to channels whose EDF units support them.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the all-file channel-QC audit and write local versioned evidence."""

    args = parse_args()
    repo_root = repo_root_from_script()
    raw_dir = repo_root / RAW_EEG_DIR
    output_dir = repo_root / OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = discover_edf_files(raw_dir, args.max_files)

    all_rows: list[dict[str, Any]] = []
    for file_index, path in enumerate(paths, start=1):
        print(f"[{file_index:02d}/{len(paths):02d}] scanning {path.name}", flush=True)
        all_rows.extend(scan_raw_file(path, repo_root))
    metrics = add_descriptive_ranks(pd.DataFrame(all_rows))
    candidates = build_candidate_table(metrics)
    file_summary = build_file_summary(metrics, candidates)
    stim_summary = metrics[metrics["configured_channel_type"].eq("stim")].copy()

    metrics.to_csv(output_dir / CHANNEL_METRICS_FILENAME, index=False)
    file_summary.to_csv(output_dir / FILE_SUMMARY_FILENAME, index=False)
    candidates.to_csv(output_dir / CANDIDATE_FILENAME, index=False)
    stim_summary.to_csv(output_dir / STIM_SUMMARY_FILENAME, index=False)
    write_summary(output_dir / SUMMARY_FILENAME, metrics, file_summary, candidates)

    script_path = Path(__file__).resolve()
    helper_path = script_path.parent / "channel_qc.py"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_timestamp_local": datetime.now().astimezone().isoformat(),
        "script_path": script_path.relative_to(repo_root).as_posix(),
        "script_sha256": sha256_file(script_path),
        "helper_path": helper_path.relative_to(repo_root).as_posix(),
        "helper_sha256": sha256_file(helper_path),
        "git_commit": git_text(repo_root, ["rev-parse", "HEAD"]),
        "git_status_short_at_run": git_text(repo_root, ["status", "--short"]),
        "runtime_versions": package_versions(),
        "input_file_count": len(paths),
        "input_inventory_filename_size_mtime_sha256": input_inventory_digest(paths, repo_root),
        "raw_content_hash_note": "Raw EDF content was not hashed because the local corpus is approximately 32 GB; filenames, sizes, and nanosecond mtimes are hashed for run identity.",
        "full_scan_chunk_seconds": FULL_SCAN_CHUNK_SECONDS,
        "maximum_robust_samples_per_channel": MAXIMUM_ROBUST_SAMPLES_PER_CHANNEL,
        "welch_nperseg_seconds": WELCH_NPERSEG_SECONDS,
        "welch_overlap_proportion": WELCH_OVERLAP_PROPORTION,
        "candidate_extreme_percentile": CANDIDATE_EXTREME_PERCENTILE,
        "channel_row_count": len(metrics),
        "file_summary_row_count": len(file_summary),
        "candidate_row_count": len(candidates),
        "stim_summary_row_count": len(stim_summary),
        "value_status_counts": metrics["value_status"].value_counts().to_dict(),
        "safety": {
            "raw_files_modified": False,
            "filtering": False,
            "notch_filtering": False,
            "rereferencing": False,
            "automatic_bad_channel_assignment": False,
            "interpolation": False,
            "ica": False,
            "csd": False,
            "cleaned_derivative_writing": False,
            "epoch_construction": False,
            "participant_file_trial_decisions": False,
        },
        "outputs": sorted(
            path.name for path in output_dir.iterdir() if path.is_file()
        ),
    }
    (output_dir / RUN_MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"Wrote {len(metrics)} channel rows, {len(candidates)} candidate rows, "
        f"and {len(file_summary)} file rows to {output_dir.relative_to(repo_root)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
