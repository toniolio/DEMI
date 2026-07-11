"""Synthetic tests for DEMI unit-aware per-channel raw QC helpers.

No test opens private/raw DEMI data.  Unit classifications are small synthetic
cases, and the EDF parser fixture contains only a fabricated fixed-width
header with no signal samples.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_ANALYSIS_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_ANALYSIS_DIR))

from channel_qc import (  # noqa: E402
    classify_unit_status,
    deterministic_chunk_bounds,
    deterministic_sample_stride,
    longest_true_run,
    mne_scaling_for_unit_status,
    read_edf_signal_headers,
    summarize_signal_metrics,
    summarize_spectral_metrics,
)


def _field(values: list[str], width: int) -> bytes:
    """Encode one repeated synthetic EDF signal-header field."""

    return b"".join(value.encode("ascii").ljust(width) for value in values)


def _synthetic_edf_header() -> bytes:
    """Return a complete two-signal EDF header without any sample data."""

    signal_count = 2
    fixed = bytearray(b" " * 256)
    fixed[0:8] = b"0       "
    fixed[184:192] = str(256 + 256 * signal_count).encode("ascii").ljust(8)
    fixed[252:256] = str(signal_count).encode("ascii").ljust(4)
    signal = b"".join(
        [
            _field(["Cz", "Trigger"], 16),
            _field(["", ""], 80),
            _field(["uV", ""], 8),
            _field(["-100", "65280"], 8),
            _field(["100", "65535"], 8),
            _field(["-32768", "-32768"], 8),
            _field(["32767", "32767"], 8),
            _field(["", ""], 80),
            _field(["1000", "1000"], 8),
            _field(["", ""], 32),
        ]
    )
    assert len(signal) == signal_count * 256
    return bytes(fixed) + signal


def test_microvolt_eeg_is_calibrated_and_mne_scaled_to_volts() -> None:
    """Declared microvolts support calibrated voltage QC after MNE scaling."""

    status = classify_unit_status("eeg", "µV")
    assert status.value_status == "mne_scaled_calibrated_voltage"
    assert status.in_memory_value_unit == "V"
    assert status.values_are_calibrated_physical_units
    assert status.values_are_mne_scaled
    assert status.absolute_cross_file_magnitude_interpretable
    assert mne_scaling_for_unit_status(status) == pytest.approx(1e-6)


@pytest.mark.parametrize("channel_type", ["eog", "emg"])
@pytest.mark.parametrize("unknown_unit", ["", "n/a", "unknown", None])
def test_unknown_unit_auxiliary_is_never_labelled_calibrated_si(
    channel_type: str,
    unknown_unit: object,
) -> None:
    """A type override must not invent calibration for auxiliary channels."""

    status = classify_unit_status(channel_type, unknown_unit)
    assert status.value_status == "arbitrary_unknown_acquisition_units"
    assert status.in_memory_value_unit == "unknown_acquisition_unit"
    assert not status.values_are_calibrated_physical_units
    assert not status.values_are_mne_scaled
    assert not status.absolute_cross_file_magnitude_interpretable
    assert "within_channel_file" in status.metric_semantics


def test_trigger_is_digital_not_a_continuous_physical_channel() -> None:
    """Trigger semantics remain digital even if a header declared a unit."""

    for original_unit in ("", "uV"):
        status = classify_unit_status("stim", original_unit)
        assert status.value_status == "digital_stim_values"
        assert status.in_memory_value_unit == "digital_code"
        assert not status.values_are_calibrated_physical_units
        assert not status.absolute_cross_file_magnitude_interpretable
        assert status.metric_semantics == "event_transitions_and_digital_states_only"


def test_synthetic_edf_header_preserves_ranges_and_calibration(tmp_path: Path) -> None:
    """The direct parser retains original dimensions and digital/physical ranges."""

    path = tmp_path / "synthetic.edf"
    path.write_bytes(_synthetic_edf_header())
    rows = read_edf_signal_headers(path)

    assert [row.label for row in rows] == ["Cz", "Trigger"]
    assert rows[0].physical_dimension == "uV"
    assert rows[1].physical_dimension == ""
    assert rows[0].digital_minimum == -32768
    assert rows[0].digital_maximum == 32767
    assert rows[0].calibration_slope() == pytest.approx(200 / 65535)
    assert rows[0].as_qc_dict()[
        "edf_calibration_slope_physical_units_per_digital_code"
    ] == pytest.approx(200 / 65535)


def test_inconsistent_edf_header_length_fails(tmp_path: Path) -> None:
    """A malformed declared header size cannot be silently parsed."""

    payload = bytearray(_synthetic_edf_header())
    payload[184:192] = b"999     "
    path = tmp_path / "bad.edf"
    path.write_bytes(payload)
    with pytest.raises(ValueError, match="header byte count is inconsistent"):
        read_edf_signal_headers(path)


def test_deterministic_chunks_cover_every_sample_once() -> None:
    """Chunk bounds are ordered, gap-free, non-overlapping, and complete."""

    bounds = deterministic_chunk_bounds(23, 7)
    assert bounds == [(0, 7), (7, 14), (14, 21), (21, 23)]
    covered = [sample for start, stop in bounds for sample in range(start, stop)]
    assert covered == list(range(23))
    assert deterministic_chunk_bounds(0, 7) == []
    assert deterministic_sample_stride(1_000_001, 200_000) == 6


def test_longest_true_run_handles_edges_and_empty_values() -> None:
    """Flat-run helper counts runs at the start, middle, and end."""

    assert longest_true_run(np.array([], dtype=bool)) == 0
    assert longest_true_run(np.array([False, False])) == 0
    assert longest_true_run(np.array([True, True, False, True])) == 2
    assert longest_true_run(np.array([False, True, True, True])) == 3


def test_flat_signal_metrics_surface_zero_scale_and_full_flat_run() -> None:
    """A fully flat signal remains descriptive rather than a bad-channel label."""

    metrics = summarize_signal_metrics(np.ones(1_000), 1_000.0, digital_step_in_memory_units=1.0)
    assert metrics["robust_scale_mad"] == 0.0
    assert metrics["variance"] == 0.0
    assert metrics["unchanged_difference_proportion"] == 1.0
    assert metrics["longest_near_flat_run_samples"] == 1_000
    assert metrics["longest_near_flat_run_seconds"] == pytest.approx(1.0)


def test_noisy_and_normal_signals_have_nonzero_robust_scale() -> None:
    """Deterministic noise and a normal waveform do not look flat."""

    rng = np.random.default_rng(97)
    noisy = rng.normal(size=10_000)
    normal = np.sin(2 * np.pi * 10 * np.arange(10_000) / 1_000)
    noisy_metrics = summarize_signal_metrics(noisy, 1_000.0)
    normal_metrics = summarize_signal_metrics(normal, 1_000.0)
    assert noisy_metrics["robust_scale_mad"] > 0
    assert normal_metrics["robust_scale_mad"] > 0
    assert noisy_metrics["unchanged_difference_proportion"] == 0.0
    assert normal_metrics["nonfinite_count"] == 0


def test_stepped_signal_has_larger_extreme_step_than_normal_signal() -> None:
    """An injected step is visible in absolute and normalized change metrics."""

    time = np.arange(20_000) / 1_000
    normal = np.sin(2 * np.pi * 5 * time)
    stepped = normal.copy()
    stepped[10_000:] += 100.0
    normal_metrics = summarize_signal_metrics(normal, 1_000.0)
    stepped_metrics = summarize_signal_metrics(stepped, 1_000.0)
    assert stepped_metrics["maximum_absolute_step"] > 50
    assert stepped_metrics["maximum_absolute_step"] > normal_metrics["maximum_absolute_step"]
    assert stepped_metrics["extreme_step_ratio"] > normal_metrics["extreme_step_ratio"]


def test_clipped_signal_counts_metadata_supported_rail_samples() -> None:
    """Rail counts use explicit range metadata and do not infer calibration."""

    signal = np.linspace(-1.0, 1.0, 1_000)
    signal[:50] = -1.0
    signal[-75:] = 1.0
    metrics = summarize_signal_metrics(
        signal,
        1_000.0,
        digital_step_in_memory_units=2.0 / 65_535,
        rail_minimum_in_memory_units=-1.0,
        rail_maximum_in_memory_units=1.0,
    )
    assert metrics["rail_sample_count"] >= 125
    assert metrics["rail_sample_proportion"] >= 0.125


def test_nonfinite_samples_are_counted_as_dropouts() -> None:
    """NaN and infinite values are reported without contaminating finite metrics."""

    signal = np.array([0.0, 1.0, np.nan, 2.0, np.inf, 3.0])
    metrics = summarize_signal_metrics(signal, 1_000.0)
    assert metrics["nonfinite_count"] == 2
    assert metrics["finite_proportion"] == pytest.approx(4 / 6)
    assert metrics["minimum"] == 0.0
    assert metrics["maximum"] == 3.0


def test_spectral_metrics_detect_a_60_hz_component() -> None:
    """The scale-free line ratio distinguishes 60-Hz-rich synthetic data."""

    sfreq = 1_000.0
    time = np.arange(30_000) / sfreq
    base = np.sin(2 * np.pi * 10 * time)
    with_line = base + 4 * np.sin(2 * np.pi * 60 * time)
    base_metrics = summarize_spectral_metrics(base, sfreq)
    line_metrics = summarize_spectral_metrics(with_line, sfreq)
    assert line_metrics["line_noise_ratio_60hz"] > base_metrics["line_noise_ratio_60hz"]
    assert line_metrics["mean_psd_59_61_hz"] > base_metrics["mean_psd_59_61_hz"]
