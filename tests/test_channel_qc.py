"""Synthetic tests for DEMI unit-aware per-channel raw QC helpers.

No test opens private/raw DEMI data.  Unit classifications are small synthetic
cases, and the EDF parser fixture contains only a fabricated fixed-width
header with no signal samples.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_ANALYSIS_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_ANALYSIS_DIR))

from channel_qc import (  # noqa: E402
    classify_unit_status,
    mne_scaling_for_unit_status,
    read_edf_signal_headers,
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
