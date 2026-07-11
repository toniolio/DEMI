"""Pure unit and calibration contracts for DEMI raw-channel QC.

This module exists because an EDF physical-dimension field and an MNE channel
unit are not interchangeable evidence.  In the DEMI recordings, the 32 scalp
channels declare microvolts in the EDF header, while HEO, VEO, EMG-L, EMG-A,
and Trigger have blank physical-dimension fields.  MNE can still expose those
blank-dimension samples through ``Raw.get_data()`` and can assign a voltage
unit when channel types are changed in memory; neither action establishes an
acquisition calibration that is absent from the source EDF.

The helpers here:

* parse the fixed-width EDF signal header without loading signal samples;
* retain original physical/digital ranges and calibration coefficients;
* classify values as calibrated/MNE-scaled voltage, calibrated physical
  values, unknown acquisition units, or digital/stim values; and
* state whether absolute cross-file magnitude comparisons are supported.

The module performs no filtering, referencing, bad-channel assignment,
interpolation, ICA, CSD, event repair, epoch construction, or file writing.
``read_edf_signal_headers`` opens an EDF read-only and reads only its header.
It is used by the separate channel-QC driver, which keeps generated evidence
under the ignored ``_Data/`` tree.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import BinaryIO


UNKNOWN_EDF_UNITS = {"", "N/A", "NA", "NONE", "UNKNOWN", "?"}
VOLTAGE_UNIT_SCALE_TO_VOLTS = {
    "V": 1.0,
    "MV": 1e-3,
    "UV": 1e-6,
    "NV": 1e-9,
}


@dataclass(frozen=True)
class EdfSignalHeader:
    """One signal's original fixed-width EDF header fields.

    Attributes:
        signal_index: Zero-based position in the EDF signal header.
        label: Original EDF signal label.
        transducer_type: Original transducer description, often blank.
        physical_dimension: Original physical dimension, such as ``uV``.
        physical_minimum: Physical value represented by ``digital_minimum``.
        physical_maximum: Physical value represented by ``digital_maximum``.
        digital_minimum: Smallest stored digital code.
        digital_maximum: Largest stored digital code.
        prefiltering: Original EDF prefiltering description, often blank.
        samples_per_record: Samples for this signal in each EDF data record.
        reserved: Per-signal EDF reserved text.
    """

    signal_index: int
    label: str
    transducer_type: str
    physical_dimension: str
    physical_minimum: float
    physical_maximum: float
    digital_minimum: int
    digital_maximum: int
    prefiltering: str
    samples_per_record: int
    reserved: str

    def calibration_slope(self) -> float:
        """Return physical units per digital code from EDF range fields."""

        digital_span = self.digital_maximum - self.digital_minimum
        if digital_span <= 0:
            raise ValueError(
                f"EDF signal {self.label!r} has a non-positive digital range"
            )
        return (self.physical_maximum - self.physical_minimum) / digital_span

    def calibration_offset(self) -> float:
        """Return the physical-unit intercept for ``physical=slope*digital+offset``."""

        return self.physical_minimum - self.calibration_slope() * self.digital_minimum

    def as_qc_dict(self) -> dict[str, object]:
        """Return header and derived calibration fields for a QC table row."""

        row = asdict(self)
        row["edf_calibration_slope_physical_units_per_digital_code"] = (
            self.calibration_slope()
        )
        row["edf_calibration_offset_physical_units"] = self.calibration_offset()
        return row


@dataclass(frozen=True)
class UnitStatus:
    """Scientific interpretation of values returned for one raw channel.

    Attributes:
        value_status: Stable machine-readable classification.
        edf_original_unit_status: Whether the EDF dimension is known.
        edf_original_unit_normalized: Normalized original dimension.
        in_memory_value_unit: Unit to use when labelling MNE-returned values.
        values_are_calibrated_physical_units: Whether the EDF declares a
            physical dimension and range calibration.
        values_are_mne_scaled: Whether MNE converts a declared voltage unit to
            volts in memory.
        absolute_cross_file_magnitude_interpretable: Whether an absolute
            magnitude has a common declared physical scale across files.
        metric_semantics: Permitted interpretation of descriptive metrics.
        interpretation_note: Concise human-readable boundary.
    """

    value_status: str
    edf_original_unit_status: str
    edf_original_unit_normalized: str
    in_memory_value_unit: str
    values_are_calibrated_physical_units: bool
    values_are_mne_scaled: bool
    absolute_cross_file_magnitude_interpretable: bool
    metric_semantics: str
    interpretation_note: str

    def as_dict(self) -> dict[str, object]:
        """Return the classification as a plain output-row dictionary."""

        return asdict(self)


def _read_ascii_fields(
    signal_header: bytes,
    offset: int,
    width: int,
    count: int,
) -> tuple[list[str], int]:
    """Decode one repeated fixed-width EDF signal-header field.

    Args:
        signal_header: Concatenated per-signal EDF header bytes.
        offset: Starting byte offset within ``signal_header``.
        width: Field width per signal.
        count: Number of signals.

    Returns:
        A list of stripped Latin-1 strings and the next byte offset.

    Side effects:
        None.
    """

    values = []
    for index in range(count):
        start = offset + index * width
        stop = start + width
        values.append(signal_header[start:stop].decode("latin-1").strip())
    return values, offset + width * count


def _required_number(text: str, field: str, label: str, number_type: type):
    """Parse a required numeric EDF field with a useful error message."""

    try:
        return number_type(text)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"EDF signal {label!r} has invalid {field} value {text!r}"
        ) from error


def _read_fixed_header(stream: BinaryIO, path: Path) -> tuple[bytes, int]:
    """Read and validate the 256-byte EDF fixed header."""

    fixed_header = stream.read(256)
    if len(fixed_header) != 256:
        raise ValueError(f"EDF header is truncated: {path}")
    signal_count_text = fixed_header[252:256].decode("latin-1").strip()
    signal_count = _required_number(
        signal_count_text,
        "number of signals",
        path.name,
        int,
    )
    if signal_count <= 0:
        raise ValueError(f"EDF has no signals: {path}")
    expected_header_bytes = 256 + signal_count * 256
    declared_header_text = fixed_header[184:192].decode("latin-1").strip()
    declared_header_bytes = _required_number(
        declared_header_text,
        "header byte count",
        path.name,
        int,
    )
    if declared_header_bytes != expected_header_bytes:
        raise ValueError(
            f"EDF header byte count is inconsistent for {path.name}: "
            f"declared {declared_header_bytes}, expected {expected_header_bytes}"
        )
    return fixed_header, signal_count


def read_edf_signal_headers(path: Path) -> list[EdfSignalHeader]:
    """Read original per-signal metadata directly from an EDF header.

    Args:
        path: Raw EDF path.

    Returns:
        One :class:`EdfSignalHeader` per EDF signal, including an EDF
        annotation signal when present.

    Side effects:
        Opens ``path`` in binary read-only mode and reads only header bytes.
    """

    with path.open("rb") as stream:
        _, signal_count = _read_fixed_header(stream, path)
        signal_header = stream.read(signal_count * 256)
    if len(signal_header) != signal_count * 256:
        raise ValueError(f"EDF signal header is truncated: {path}")

    offset = 0
    labels, offset = _read_ascii_fields(signal_header, offset, 16, signal_count)
    transducers, offset = _read_ascii_fields(signal_header, offset, 80, signal_count)
    dimensions, offset = _read_ascii_fields(signal_header, offset, 8, signal_count)
    physical_minima, offset = _read_ascii_fields(signal_header, offset, 8, signal_count)
    physical_maxima, offset = _read_ascii_fields(signal_header, offset, 8, signal_count)
    digital_minima, offset = _read_ascii_fields(signal_header, offset, 8, signal_count)
    digital_maxima, offset = _read_ascii_fields(signal_header, offset, 8, signal_count)
    prefiltering, offset = _read_ascii_fields(signal_header, offset, 80, signal_count)
    samples_per_record, offset = _read_ascii_fields(signal_header, offset, 8, signal_count)
    reserved, offset = _read_ascii_fields(signal_header, offset, 32, signal_count)
    if offset != len(signal_header):
        raise AssertionError("internal EDF signal-header parser offset mismatch")

    rows: list[EdfSignalHeader] = []
    for index, label in enumerate(labels):
        rows.append(
            EdfSignalHeader(
                signal_index=index,
                label=label,
                transducer_type=transducers[index],
                physical_dimension=dimensions[index],
                physical_minimum=_required_number(
                    physical_minima[index], "physical minimum", label, float
                ),
                physical_maximum=_required_number(
                    physical_maxima[index], "physical maximum", label, float
                ),
                digital_minimum=_required_number(
                    digital_minima[index], "digital minimum", label, int
                ),
                digital_maximum=_required_number(
                    digital_maxima[index], "digital maximum", label, int
                ),
                prefiltering=prefiltering[index],
                samples_per_record=_required_number(
                    samples_per_record[index], "samples per record", label, int
                ),
                reserved=reserved[index],
            )
        )
    return rows


def normalize_edf_unit(value: object) -> str:
    """Normalize common EDF voltage symbols without inventing a dimension."""

    text = str(value or "").strip()
    if text.upper() in UNKNOWN_EDF_UNITS:
        return ""
    normalized = text.replace("µ", "u").replace("μ", "u").replace(" ", "")
    if normalized.upper() in {"V", "MV", "UV", "NV"}:
        return normalized.upper()
    return normalized


def classify_unit_status(
    configured_channel_type: str,
    edf_original_physical_dimension: object,
) -> UnitStatus:
    """Classify the scientific unit semantics for one channel.

    Args:
        configured_channel_type: Reviewed type: ``eeg``, ``eog``, ``emg``,
            ``stim``, or another explicit MNE type.
        edf_original_physical_dimension: Original EDF dimension text. Blank,
            ``n/a``, and equivalent values remain unknown.

    Returns:
        A :class:`UnitStatus` that controls output labels and comparison
        boundaries.

    Side effects:
        None.
    """

    channel_type = configured_channel_type.strip().lower()
    original_unit = normalize_edf_unit(edf_original_physical_dimension)
    known_unit = bool(original_unit)

    # Stim samples encode event/digital states.  A physical-looking MNE unit or
    # EDF range must never turn them into a continuous amplitude measurement.
    if channel_type == "stim":
        return UnitStatus(
            value_status="digital_stim_values",
            edf_original_unit_status="known_but_not_used_for_stim" if known_unit else "unknown",
            edf_original_unit_normalized=original_unit,
            in_memory_value_unit="digital_code",
            values_are_calibrated_physical_units=False,
            values_are_mne_scaled=False,
            absolute_cross_file_magnitude_interpretable=False,
            metric_semantics="event_transitions_and_digital_states_only",
            interpretation_note=(
                "Trigger values are digital event states; continuous amplitude and PSD "
                "interpretation is not permitted."
            ),
        )

    if not known_unit:
        return UnitStatus(
            value_status="arbitrary_unknown_acquisition_units",
            edf_original_unit_status="unknown",
            edf_original_unit_normalized="",
            in_memory_value_unit="unknown_acquisition_unit",
            values_are_calibrated_physical_units=False,
            values_are_mne_scaled=False,
            absolute_cross_file_magnitude_interpretable=False,
            metric_semantics="within_channel_file_morphology_frequency_and_flatness_only",
            interpretation_note=(
                "The EDF physical dimension is unknown. Morphology, frequency shape, "
                "flatness, rail, and within-channel/file ratios remain descriptive; "
                "absolute magnitudes must not be compared across files."
            ),
        )

    if original_unit in VOLTAGE_UNIT_SCALE_TO_VOLTS:
        return UnitStatus(
            value_status="mne_scaled_calibrated_voltage",
            edf_original_unit_status="known_voltage",
            edf_original_unit_normalized=original_unit,
            in_memory_value_unit="V",
            values_are_calibrated_physical_units=True,
            values_are_mne_scaled=original_unit != "V",
            absolute_cross_file_magnitude_interpretable=True,
            metric_semantics="calibrated_physical_voltage_qc",
            interpretation_note=(
                "The EDF declares a voltage dimension and digital-to-physical range; "
                "MNE returns volts. Cross-file magnitude QC is technically supported, "
                "subject to the common acquisition context."
            ),
        )

    return UnitStatus(
        value_status="calibrated_physical_units",
        edf_original_unit_status="known_non_voltage",
        edf_original_unit_normalized=original_unit,
        in_memory_value_unit=original_unit,
        values_are_calibrated_physical_units=True,
        values_are_mne_scaled=False,
        absolute_cross_file_magnitude_interpretable=True,
        metric_semantics="calibrated_physical_unit_qc",
        interpretation_note=(
            "The EDF declares a physical dimension and digital-to-physical range. "
            "MNE scaling for this uncommon unit must be verified before mutation."
        ),
    )


def mne_scaling_for_unit_status(status: UnitStatus) -> float:
    """Return the expected EDF-physical to MNE-memory scale for a unit status.

    Unknown auxiliary values deliberately retain scale 1.0: this describes
    observed MNE behavior without assigning those numbers a physical unit.
    Stim values also retain their digital numeric scale.
    """

    if status.value_status == "mne_scaled_calibrated_voltage":
        return VOLTAGE_UNIT_SCALE_TO_VOLTS[status.edf_original_unit_normalized]
    return 1.0
