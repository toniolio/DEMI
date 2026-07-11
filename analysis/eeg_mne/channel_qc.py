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

import numpy as np
from scipy.signal import welch


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


def deterministic_chunk_bounds(
    n_times: int,
    chunk_samples: int,
) -> list[tuple[int, int]]:
    """Partition a signal into deterministic, contiguous full-coverage chunks.

    Args:
        n_times: Total sample count.
        chunk_samples: Maximum samples per chunk.

    Returns:
        Ordered ``(start, stop)`` bounds that cover each sample exactly once.

    Side effects:
        None.
    """

    if n_times < 0:
        raise ValueError("n_times must be non-negative")
    if chunk_samples <= 0:
        raise ValueError("chunk_samples must be positive")
    return [
        (start, min(n_times, start + chunk_samples))
        for start in range(0, n_times, chunk_samples)
    ]


def deterministic_sample_stride(n_times: int, maximum_samples: int) -> int:
    """Return a stable stride that samples across the complete duration.

    Args:
        n_times: Total sample count.
        maximum_samples: Desired upper bound for deterministic quantile
            samples.

    Returns:
        Integer stride of at least one.

    Side effects:
        None.
    """

    if n_times < 0:
        raise ValueError("n_times must be non-negative")
    if maximum_samples <= 0:
        raise ValueError("maximum_samples must be positive")
    return max(1, int(np.ceil(n_times / maximum_samples)))


def longest_true_run(values: np.ndarray) -> int:
    """Return the longest consecutive run of true values in one dimension."""

    array = np.asarray(values, dtype=bool)
    if array.ndim != 1:
        raise ValueError("longest_true_run expects a one-dimensional array")
    if array.size == 0 or not np.any(array):
        return 0
    padded = np.concatenate(([False], array, [False]))
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    return int(np.max(changes[1::2] - changes[::2]))


def _finite_float(value: float) -> float | None:
    """Convert non-finite numerical output to ``None`` for table safety."""

    return float(value) if np.isfinite(value) else None


def summarize_signal_metrics(
    signal: np.ndarray,
    sfreq: float,
    *,
    digital_step_in_memory_units: float | None = None,
    rail_minimum_in_memory_units: float | None = None,
    rail_maximum_in_memory_units: float | None = None,
) -> dict[str, float | int | None]:
    """Compute deterministic descriptive metrics for one synthetic/full signal.

    This helper is intentionally descriptive. It does not decide whether a
    channel is bad. The full driver uses the same definitions while scanning
    EDFs in bounded-memory chunks.

    Args:
        signal: One-dimensional numeric samples.
        sfreq: Sampling frequency in Hz.
        digital_step_in_memory_units: Optional value represented by one EDF
            digital-code increment after MNE scaling.
        rail_minimum_in_memory_units: Optional calibrated lower EDF rail.
        rail_maximum_in_memory_units: Optional calibrated upper EDF rail.

    Returns:
        Dictionary containing finite/dropout, location/scale, range, repeated
        sample, near-flat run, rail, and extreme-step metrics.

    Side effects:
        None.
    """

    values = np.asarray(signal, dtype=float)
    if values.ndim != 1:
        raise ValueError("summarize_signal_metrics expects a one-dimensional signal")
    if sfreq <= 0:
        raise ValueError("sfreq must be positive")

    n_samples = int(values.size)
    finite_mask = np.isfinite(values)
    finite = values[finite_mask]
    nonfinite_count = int(n_samples - finite.size)
    if finite.size == 0:
        return {
            "n_samples": n_samples,
            "nonfinite_count": nonfinite_count,
            "finite_proportion": 0.0 if n_samples else None,
            "median": None,
            "robust_mad": None,
            "robust_scale_mad": None,
            "robust_q01": None,
            "robust_q99": None,
            "robust_peak_to_peak_q01_q99": None,
            "minimum": None,
            "maximum": None,
            "peak_to_peak": None,
            "variance": None,
            "unchanged_difference_count": 0,
            "unchanged_difference_proportion": None,
            "near_flat_tolerance": None,
            "longest_near_flat_run_samples": 0,
            "longest_near_flat_run_seconds": 0.0,
            "rail_sample_count": 0,
            "rail_sample_proportion": None,
            "maximum_absolute_step": None,
            "median_absolute_step": None,
            "p999_absolute_step": None,
            "extreme_step_ratio": None,
        }

    median = float(np.median(finite))
    mad = float(np.median(np.abs(finite - median)))
    q01, q99 = np.quantile(finite, [0.01, 0.99])
    scale_reference = max(abs(median), mad, float(np.max(np.abs(finite))), 1.0)
    quantization_tolerance = (
        abs(float(digital_step_in_memory_units)) * 0.5
        if digital_step_in_memory_units is not None
        and np.isfinite(digital_step_in_memory_units)
        else 0.0
    )
    near_flat_tolerance = max(
        quantization_tolerance,
        np.finfo(float).eps * scale_reference * 8.0,
    )

    if n_samples >= 2:
        adjacent_finite = finite_mask[:-1] & finite_mask[1:]
        differences = np.diff(values)
        valid_differences = differences[adjacent_finite]
        absolute_steps = np.abs(valid_differences)
        unchanged_count = int(np.count_nonzero(valid_differences == 0.0))
        unchanged_proportion = (
            float(unchanged_count / valid_differences.size)
            if valid_differences.size
            else None
        )
        near_flat = np.zeros(n_samples - 1, dtype=bool)
        near_flat[adjacent_finite] = absolute_steps <= near_flat_tolerance
        longest_flat_samples = longest_true_run(near_flat) + 1 if np.any(near_flat) else 1
        maximum_step = float(np.max(absolute_steps)) if absolute_steps.size else None
        median_step = float(np.median(absolute_steps)) if absolute_steps.size else None
        p999_step = float(np.quantile(absolute_steps, 0.999)) if absolute_steps.size else None
    else:
        unchanged_count = 0
        unchanged_proportion = None
        longest_flat_samples = n_samples
        maximum_step = None
        median_step = None
        p999_step = None

    rail_count = 0
    rail_metadata = (
        rail_minimum_in_memory_units is not None
        and rail_maximum_in_memory_units is not None
        and np.isfinite(rail_minimum_in_memory_units)
        and np.isfinite(rail_maximum_in_memory_units)
    )
    if rail_metadata:
        rail_count = int(
            np.count_nonzero(
                (np.abs(finite - float(rail_minimum_in_memory_units)) <= near_flat_tolerance)
                | (np.abs(finite - float(rail_maximum_in_memory_units)) <= near_flat_tolerance)
            )
        )

    step_denominator = max(
        median_step or 0.0,
        abs(float(digital_step_in_memory_units or 0.0)),
        np.finfo(float).eps * scale_reference,
    )
    extreme_step_ratio = (
        float(maximum_step / step_denominator) if maximum_step is not None else None
    )

    return {
        "n_samples": n_samples,
        "nonfinite_count": nonfinite_count,
        "finite_proportion": float(finite.size / n_samples) if n_samples else None,
        "median": median,
        "robust_mad": mad,
        "robust_scale_mad": float(1.4826 * mad),
        "robust_q01": float(q01),
        "robust_q99": float(q99),
        "robust_peak_to_peak_q01_q99": float(q99 - q01),
        "minimum": float(np.min(finite)),
        "maximum": float(np.max(finite)),
        "peak_to_peak": float(np.ptp(finite)),
        "variance": float(np.var(finite)),
        "unchanged_difference_count": unchanged_count,
        "unchanged_difference_proportion": unchanged_proportion,
        "near_flat_tolerance": float(near_flat_tolerance),
        "longest_near_flat_run_samples": int(longest_flat_samples),
        "longest_near_flat_run_seconds": float(longest_flat_samples / sfreq),
        "rail_sample_count": rail_count,
        "rail_sample_proportion": float(rail_count / finite.size) if rail_metadata else None,
        "maximum_absolute_step": maximum_step,
        "median_absolute_step": median_step,
        "p999_absolute_step": p999_step,
        "extreme_step_ratio": extreme_step_ratio,
    }


def _mean_psd_in_band(
    frequencies: np.ndarray,
    psd: np.ndarray,
    low_hz: float,
    high_hz: float,
) -> float | None:
    """Return mean PSD in an inclusive band or ``None`` when unavailable."""

    mask = (frequencies >= low_hz) & (frequencies <= high_hz)
    if not np.any(mask):
        return None
    finite = psd[mask][np.isfinite(psd[mask])]
    return float(np.mean(finite)) if finite.size else None


def summarize_spectral_metrics(
    signal: np.ndarray,
    sfreq: float,
    *,
    nperseg_seconds: float = 4.0,
) -> dict[str, float | int | None]:
    """Compute scale-aware PSD summaries and scale-free frequency ratios.

    Args:
        signal: One-dimensional continuous physical/auxiliary signal. Stim
            channels must be skipped by the caller.
        sfreq: Sampling frequency in Hz.
        nperseg_seconds: Welch segment duration.

    Returns:
        Welch configuration, broadband/band summaries, a 60-Hz line-to-
        sideband ratio, and a low/high-frequency ratio. Absolute PSD values
        inherit the channel's unit status; ratios remain useful for unknown-
        unit within-channel/file QC.

    Side effects:
        None.
    """

    values = np.asarray(signal, dtype=float)
    if values.ndim != 1:
        raise ValueError("summarize_spectral_metrics expects a one-dimensional signal")
    if sfreq <= 0 or nperseg_seconds <= 0:
        raise ValueError("sfreq and nperseg_seconds must be positive")
    if values.size < 2 or not np.all(np.isfinite(values)):
        return {
            "welch_nperseg": None,
            "psd_frequency_resolution_hz": None,
            "median_psd_1_100_hz": None,
            "mean_psd_1_4_hz": None,
            "mean_psd_30_45_hz": None,
            "mean_psd_59_61_hz": None,
            "mean_psd_60hz_sidebands": None,
            "line_noise_ratio_60hz": None,
            "low_high_frequency_ratio": None,
        }

    nperseg = max(2, min(values.size, int(round(nperseg_seconds * sfreq))))
    frequencies, psd = welch(
        values,
        fs=sfreq,
        nperseg=nperseg,
        noverlap=nperseg // 2,
        detrend="constant",
        scaling="density",
    )
    broadband_mask = (frequencies >= 1.0) & (frequencies <= min(100.0, sfreq / 2.0))
    broadband = psd[broadband_mask]
    low = _mean_psd_in_band(frequencies, psd, 1.0, 4.0)
    high = _mean_psd_in_band(frequencies, psd, 30.0, 45.0)
    line = _mean_psd_in_band(frequencies, psd, 59.0, 61.0)
    side_mask = (
        ((frequencies >= 55.0) & (frequencies <= 58.0))
        | ((frequencies >= 62.0) & (frequencies <= 65.0))
    )
    side_values = psd[side_mask]
    side = float(np.mean(side_values)) if side_values.size else None
    line_ratio = (
        float(line / side)
        if line is not None and side is not None and side > 0
        else None
    )
    low_high_ratio = (
        float(low / high)
        if low is not None and high is not None and high > 0
        else None
    )
    frequency_resolution = frequencies[1] - frequencies[0] if frequencies.size > 1 else np.nan
    return {
        "welch_nperseg": nperseg,
        "psd_frequency_resolution_hz": _finite_float(float(frequency_resolution)),
        "median_psd_1_100_hz": (
            float(np.median(broadband[np.isfinite(broadband)]))
            if np.any(np.isfinite(broadband))
            else None
        ),
        "mean_psd_1_4_hz": low,
        "mean_psd_30_45_hz": high,
        "mean_psd_59_61_hz": line,
        "mean_psd_60hz_sidebands": side,
        "line_noise_ratio_60hz": line_ratio,
        "low_high_frequency_ratio": low_high_ratio,
    }
