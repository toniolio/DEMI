"""Signal stages through reference, interpolation, and validation.

This module contains the production implementations of the fixed line-noise
branch, analysis/ICA filters, criterion-level PyPREP detection, accepted
published-criterion union, explicit reference calculation/application, and
MNE spherical interpolation. Each operation returns manifest evidence instead
of hiding multiple scientific decisions inside one opaque call.

Inputs are typed, montaged, full-session MNE ``Raw`` objects. Signal mutation
is confined to caller-owned in-memory copies. The module writes no files,
changes no EDF, constructs no epochs, runs no AutoReject, and computes no CSD.
"""

from __future__ import annotations

import gc
from typing import Any, Mapping, Sequence

import mne
import numpy as np

from .contracts import (
    EEG_TARGET_CHANNELS,
    MASTOID_CHANNELS,
    PRIMARY_PYPREP_CRITERIA,
    REPORT_ONLY_PYPREP_CRITERIA,
    SCALP_SOURCE_CHANNELS,
    sha256_bytes,
)


class ObjectiveStop(RuntimeError):
    """A scientifically accepted per-recording stop that is not a run failure."""

    def __init__(self, stage: str, code: str, detail: str, evidence: Mapping[str, Any] | None = None):
        """Create an objective stop with stable code and optional evidence."""

        super().__init__(detail)
        self.stage = stage
        self.code = code
        self.detail = detail
        self.evidence = dict(evidence or {})


def _picks_for_types(raw: mne.io.BaseRaw, channel_types: Sequence[str]) -> list[str]:
    """Return channel names matching explicit MNE channel types."""

    allowed = set(channel_types)
    return [
        name
        for name, channel_type in zip(raw.ch_names, raw.get_channel_types())
        if channel_type in allowed
    ]


def apply_line_noise_branch(raw: mne.io.BaseRaw, config: Mapping[str, Any]) -> dict[str, Any]:
    """Apply the fixed run-level 60-Hz spectrum-fit branch when enabled.

    Args:
        raw: Caller-owned preloaded branch.
        config: Validated complete production configuration.

    Returns:
        Exact toggle, method, frequency, picks, and edge-handling evidence.

    Side effects:
        When enabled, changes signal samples in ``raw`` in memory.
    """

    line = config["line_noise"]
    picks = _picks_for_types(raw, line["picks"])
    evidence = {
        "enabled": bool(line["enabled"]),
        "frequency_hz": float(line["frequency_hz"]),
        "method": line["method"],
        "filter_length": line["filter_length"],
        "multitaper_bandwidth_hz": float(line["multitaper_bandwidth_hz"]),
        "p_value": float(line["p_value"]),
        "picks": picks,
        "fixed_for_run": True,
        "participant_specific_tuning": False,
    }
    if line["enabled"]:
        raw.notch_filter(
            freqs=[float(line["frequency_hz"])],
            picks=picks,
            method=line["method"],
            filter_length=line["filter_length"],
            mt_bandwidth=float(line["multitaper_bandwidth_hz"]),
            p_value=float(line["p_value"]),
            verbose="ERROR",
        )
    return evidence


def apply_fir_filter(
    raw: mne.io.BaseRaw,
    *,
    high_pass_hz: float,
    low_pass_hz: float,
    settings: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply and describe one zero-phase FIR filter at the native sample rate.

    Args:
        raw: Caller-owned preloaded branch.
        high_pass_hz: High-pass cutoff.
        low_pass_hz: Low-pass cutoff.
        settings: Filter method, design, padding, skip annotations, and picks.

    Returns:
        Exact design and edge-handling metadata, including realized FIR length.

    Side effects:
        Filters selected channels in ``raw`` in memory.
    """

    sfreq = float(raw.info["sfreq"])
    picks = _picks_for_types(raw, settings["picks"])
    taps = mne.filter.create_filter(
        None,
        sfreq,
        l_freq=high_pass_hz,
        h_freq=low_pass_hz,
        method=settings["method"],
        phase=settings["phase"],
        fir_design=settings["fir_design"],
        verbose="ERROR",
    )
    raw.filter(
        high_pass_hz,
        low_pass_hz,
        picks=picks,
        method=settings["method"],
        phase=settings["phase"],
        fir_design=settings["fir_design"],
        pad=settings["pad"],
        skip_by_annotation=tuple(settings["skip_by_annotation"]),
        verbose="ERROR",
    )
    return {
        "high_pass_hz": float(high_pass_hz),
        "low_pass_hz": float(low_pass_hz),
        "method": settings["method"],
        "phase": settings["phase"],
        "fir_design": settings["fir_design"],
        "sampling_frequency_hz": sfreq,
        "filter_length_samples": int(len(taps)),
        "filter_length_seconds": float(len(taps) / sfreq),
        "picks": picks,
        "edge_handling": {
            "pad": settings["pad"],
            "skip_by_annotation": list(settings["skip_by_annotation"]),
            "manifest_only_edge_metadata": True,
            "filter_edge_annotation_added": False,
            "half_filter_length_seconds": float(len(taps) / sfreq / 2.0),
        },
    }


def prepare_detector_input(raw: mne.io.BaseRaw, config: Mapping[str, Any]) -> tuple[mne.io.BaseRaw, dict[str, Any]]:
    """Create the accepted unfiltered 30-scalp full-session detector input.

    Args:
        raw: Typed and montaged source branch before line/filter operations.
        config: Validated production configuration.

    Returns:
        Independent EEG-only detector ``Raw`` and surface evidence.

    Side effects:
        Allocates an in-memory copy. It writes nothing.
    """

    channels = list(config["channels"]["detector_source_channels"])
    if channels != list(SCALP_SOURCE_CHANNELS):
        raise ValueError("Detector source channels differ from the accepted 30-scalp contract.")
    detector_raw = raw.copy().pick(channels)
    if detector_raw.ch_names != channels:
        raise ValueError("Detector input order differs from the accepted source pool.")
    if detector_raw.get_channel_types(unique=True) != ["eeg"]:
        raise ValueError("Detector input must contain only EEG channels.")
    if detector_raw.n_times != raw.n_times or detector_raw.info["sfreq"] != raw.info["sfreq"]:
        raise ValueError("Detector input must retain the full native continuous session.")
    return detector_raw, {
        "input_branch": config["detector"]["input_branch"],
        "source_channels": channels,
        "source_channel_count": len(channels),
        "m1_m2_source_excluded": list(MASTOID_CHANNELS),
        "full_session_sample_count": int(detector_raw.n_times),
        "sampling_frequency_hz": float(detector_raw.info["sfreq"]),
        "resampled": False,
        "referenced": False,
        "notched": False,
        "analysis_filtered": False,
        "montage": "standard_1005",
    }


def _criterion_scores(noisy: Any, channel_names: Sequence[str]) -> tuple[list[dict[str, Any]], list[str]]:
    """Extract available PyPREP diagnostics and surface nonfinite numeric values."""

    score_keys = {
        "bad_by_deviation": ("robust_channel_deviations",),
        "bad_by_hf_noise": ("hf_noise_zscores",),
        "bad_by_correlation": ("bad_window_fractions", "median_max_correlations"),
        "bad_by_dropout": ("bad_window_fractions",),
        "bad_by_psd": ("psd_zscore", "zscore_low", "zscore_mid", "zscore_high"),
        "bad_by_ransac": ("bad_window_fractions",),
    }
    rows = [{"channel": channel} for channel in channel_names]
    nonfinite: list[str] = []
    for criterion, keys in score_keys.items():
        info = noisy._extra_info.get(criterion, {})  # PyPREP has no public score API.
        for key in keys:
            values = info.get(key)
            if values is None:
                continue
            array = np.asarray(values)
            if array.ndim != 1 or array.size != len(channel_names):
                continue
            for index, value in enumerate(array.astype(float)):
                rows[index][f"{criterion}_{key}"] = float(value)
                if criterion not in REPORT_ONLY_PYPREP_CRITERIA and not np.isfinite(value):
                    nonfinite.append(f"{criterion}.{key}.{channel_names[index]}")
    return rows, sorted(nonfinite)


def _run_noisy_channels_once(raw: mne.io.BaseRaw, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Run one seeded PyPREP detector repetition and return criterion evidence."""

    from pyprep.find_noisy_channels import NoisyChannels

    noisy = NoisyChannels(
        raw.copy(),
        do_detrend=bool(settings["do_detrend"]),
        random_state=int(settings["random_seed"]),
        matlab_strict=bool(settings["matlab_strict"]),
        ransac=bool(settings["ransac"]),
        correlation=bool(settings["correlation"]),
        reject_by_annotation=settings.get("reject_by_annotation"),
    )
    noisy.find_all_bads(
        ransac=bool(settings["ransac"]),
        channel_wise=bool(settings["channel_wise"]),
        max_chunk_size=settings.get("max_chunk_size"),
        correlation=bool(settings["correlation"]),
    )
    raw_criteria = noisy.get_bads(as_dict=True)
    keys = tuple(PRIMARY_PYPREP_CRITERIA) + tuple(REPORT_ONLY_PYPREP_CRITERIA)
    criteria = {key: sorted(set(raw_criteria.get(key, []))) for key in keys}
    rows, nonfinite = _criterion_scores(noisy, raw.ch_names)
    for row in rows:
        channel = row["channel"]
        for criterion in keys:
            row[criterion] = channel in criteria[criterion]
    del noisy
    gc.collect()
    return {"criteria": criteria, "channel_details": rows, "nonfinite_required_metrics": nonfinite}


def run_pyprep_detection(raw: mne.io.BaseRaw, config: Mapping[str, Any]) -> dict[str, Any]:
    """Run repeated fixed PyPREP detection and enforce exact determinism.

    Args:
        raw: Accepted 30-channel detector input.
        config: Validated complete configuration.

    Returns:
        Criterion lists, channel diagnostics, and repeated-result hashes.

    Raises:
        ObjectiveStop: If required metrics are nonfinite or repeated criterion
            lists are not identical.

    Side effects:
        Allocates independent PyPREP copies. It writes nothing.
    """

    settings = config["detector"]
    repeats = [_run_noisy_channels_once(raw, settings) for _ in range(settings["determinism_repeats"])]
    criterion_payloads = [repeat["criteria"] for repeat in repeats]
    hashes = [sha256_bytes(str(payload).encode("utf-8")) for payload in criterion_payloads]
    if any(payload != criterion_payloads[0] for payload in criterion_payloads[1:]):
        raise ObjectiveStop(
            "pyprep_criterion_detection",
            "nondeterministic_repeated_detector_result",
            "Seeded PyPREP criterion lists differed across repeated runs.",
            {"repeat_hashes": hashes, "repeat_criteria": criterion_payloads},
        )
    nonfinite = sorted(set().union(*(set(repeat["nonfinite_required_metrics"]) for repeat in repeats)))
    if nonfinite:
        raise ObjectiveStop(
            "pyprep_criterion_detection",
            "nonfinite_required_detector_metrics",
            f"PyPREP returned nonfinite required detector metrics: {nonfinite}",
            {"nonfinite_required_metrics": nonfinite},
        )
    return {
        "method": "pyprep.NoisyChannels",
        "parameters": {
            key: settings[key]
            for key in (
                "random_seed",
                "do_detrend",
                "correlation",
                "ransac",
                "matlab_strict",
                "channel_wise",
                "max_chunk_size",
                "reject_by_annotation",
            )
        },
        "criteria": criterion_payloads[0],
        "channel_details": repeats[0]["channel_details"],
        "determinism_repeats": int(settings["determinism_repeats"]),
        "repeat_result_hashes": hashes,
        "repeated_results_identical": True,
        "nonfinite_required_metrics": [],
    }


def accepted_global_bad_union(detection: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    """Construct only the accepted published-criterion global-bad union.

    Args:
        detection: Output of :func:`run_pyprep_detection`.
        config: Validated complete configuration.

    Returns:
        Primary bads, PSD report-only findings, and threshold evidence.

    Raises:
        ObjectiveStop: If the detector returns a non-source label. A primary
            count above the historical 25% boundary is warning evidence and
            does not stop ordinary preprocessing.

    Side effects:
        None.
    """

    criteria = detection["criteria"]
    primary = sorted(
        set().union(*(set(criteria[key]) for key in config["detector"]["primary_criteria"]))
    )
    invalid = sorted(set(primary) - set(SCALP_SOURCE_CHANNELS))
    if invalid:
        raise ObjectiveStop(
            "accepted_global_bad_union",
            "invalid_detector_channel_label",
            f"Primary detector returned labels outside the 30-channel pool: {invalid}",
        )
    psd = sorted(set().union(*(set(criteria[key]) for key in config["detector"]["report_only_criteria"])))
    psd_only = sorted(set(psd) - set(primary))
    denominator = int(config["interpolation"]["proportion_denominator"])
    proportion = len(primary) / denominator
    evidence = {
        "authoritative_criteria": list(PRIMARY_PYPREP_CRITERIA),
        "accepted_global_bads": primary,
        "accepted_global_bad_count": len(primary),
        "accepted_global_bad_proportion": proportion,
        "proportion_denominator": denominator,
        "psd_report_only_findings": psd,
        "psd_exclusive_report_only_findings": psd_only,
        "psd_entered_authoritative_union": False,
        "m1_m2_detector_source_excluded": list(MASTOID_CHANNELS),
    }
    warning_threshold = float(config["interpolation"]["warning_proportion"])
    warning = proportion > warning_threshold
    evidence["high_interpolation_warning"] = {
        "triggered": warning,
        "code": "global_bad_proportion_above_25_percent",
        "threshold_proportion": warning_threshold,
        "comparison": "strictly_greater_than",
        "observed_count": len(primary),
        "observed_proportion": proportion,
        "denominator": denominator,
        "action": "continue_preprocessing_and_mark_complete_with_qc_warning",
        "participant_event_epoch_or_analytic_decision": False,
    }
    return evidence


def calculate_reference_sources(global_bads: Sequence[str], config: Mapping[str, Any]) -> dict[str, Any]:
    """Select non-bad members of the fixed scalp average-reference pool."""

    bad_set = set(global_bads)
    sources = [channel for channel in SCALP_SOURCE_CHANNELS if channel not in bad_set]
    minimum = int(config["reference"]["minimum_usable_scalp_channels"])
    evidence = {
        "reference_kind": config["reference"]["kind"],
        "reference_sources": sources,
        "reference_source_count": len(sources),
        "excluded_accepted_global_bads": sorted(bad_set),
        "m1_m2_source_excluded": list(MASTOID_CHANNELS),
        "reference_targets": list(EEG_TARGET_CHANNELS),
        "reference_target_count": len(EEG_TARGET_CHANNELS),
        "m1_m2_target_retained": list(MASTOID_CHANNELS),
        "minimum_usable_scalp_channels": minimum,
    }
    if len(sources) < minimum:
        raise ObjectiveStop(
            "reference_estimation",
            "too_few_usable_scalp_channels_for_reference",
            f"Only {len(sources)} accepted scalp reference sources remain; minimum is {minimum}.",
            evidence,
        )
    return evidence


def apply_reference(raw: mne.io.BaseRaw, reference: Mapping[str, Any]) -> dict[str, Any]:
    """Apply one explicit scalp average to every retained EEG target.

    Args:
        raw: Analysis or ICA-fitting branch.
        reference: Output of :func:`calculate_reference_sources`.

    Returns:
        Reference signal hash and target verification evidence.

    Side effects:
        Subtracts the explicit reference signal from all EEG channels in
        ``raw``. EOG, EMG, and stim channels are unchanged.
    """

    sources = list(reference["reference_sources"])
    targets = [
        name
        for name, channel_type in zip(raw.ch_names, raw.get_channel_types())
        if channel_type == "eeg"
    ]
    if targets != list(EEG_TARGET_CHANNELS):
        raise ValueError(f"Reference targets differ from accepted EEG targets: {targets}")
    _, reference_data = mne.set_eeg_reference(
        raw,
        ref_channels=sources,
        copy=False,
        projection=False,
        ch_type="eeg",
        verbose="ERROR",
    )
    if reference_data is None or reference_data.shape != (raw.n_times,):
        raise RuntimeError("MNE did not return the expected explicit reference signal.")
    if not np.isfinite(reference_data).all():
        raise ObjectiveStop(
            "reference_application",
            "nonfinite_post_reference_signal",
            "The calculated scalp reference contains nonfinite samples.",
        )
    return {
        **dict(reference),
        "applied": True,
        "reference_signal_sha256_float64": sha256_bytes(
            np.ascontiguousarray(reference_data, dtype=np.float64).tobytes()
        ),
        "mne_custom_ref_applied": int(raw.info["custom_ref_applied"]),
        "all_32_eeg_targets_received_reference": targets == list(EEG_TARGET_CHANNELS),
    }


def interpolate_global_bads(raw: mne.io.BaseRaw, global_bads: Sequence[str], config: Mapping[str, Any]) -> dict[str, Any]:
    """Interpolate accepted bad scalp channels with MNE spherical splines.

    Args:
        raw: Referenced analysis or ICA-fitting branch.
        global_bads: Accepted primary scalp bads.
        config: Validated complete configuration.

    Returns:
        Candidate, success, proportion, and operational bad-reset evidence.

    Side effects:
        Replaces bad-channel samples in ``raw`` and resets ``raw.info['bads']``
        after successful interpolation.
    """

    candidates = sorted(global_bads)
    if set(candidates) & set(MASTOID_CHANNELS):
        raise ValueError("M1/M2 cannot be primary interpolation candidates.")
    if not set(candidates).issubset(SCALP_SOURCE_CHANNELS):
        raise ValueError("Interpolation candidates must belong to the 30-channel scalp pool.")
    raw.info["bads"] = candidates
    if candidates:
        raw.interpolate_bads(
            reset_bads=bool(config["interpolation"]["reset_bads_after_success"]),
            mode="accurate",
            origin="auto",
            method=dict(eeg="spline"),
            on_bad_position="raise",
            verbose="ERROR",
        )
    elif config["interpolation"]["reset_bads_after_success"]:
        raw.info["bads"] = []
    if raw.info["bads"]:
        raise ObjectiveStop(
            "interpolation",
            "operational_bad_list_not_reset",
            f"Operational bad list remained after interpolation: {raw.info['bads']}",
        )
    return {
        "detected_global_bads": candidates,
        "interpolation_candidates": candidates,
        "successfully_interpolated_channels": candidates,
        "interpolation_count": len(candidates),
        "interpolation_proportion": len(candidates) / len(SCALP_SOURCE_CHANNELS),
        "proportion_denominator": len(SCALP_SOURCE_CHANNELS),
        "method": "MNE spherical spline interpolation",
        "montage": "standard_1005",
        "m1_m2_primary_interpolation_candidates": False,
        "reset_bads_after_success": True,
        "operational_bad_list_after": [],
    }


def validate_finite_continuous(raw: mne.io.BaseRaw, chunk_seconds: float) -> dict[str, Any]:
    """Scan a full continuous branch for nonfinite samples in bounded chunks.

    Args:
        raw: Preloaded or file-backed continuous MNE object.
        chunk_seconds: Maximum scan duration per chunk.

    Returns:
        Coverage and channel/montage/bad-list validation evidence.

    Raises:
        ObjectiveStop: If any nonfinite sample is found.

    Side effects:
        Reads all signal samples; writes nothing.
    """

    step = max(1, int(round(chunk_seconds * raw.info["sfreq"])))
    nonfinite_count = 0
    chunks = 0
    for start in range(0, raw.n_times, step):
        data = raw.get_data(start=start, stop=min(raw.n_times, start + step))
        nonfinite_count += int(np.size(data) - np.isfinite(data).sum())
        chunks += 1
    if nonfinite_count:
        raise ObjectiveStop(
            "post_interpolation_validation",
            "nonfinite_post_reference_or_interpolation_signal",
            f"Continuous validation found {nonfinite_count} nonfinite samples.",
        )
    montage = raw.get_montage()
    montage_channels = set(montage.ch_names) if montage else set()
    missing_positions = sorted(set(EEG_TARGET_CHANNELS) - montage_channels)
    if missing_positions:
        raise ObjectiveStop(
            "post_interpolation_validation",
            "montage_missing_after_interpolation",
            f"Saved surface is missing montage channels: {missing_positions}",
        )
    return {
        "finite_sample_validation": "passed",
        "nonfinite_sample_count": 0,
        "sample_count": int(raw.n_times),
        "channel_count": len(raw.ch_names),
        "chunks_scanned": chunks,
        "operational_bad_list": list(raw.info["bads"]),
        "operational_bad_list_reset": not raw.info["bads"],
        "m1_m2_retained": all(channel in raw.ch_names for channel in MASTOID_CHANNELS),
        "montage_available": montage is not None,
    }
