"""Reusable helpers for the DEMI preprocessing-parameter evidence audit.

This module supports the final read-only evidence pass before production EEG
preprocessing. It prepares full continuous EEG recordings consistently, runs
global bad-channel candidates without changing an authoritative derivative,
summarizes PREP reference/interpolation state, and computes compact filter and
spectral comparison metrics.

Inputs are raw EDF recordings plus the public YAML contract and a private file
selection/configuration supplied by the calling scripts. Outputs are ordinary
Python records that the calling scripts may write only as CSV, JSON, or
Markdown audit evidence under ``_Data/``.

This module explicitly does not write EEG signals, modify EDF files, construct
epochs, fit or apply ICA, compute CSD, or implement the production preprocessing
pipeline. All signal changes occur only on in-memory copies used for comparison.

The command-line audit scripts import this module. It is also deliberately
small enough to test with synthetic ``mne.io.RawArray`` objects.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import mne
import numpy as np
import yaml
from scipy import signal


CRITERION_KEYS = (
    "bad_by_nan",
    "bad_by_flat",
    "bad_by_deviation",
    "bad_by_hf_noise",
    "bad_by_correlation",
    "bad_by_SNR",
    "bad_by_dropout",
    "bad_by_psd",
    "bad_by_ransac",
    "bad_by_manual",
)


def load_audit_contract(path: Path) -> dict[str, Any]:
    """Load and fail-closed validate the public preprocessing audit contract.

    Parameters
    ----------
    path
        YAML contract path.

    Returns
    -------
    dict
        Parsed and minimally validated contract.

    Raises
    ------
    ValueError
        If any safety boundary permits signal derivatives, epochs, or source
        mutation, or if the approved montage/input surface is not explicit.

    Side Effects
    ------------
    Reads one YAML file. It writes nothing.
    """

    with path.open("r", encoding="utf-8") as handle:
        contract = yaml.safe_load(handle)
    if not isinstance(contract, dict):
        raise ValueError("Audit contract must be a YAML mapping.")

    safety = contract.get("safety", {})
    required_false = (
        "write_signal_derivatives",
        "construct_epochs",
        "mutate_source_edf",
    )
    if safety.get("audit_only") is not True:
        raise ValueError("Audit contract must set safety.audit_only=true.")
    for key in required_false:
        if safety.get(key) is not False:
            raise ValueError(f"Audit contract must set safety.{key}=false.")
    if set(safety.get("allowed_output_suffixes", [])) != {".csv", ".json", ".md"}:
        raise ValueError("Audit outputs must be restricted to CSV, JSON, and Markdown.")

    prepared = contract.get("prepared_continuous_input", {})
    if prepared.get("montage") != "standard_1005":
        raise ValueError("The approved audit montage must be standard_1005.")
    if prepared.get("eeg_only") is not True:
        raise ValueError("Global detector input must be EEG-only.")
    for key in ("apply_reference", "apply_notch", "apply_analysis_filter"):
        if prepared.get(key) is not False:
            raise ValueError(f"Prepared detector input must set {key}=false.")
    if float(prepared.get("resample_hz", 0.0)) <= 120.0:
        raise ValueError("Audit resampling must retain the 60-Hz line frequency.")

    reference = contract.get("approved_reference_policy", {})
    expected_mastoids = ["M1", "M2"]
    if reference.get("decision_status") != "accepted":
        raise ValueError("The M1/M2 reference policy must be explicitly accepted.")
    for key in (
        "excluded_source_channels",
        "reference_target_channels",
        "retain_as_recorded_provenance_channels",
        "primary_feature_excluded_channels",
    ):
        if reference.get(key) != expected_mastoids:
            raise ValueError(f"approved_reference_policy.{key} must be [M1, M2].")
    if reference.get("primary_detector_source_pool") != "scalp_30":
        raise ValueError("The approved detector source pool must be scalp_30.")
    if reference.get("average_reference_source_pool") != "scalp_30":
        raise ValueError("The approved reference source pool must be scalp_30.")
    if reference.get("later_mastoid_sensitivity_permitted") is not True:
        raise ValueError("The approved policy must preserve mastoid sensitivity analyses.")
    if reference.get("deletion_permitted") is not False:
        raise ValueError("The approved policy must prohibit deleting M1/M2.")
    return contract


def validate_audit_output_path(path: Path, output_root: Path, suffixes: Sequence[str]) -> None:
    """Reject output paths outside the local evidence directory or with unsafe suffixes.

    Parameters
    ----------
    path
        Proposed output file.
    output_root
        Local audit evidence root.
    suffixes
        Explicit allowed suffixes from the contract.

    Returns
    -------
    None

    Side Effects
    ------------
    None. This function only validates paths.
    """

    resolved_root = output_root.resolve()
    resolved_path = path.resolve()
    if resolved_path.parent != resolved_root:
        raise ValueError(f"Audit output must be directly inside {resolved_root}: {path}")
    if resolved_path.suffix not in set(suffixes):
        raise ValueError(f"Unsafe audit output suffix: {resolved_path.suffix}")
    blocked_terms = ("epoch", "cleaned", "preprocessed", "derivative")
    if any(term in resolved_path.stem.lower() for term in blocked_terms):
        raise ValueError(f"Audit output name resembles a prohibited signal product: {path}")


def make_uppercase_standard_1005() -> mne.channels.DigMontage:
    """Return MNE's standard_1005 montage with DEMI-compatible uppercase labels.

    MNE's standard montage uses labels such as ``Fp1`` and ``Fz`` while DEMI
    EDF labels are uppercase. PyPREP calls ``set_montage`` internally without a
    case-insensitive option, so the montage itself must carry the exact labels.

    Returns
    -------
    mne.channels.DigMontage
        A copy of ``standard_1005`` whose channel labels are uppercase.

    Side Effects
    ------------
    None.
    """

    montage = mne.channels.make_standard_montage("standard_1005")
    rename = {name: name.upper() for name in montage.ch_names}
    if len(set(rename.values())) != len(rename):
        raise ValueError("Uppercasing standard_1005 unexpectedly creates duplicate labels.")
    montage.rename_channels(rename)
    return montage


def prepare_continuous_eeg(
    edf_path: Path,
    *,
    channel_types: Mapping[str, str],
    resample_hz: float,
) -> mne.io.BaseRaw:
    """Read and explicitly prepare one full continuous EEG detector input.

    Parameters
    ----------
    edf_path
        Source raw EDF path. The file is opened read-only.
    channel_types
        Required auxiliary/stim channel type overrides.
    resample_hz
        Audit-only target sampling rate. It must remain above 120 Hz.

    Returns
    -------
    mne.io.BaseRaw
        Preloaded EEG-only full-session data with the uppercase
        ``standard_1005`` montage and no reference/notch/analysis filter.

    Side Effects
    ------------
    Reads the EDF into memory. It never writes to or mutates the source file.
    """

    if resample_hz <= 120.0:
        raise ValueError("resample_hz must preserve 60-Hz content.")
    raw = mne.io.read_raw_edf(edf_path, preload=False, verbose="ERROR")
    missing = sorted(set(channel_types) - set(raw.ch_names))
    if missing:
        raise ValueError(f"Missing required channel-type labels in {edf_path.name}: {missing}")
    raw.set_channel_types(dict(channel_types), verbose="ERROR")
    raw.pick("eeg")
    if len(raw.ch_names) < 4:
        raise ValueError(f"Too few EEG channels after typing: {raw.ch_names}")
    raw.load_data(verbose="ERROR")
    raw.set_montage(make_uppercase_standard_1005(), on_missing="raise", verbose="ERROR")
    if not np.isclose(raw.info["sfreq"], resample_hz):
        raw.resample(resample_hz, verbose="ERROR")
    return raw


def parse_python_assignments(path: Path, names: Iterable[str]) -> dict[str, Any]:
    """Extract literal top-level Python assignments for a historical code audit.

    Parameters
    ----------
    path
        Python source path.
    names
        Assignment names to recover.

    Returns
    -------
    dict
        Recovered literal values keyed by requested assignment name.

    Side Effects
    ------------
    Reads source text. It does not import or execute the historical script.
    """

    requested = set(names)
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if isinstance(target, ast.Name) and target.id in requested:
            found[target.id] = ast.literal_eval(node.value)
    return found


def parse_r_scalar_assignments(path: Path, names: Iterable[str]) -> dict[str, Any]:
    """Extract simple logical/numeric scalar assignments from an R settings file.

    Parameters
    ----------
    path
        R source path.
    names
        Setting names to recover.

    Returns
    -------
    dict
        Boolean, integer, float, or raw string values.

    Side Effects
    ------------
    Reads source text only.
    """

    requested = set(names)
    found: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        code = line.split("#", 1)[0].strip()
        if "<-" not in code:
            continue
        key, value = (part.strip() for part in code.split("<-", 1))
        if key not in requested:
            continue
        if value in {"TRUE", "FALSE"}:
            found[key] = value == "TRUE"
        else:
            try:
                found[key] = int(value)
            except ValueError:
                try:
                    found[key] = float(value)
                except ValueError:
                    found[key] = value
    return found


def channel_list_json(channels: Iterable[str]) -> str:
    """Serialize a channel collection deterministically for audit CSV fields."""

    return json.dumps(sorted(set(channels)), separators=(",", ":"))


def stable_array_hash(values: np.ndarray) -> str:
    """Return a SHA-256 digest used only to compare repeated in-memory results."""

    contiguous = np.ascontiguousarray(values, dtype=np.float64)
    return hashlib.sha256(contiguous.tobytes()).hexdigest()


def run_lof(
    raw: mne.io.BaseRaw,
    *,
    n_neighbors: int,
    metric: str,
    threshold: float,
) -> dict[str, Any]:
    """Run MNE-native LOF and return bad labels plus channel-level scores.

    Parameters are passed directly to
    :func:`mne.preprocessing.find_bad_channels_lof`. The input is not modified.
    """

    bads, scores = mne.preprocessing.find_bad_channels_lof(
        raw,
        n_neighbors=n_neighbors,
        metric=metric,
        threshold=threshold,
        return_scores=True,
        verbose="ERROR",
    )
    return {
        "bad_all": sorted(bads),
        "scores": {name: float(score) for name, score in zip(raw.ch_names, scores)},
    }


def _criterion_score(noisy: Any, criterion: str, channel_index: int) -> dict[str, float]:
    """Extract available PyPREP diagnostic scores without assuming every key exists."""

    info = noisy._extra_info.get(criterion, {})  # PyPREP exposes no public score API.
    score_keys = {
        "bad_by_deviation": ("robust_channel_deviations",),
        "bad_by_hf_noise": ("hf_noise_zscores",),
        "bad_by_correlation": ("bad_window_fractions", "median_max_correlations"),
        "bad_by_dropout": ("bad_window_fractions",),
        "bad_by_psd": ("psd_zscore", "zscore_low", "zscore_mid", "zscore_high"),
        "bad_by_ransac": ("bad_window_fractions",),
    }.get(criterion, ())
    output: dict[str, float] = {}
    for key in score_keys:
        values = info.get(key)
        if values is None or len(values) <= channel_index:
            continue
        value = float(values[channel_index])
        output[key] = value if np.isfinite(value) else np.nan
    return output


def run_noisy_channels(raw: mne.io.BaseRaw, settings: Mapping[str, Any], seed: int) -> dict[str, Any]:
    """Run full-session PyPREP ``NoisyChannels`` without signal interpolation.

    Parameters
    ----------
    raw
        Prepared full-session EEG-only input.
    settings
        Explicit detrending, correlation, RANSAC, PSD, chunking, and MATLAB
        compatibility settings.
    seed
        Fixed RANSAC seed.

    Returns
    -------
    dict
        Criterion lists, union, and available per-channel diagnostics.

    Side Effects
    ------------
    PyPREP operates on an independent in-memory copy. No signal is written and
    the caller's raw object is not changed.
    """

    from pyprep.find_noisy_channels import NoisyChannels

    if bool(settings["psd_enabled"]) == bool(settings["matlab_strict"]):
        raise ValueError("PyPREP PSD is enabled exactly when matlab_strict is false.")
    noisy = NoisyChannels(
        raw.copy(),
        do_detrend=bool(settings["do_detrend"]),
        random_state=seed,
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
    criteria = noisy.get_bads(as_dict=True)
    details: list[dict[str, Any]] = []
    for index, channel in enumerate(raw.ch_names):
        row: dict[str, Any] = {"channel": channel}
        for criterion in CRITERION_KEYS:
            row[criterion] = channel in criteria.get(criterion, [])
            for score_name, score in _criterion_score(noisy, criterion, index).items():
                row[f"{criterion}_{score_name}"] = score
        details.append(row)
    return {
        "criteria": {key: sorted(criteria.get(key, [])) for key in CRITERION_KEYS},
        "bad_all": sorted(criteria.get("bad_all", [])),
        "channel_details": details,
    }


def run_full_prep(raw: mne.io.BaseRaw, settings: Mapping[str, Any], seed: int) -> dict[str, Any]:
    """Evaluate the integrated PyPREP pipeline on an independent signal copy.

    The returned record separates original detector findings, unusable/bad
    channels before interpolation, interpolated channels, still-noisy channels,
    and robust-reference hashes. No processed signal is returned or written.
    """

    from pyprep.prep_pipeline import PrepPipeline

    if settings.get("interpolate_bads") is not True:
        raise ValueError("The integrated comparison must evaluate PREP interpolation.")
    montage = make_uppercase_standard_1005()
    prep_params = {
        "ref_chs": settings["ref_chs"],
        "reref_chs": settings["reref_chs"],
        "line_freqs": [float(value) for value in settings["line_freqs_hz"]],
        "max_iterations": int(settings["max_iterations"]),
    }
    pipeline = PrepPipeline(
        raw.copy(),
        prep_params,
        montage,
        ransac=bool(settings["ransac"]),
        channel_wise=bool(settings["channel_wise"]),
        max_chunk_size=settings.get("max_chunk_size"),
        random_state=seed,
        matlab_strict=bool(settings["matlab_strict"]),
    )
    pipeline.fit()
    original = pipeline.noisy_channels_original or {}
    before = pipeline.noisy_channels_before_interpolation or {}
    after = pipeline.noisy_channels_after_interpolation or {}
    interpolated = sorted(pipeline.interpolated_channels or [])
    still_noisy = sorted(pipeline.still_noisy_channels or [])
    bad_before = sorted(pipeline.bad_before_interpolation or [])
    unusable = sorted(
        set(before.get("bad_by_nan", []))
        | set(before.get("bad_by_flat", []))
        | set(before.get("bad_by_manual", []))
    )
    return {
        "original_criteria": {key: sorted(original.get(key, [])) for key in CRITERION_KEYS},
        "original_bad_all": sorted(original.get("bad_all", [])),
        "before_interpolation_criteria": {
            key: sorted(before.get(key, [])) for key in CRITERION_KEYS
        },
        "after_interpolation_criteria": {
            key: sorted(after.get(key, [])) for key in CRITERION_KEYS
        },
        "unusable_channels": unusable,
        "bad_before_interpolation": bad_before,
        "interpolated_channels": interpolated,
        "still_noisy_channels": still_noisy,
        "interpolated_count": len(interpolated),
        "interpolated_proportion": len(interpolated) / len(raw.ch_names),
        "reference_kind": "iterative_robust_average_reference",
        "reference_channels": list(raw.ch_names),
        "rereference_channels": list(raw.ch_names),
        "reference_before_hash": stable_array_hash(pipeline.reference_before_interpolation),
        "reference_after_hash": stable_array_hash(pipeline.reference_after_interpolation),
        "reference_before_rms_v": float(
            np.sqrt(np.mean(np.square(pipeline.reference_before_interpolation)))
        ),
        "reference_after_rms_v": float(
            np.sqrt(np.mean(np.square(pipeline.reference_after_interpolation)))
        ),
    }


def welch_band_powers(
    data: np.ndarray,
    sfreq: float,
    *,
    bands: Mapping[str, Sequence[float]],
    n_fft: int,
    n_overlap: int,
) -> dict[str, np.ndarray]:
    """Compute channel-wise Welch mean PSD within named frequency bands.

    Parameters
    ----------
    data
        Channels by samples array.
    sfreq
        Sampling rate in Hz.
    bands
        Named inclusive frequency limits.
    n_fft, n_overlap
        Welch segment parameters.

    Returns
    -------
    dict
        One channel-length mean-PSD array per band.
    """

    nperseg = min(int(n_fft), data.shape[1])
    noverlap = min(int(n_overlap), nperseg - 1)
    freqs, psd = signal.welch(
        data,
        fs=sfreq,
        nperseg=nperseg,
        noverlap=noverlap,
        axis=-1,
        detrend="constant",
    )
    output: dict[str, np.ndarray] = {}
    for name, limits in bands.items():
        low, high = (float(limits[0]), float(limits[1]))
        mask = (freqs >= low) & (freqs <= high)
        if not np.any(mask):
            raise ValueError(f"No Welch bins for {name}: {low}-{high} Hz")
        output[name] = np.mean(psd[:, mask], axis=1)
    return output


def analysis_filter_response(
    sfreq: float,
    *,
    high_pass_hz: float,
    low_pass_hz: float,
) -> dict[str, Any]:
    """Create the proposed MNE FIR and summarize response at relevant frequencies."""

    taps = mne.filter.create_filter(
        None,
        sfreq,
        l_freq=high_pass_hz,
        h_freq=low_pass_hz,
        method="fir",
        phase="zero",
        fir_design="firwin",
        verbose="ERROR",
    )
    frequencies, response = signal.freqz(taps, worN=32768, fs=sfreq)
    response_db = 20.0 * np.log10(np.maximum(np.abs(response), np.finfo(float).tiny))

    def at(freq: float) -> float:
        return float(response_db[np.argmin(np.abs(frequencies - freq))])

    return {
        "filter_length_samples": int(len(taps)),
        "filter_length_seconds": float(len(taps) / sfreq),
        "response_db_0_1_hz": at(0.1),
        "response_db_0_5_hz": at(0.5),
        "response_db_4_hz": at(4.0),
        "response_db_8_hz": at(8.0),
        "response_db_12_hz": at(12.0),
        "response_db_30_hz": at(30.0),
        "response_db_45_hz": at(45.0),
        "response_db_60_hz": at(60.0),
    }


def compare_line_noise_branches(
    raw: mne.io.BaseRaw,
    settings: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Compare no-notch and explicit-60-Hz-notch branches in memory.

    Both branches receive the same proposed 0.5--45 Hz MNE FIR after the
    optional notch decision. Metrics exclude a fixed edge margin. The function
    returns channel-level spectral/change rows and the exact FIR response; it
    writes no signal derivative.
    """

    sfreq = float(raw.info["sfreq"])
    high_pass = float(settings["analysis_high_pass_hz"])
    low_pass = float(settings["analysis_low_pass_hz"])
    edge = int(round(float(settings["measurement_edge_exclusion_seconds"]) * sfreq))
    if raw.n_times <= edge * 2:
        raise ValueError("Recording is too short for the configured edge exclusion.")
    bands = settings["bands_hz"]
    n_fft = int(settings["welch_n_fft"])
    n_overlap = int(settings["welch_n_overlap"])

    raw_data = raw.get_data()[:, edge:-edge]
    raw_powers = welch_band_powers(
        raw_data, sfreq, bands=bands, n_fft=n_fft, n_overlap=n_overlap
    )

    branch_a = raw.copy().filter(
        high_pass,
        low_pass,
        method="fir",
        phase=settings["analysis_phase"],
        fir_design=settings["analysis_fir_design"],
        verbose="ERROR",
    )
    data_a = branch_a.get_data()[:, edge:-edge]
    powers_a = welch_band_powers(data_a, sfreq, bands=bands, n_fft=n_fft, n_overlap=n_overlap)

    branch_b = raw.copy().notch_filter(
        freqs=[float(settings["line_frequency_hz"])],
        method=settings["notch_method"],
        filter_length=settings["notch_filter_length"],
        mt_bandwidth=float(settings["notch_mt_bandwidth_hz"]),
        p_value=float(settings["notch_p_value"]),
        verbose="ERROR",
    )
    notch_data = branch_b.get_data()[:, edge:-edge]
    notch_powers = welch_band_powers(
        notch_data, sfreq, bands=bands, n_fft=n_fft, n_overlap=n_overlap
    )
    branch_b.filter(
        high_pass,
        low_pass,
        method="fir",
        phase=settings["analysis_phase"],
        fir_design=settings["analysis_fir_design"],
        verbose="ERROR",
    )
    data_b = branch_b.get_data()[:, edge:-edge]
    powers_b = welch_band_powers(data_b, sfreq, bands=bands, n_fft=n_fft, n_overlap=n_overlap)

    rows: list[dict[str, Any]] = []
    for index, channel in enumerate(raw.ch_names):
        rms_a = float(np.sqrt(np.mean(np.square(data_a[index]))))
        difference = data_b[index] - data_a[index]
        rms_difference = float(np.sqrt(np.mean(np.square(difference))))
        correlation = float(np.corrcoef(data_a[index], data_b[index])[0, 1])
        row: dict[str, Any] = {
            "channel": channel,
            "analysis_a_rms_v": rms_a,
            "analysis_b_minus_a_rms_v": rms_difference,
            "analysis_b_minus_a_rms_proportion": rms_difference / rms_a if rms_a else np.nan,
            "analysis_a_b_correlation": correlation,
        }
        for band in bands:
            raw_value = float(raw_powers[band][index])
            notch_value = float(notch_powers[band][index])
            a_value = float(powers_a[band][index])
            b_value = float(powers_b[band][index])
            row[f"raw_{band}_mean_psd"] = raw_value
            row[f"notch_only_{band}_mean_psd"] = notch_value
            row[f"analysis_a_{band}_mean_psd"] = a_value
            row[f"analysis_b_{band}_mean_psd"] = b_value
            row[f"analysis_b_vs_a_{band}_relative_change"] = (
                (b_value - a_value) / a_value if a_value else np.nan
            )
        rows.append(row)

    response = analysis_filter_response(
        sfreq, high_pass_hz=high_pass, low_pass_hz=low_pass
    )
    response["measurement_edge_exclusion_seconds"] = float(
        settings["measurement_edge_exclusion_seconds"]
    )
    return rows, response
