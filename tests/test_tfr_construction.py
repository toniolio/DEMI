"""Focused tests for the accepted DEMI trial-level TFR construction contract."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import mne
import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from tfr_construction import (  # noqa: E402
    ALL_FAMILIES,
    ArrayContract,
    BASELINE_TMAX,
    BASELINE_TMIN,
    CALCULATION_DTYPE,
    CROSS_FAMILY_IDENTITY_COLUMNS,
    EXPECTED_DURATION_WARNING_COUNT,
    EXPECTED_FILE49_COUNT,
    EXPECTED_STRICT_CLEAN_COUNT,
    EXPECTED_TASK_SAMPLE_COUNT,
    EXPECTED_TASK_TIMES,
    EXPECTED_TRIAL_COUNT,
    EXCLUDED_PRIMARY_CHANNELS,
    FREQUENCIES_HZ,
    INPUT_SFREQ,
    N_CYCLES,
    PERSISTED_DTYPE,
    PRIMARY_SCALP_CHANNELS,
    RESAMPLE_METHOD,
    RESAMPLE_PAD,
    RESAMPLE_WINDOW,
    TARGET_SFREQ,
    TASK_TMAX,
    TASK_TMIN,
    TFR_DECIM,
    artifact_diagnostics,
    assess_cached_recording,
    atomic_write_json,
    atomic_write_npy,
    compare_source_snapshots,
    compute_morlet_power,
    crop_task_power,
    derive_baseline_mean,
    derive_log_power,
    expected_wavelet_support,
    forbidden_output_scan,
    load_and_validate_config,
    normalize_db,
    reopen_validate_npy,
    require_ignored_output_root,
    resample_for_tfr,
    sha256_file,
    source_snapshot,
    validate_family_identity,
    validate_stored_formula,
)


CONFIG_PATH = EEG_DIR / "tfr_config_v1.yaml"
EPOCH_ROOT = REPO_ROOT / "_Data" / "eeg" / "epochs_v1"


def identity_metadata(family: str, count: int = 2) -> pd.DataFrame:
    """Return exact cross-family identity metadata plus family anchor fields."""

    rows = []
    for index in range(count):
        row = {}
        for column in CROSS_FAMILY_IDENTITY_COLUMNS:
            if column in {
                "canonical_order_index",
                "behavioural_id",
                "eeg_source_id",
                "offset_session",
                "offset_block",
                "offset_trial",
                "audit_trial_count",
                "raw_trial_sequence",
                "join_trial_count",
                "offset_trial_count",
                "interpolation_count",
                "interpolation_denominator",
            }:
                row[column] = index if column == "canonical_order_index" else 1
            elif column in {
                "physical",
                "primary_eligibility",
                "strict_clean_eligibility",
                "duration_warning_flag",
            }:
                row[column] = column != "duration_warning_flag"
            elif column == "interpolation_proportion":
                row[column] = 0.0
            elif column == "canonical_event_key":
                row[column] = f"key-{index}"
            elif column == "condition_semantics":
                row[column] = "overt_movement"
            elif column == "source_recording_filename":
                row[column] = "demi_01 Data.edf"
            else:
                row[column] = f"{column}-value"
        row["derived_anchor_type"] = family
        rows.append(row)
    return pd.DataFrame(rows)


def synthetic_epochs(
    family: str,
    *,
    count: int = 2,
    tmin: float = -1.5,
    tmax: float = 2.5,
) -> mne.EpochsArray:
    """Return accepted-shaped finite 1000-Hz synthetic epochs."""

    channels = [
        *PRIMARY_SCALP_CHANNELS[:17],
        "M1",
        *PRIMARY_SCALP_CHANNELS[17:23],
        "M2",
        *PRIMARY_SCALP_CHANNELS[23:],
        "HEO",
        "VEO",
        "EMG-L",
        "EMG-A",
        "Trigger",
    ]
    channel_types = (
        ["eeg"] * 17
        + ["eeg"]
        + ["eeg"] * 6
        + ["eeg"]
        + ["eeg"] * 7
        + ["eog", "eog", "emg", "emg", "stim"]
    )
    sample_count = int(round((tmax - tmin) * INPUT_SFREQ)) + 1
    rng = np.random.default_rng(20260720)
    data = rng.normal(scale=5e-6, size=(count, len(channels), sample_count))
    # Keep Trigger physically harmless while retaining the accepted type.
    data[:, -1] = 0.0
    info = mne.create_info(channels, INPUT_SFREQ, channel_types)
    events = np.column_stack(
        [np.arange(count, dtype=int), np.zeros(count, dtype=int), np.ones(count, dtype=int)]
    )
    return mne.EpochsArray(
        data,
        info,
        events=events,
        event_id={family: 1},
        tmin=tmin,
        baseline=None,
        metadata=identity_metadata(family, count),
        verbose="ERROR",
    )


def test_exact_configuration_frequency_cycle_channel_and_safety_contract() -> None:
    """Tracked configuration pins all accepted scientific and non-operation values."""

    config, digest = load_and_validate_config(CONFIG_PATH)
    assert len(digest) == 64
    assert np.array_equal(FREQUENCIES_HZ, np.arange(4.0, 41.0, 1.0))
    assert len(FREQUENCIES_HZ) == 37
    assert np.allclose(N_CYCLES, 4.0 + (FREQUENCIES_HZ - 4.0) / 6.0)
    assert N_CYCLES[[0, 4, 9, 26, 36]].tolist() == pytest.approx(
        [4.0, 4.666666666666667, 5.5, 8.333333333333334, 10.0]
    )
    assert len(PRIMARY_SCALP_CHANNELS) == 30
    assert not set(PRIMARY_SCALP_CHANNELS) & set(EXCLUDED_PRIMARY_CHANNELS)
    assert config["tfr"]["output"] == "power"
    assert TFR_DECIM == 1
    assert not config["safety"]["epoch_interpolation"]
    assert not config["safety"]["epoch_rejection"]
    assert not config["safety"]["run_autoreject"]


def test_resampling_is_explicit_polyphase_1000_to_100_and_preserves_grid() -> None:
    """A copied 30-channel object is anti-alias resampled before TFR."""

    epochs = synthetic_epochs("response_onset", count=1)
    before = epochs.get_data(copy=True)
    work = resample_for_tfr(epochs)
    assert RESAMPLE_METHOD == "polyphase"
    assert RESAMPLE_WINDOW == ("kaiser", 5.0)
    assert RESAMPLE_PAD == "reflect"
    assert epochs.info["sfreq"] == INPUT_SFREQ
    assert epochs.ch_names != work.ch_names
    assert np.array_equal(epochs.get_data(copy=False), before)
    assert work.info["sfreq"] == TARGET_SFREQ
    assert work.ch_names == list(PRIMARY_SCALP_CHANNELS)
    assert len(work.times) == 400
    assert work.times[0] == -1.5
    assert work.times[-1] == 2.49


def test_morlet_axis_order_time_crop_and_float64_boundary() -> None:
    """Power axes are trial/channel/frequency/time with exact retained times."""

    work = resample_for_tfr(synthetic_epochs("response_onset", count=1))
    power, times = compute_morlet_power(work)
    assert power.dtype == CALCULATION_DTYPE
    assert power.shape == (1, 30, 37, 400)
    retained, retained_times = crop_task_power(power, times)
    assert retained.shape == (1, 30, 37, EXPECTED_TASK_SAMPLE_COUNT)
    assert retained.dtype == CALCULATION_DTYPE
    assert np.allclose(retained_times, EXPECTED_TASK_TIMES, atol=1e-12, rtol=0.0)
    assert retained_times[[0, -1]].tolist() == [TASK_TMIN, TASK_TMAX]


def test_wavelet_support_keeps_task_and_baseline_edge_safe() -> None:
    """Installed-MNE kernels preserve the accepted response and red_on support."""

    support = expected_wavelet_support()
    assert support["maximum_half_support_seconds"] == pytest.approx(0.79)
    assert support["wavelet_samples"][0] == 159
    assert support["task_edge_safe_interval_seconds"][0] <= TASK_TMIN
    assert support["task_edge_safe_interval_seconds"][1] >= TASK_TMAX
    assert support["red_on_edge_safe_interval_seconds"][0] <= BASELINE_TMIN
    assert support["red_on_edge_safe_interval_seconds"][1] >= BASELINE_TMAX


def test_canonical_matching_is_one_to_one_and_provenance_exact() -> None:
    """Shared baselines cannot be matched by array order alone."""

    frames = {family: identity_metadata(family) for family in ALL_FAMILIES}
    digest = validate_family_identity(frames)
    assert len(digest) == 64
    frames["red_on"].loc[0, "canonical_event_key"] = "different"
    with pytest.raises(RuntimeError, match="identity_or_provenance_mismatch"):
        validate_family_identity(frames)
    frames = {family: identity_metadata(family) for family in ALL_FAMILIES}
    frames["red_on"].loc[1, "canonical_event_key"] = "key-0"
    with pytest.raises(RuntimeError, match="duplicate"):
        validate_family_identity(frames)


def test_baseline_is_inclusive_mean_power_before_ratio_and_shared() -> None:
    """Baseline uses 31 samples, time-mean power, then the exact dB ratio."""

    times = np.arange(-150, 81, dtype=float) / 100.0
    red = np.ones((2, 1, 1, len(times)), dtype=np.float64)
    baseline_indices = (times >= BASELINE_TMIN) & (times <= BASELINE_TMAX)
    red[0, 0, 0, baseline_indices] = np.arange(1.0, 32.0)
    red[1, 0, 0, baseline_indices] = 4.0
    baseline, evidence = derive_baseline_mean(red, times)
    assert evidence["sample_count"] == 31
    assert evidence["first_time_seconds"] == BASELINE_TMIN
    assert evidence["last_time_seconds"] == BASELINE_TMAX
    assert baseline[:, 0, 0].tolist() == [16.0, 4.0]

    task = np.full((2, 1, 1, 3), 64.0, dtype=np.float64)
    onset = normalize_db(task, baseline)
    end = normalize_db(task.copy(), baseline)
    assert np.array_equal(onset, end)
    assert onset[0, 0, 0, 0] == pytest.approx(10.0 * np.log10(64.0 / 16.0))
    pointwise_wrong = np.mean(10.0 * np.log10(64.0 / np.arange(1.0, 32.0)))
    assert onset[0, 0, 0, 0] != pytest.approx(pointwise_wrong)


@pytest.mark.parametrize(
    "bad_value,reason",
    [(0.0, "nonpositive"), (-1.0, "nonpositive"), (np.nan, "nonfinite"), (np.inf, "nonfinite")],
)
def test_invalid_baselines_fail_without_clipping(bad_value: float, reason: str) -> None:
    """Zero, negative, and non-finite denominators fail rather than receive epsilon."""

    task = np.ones((1, 1, 1, 1), dtype=np.float64)
    baseline = np.asarray([[[bad_value]]], dtype=np.float64)
    with pytest.raises(RuntimeError, match=reason):
        normalize_db(task, baseline)


def test_log_power_is_exact_deterministic_view_without_second_convolution() -> None:
    """The required sensitivity is exactly 10log10 of stored raw power."""

    raw = np.asarray([1.0, 10.0, 100.0], dtype=np.float32)
    log_power = derive_log_power(raw)
    assert log_power.dtype == np.dtype("float64")
    assert log_power.tolist() == pytest.approx([0.0, 10.0, 20.0])
    with pytest.raises(ValueError, match="strictly positive"):
        derive_log_power(np.asarray([0.0], dtype=np.float32))


def test_diagnostics_are_metadata_only_and_include_auxiliary_summaries() -> None:
    """Accepted p2p/robust/jump/flat/EOG/EMG facts never change eligibility."""

    frame = artifact_diagnostics(synthetic_epochs("response_end"), "response_end")
    assert len(frame) == 2
    assert frame["diagnostic_only"].all()
    assert not frame["changes_primary_eligibility"].any()
    required = {
        "scalp_p2p_max_uv",
        "diagnostic_robust_logp2p_z_gt_6",
        "diagnostic_jump_gt_50uv_per_native_sample",
        "diagnostic_flat_p2p_lt_1uv",
        "eog_p2p_max_uv",
        "emg_p2p_max_uv",
        "continuous_qc_warning",
        "condition_semantics",
        "canonical_event_key",
    }
    assert required <= set(frame.columns)


def test_atomic_float32_npy_reopen_formula_and_cache_reuse(tmp_path: Path) -> None:
    """Atomic arrays reopen independently and complete matching shards are reusable."""

    recording = tmp_path / "recording"
    onset = recording / "response_onset"
    end = recording / "response_end"
    red = recording / "red_on"
    shape = (2, 2, 2, 3)
    baseline_shape = shape[:3]
    raw64 = np.arange(1, np.prod(shape) + 1, dtype=np.float64).reshape(shape)
    baseline64 = np.full(baseline_shape, 2.0, dtype=np.float64)
    db64 = normalize_db(raw64, baseline64)
    arrays = {}
    for family, directory in (("response_onset", onset), ("response_end", end)):
        raw_desc = atomic_write_npy(directory / "raw_power.npy", raw64.astype(np.float32))
        db_desc = atomic_write_npy(directory / "db_power.npy", db64.astype(np.float32))
        raw_desc["relative_path"] = f"{family}/raw_power.npy"
        db_desc["relative_path"] = f"{family}/db_power.npy"
        arrays[f"{family}_raw"] = raw_desc
        arrays[f"{family}_db"] = db_desc
    baseline_desc = atomic_write_npy(red / "baseline_mean_power.npy", baseline64.astype(np.float32))
    baseline_desc["relative_path"] = "red_on/baseline_mean_power.npy"
    arrays["red_on_baseline_mean"] = baseline_desc
    contracts = {
        "response_onset_raw": ArrayContract(shape, PERSISTED_DTYPE, ("trial", "channel", "frequency", "time"), positive=True),
        "response_onset_db": ArrayContract(shape, PERSISTED_DTYPE, ("trial", "channel", "frequency", "time")),
        "response_end_raw": ArrayContract(shape, PERSISTED_DTYPE, ("trial", "channel", "frequency", "time"), positive=True),
        "response_end_db": ArrayContract(shape, PERSISTED_DTYPE, ("trial", "channel", "frequency", "time")),
        "red_on_baseline_mean": ArrayContract(baseline_shape, PERSISTED_DTYPE, ("trial", "channel", "frequency"), positive=True),
    }
    atomic_write_json(
        recording / "manifest.json",
        {"status": "complete", "fingerprint": "fingerprint", "arrays": arrays},
    )
    assert assess_cached_recording(recording, "fingerprint", contracts)["action"] == "reuse"
    validation = reopen_validate_npy(
        onset / "raw_power.npy",
        contracts["response_onset_raw"],
        expected_sha256=sha256_file(onset / "raw_power.npy"),
    )
    assert validation["dtype"] == "float32"
    formula = validate_stored_formula(
        onset / "raw_power.npy",
        onset / "db_power.npy",
        red / "baseline_mean_power.npy",
    )
    assert formula["maximum_absolute_db_difference"] < 2e-4
    with (onset / "raw_power.npy").open("ab") as handle:
        handle.write(b"x")
    assert assess_cached_recording(recording, "fingerprint", contracts)["action"] == "rebuild"


def test_source_immutability_and_output_namespace_safety(tmp_path: Path) -> None:
    """Size/mtime mutation is detected and production output must be ignored."""

    source = tmp_path / "source-epo.fif"
    source.write_bytes(b"accepted")
    before = source_snapshot([source])
    assert compare_source_snapshots(before, source_snapshot([source]))["unchanged"]
    source.write_bytes(b"changed")
    with pytest.raises(RuntimeError, match="changed"):
        compare_source_snapshots(before, source_snapshot([source]))
    require_ignored_output_root(REPO_ROOT, REPO_ROOT / "_Data/eeg/tfr_v1")
    with pytest.raises(RuntimeError, match="exactly"):
        require_ignored_output_root(REPO_ROOT, tmp_path)


def test_forbidden_output_scan_blocks_scope_expansion(tmp_path: Path) -> None:
    """No AutoReject, CSD, phase, ROI, band, model, ID86, or 54_1 output is allowed."""

    (tmp_path / "raw_power.npy").touch()
    assert forbidden_output_scan(tmp_path) == []
    (tmp_path / "autoreject_labels.csv").touch()
    (tmp_path / "id86_tfr.npy").touch()
    assert forbidden_output_scan(tmp_path) == ["autoreject_labels.csv", "id86_tfr.npy"]


@pytest.mark.skipif(not EPOCH_ROOT.is_dir(), reason="local accepted epochs unavailable")
def test_local_authority_exact_counts_channels_and_special_routes() -> None:
    """Local manifests retain 8798/8789/nine, file 49, and exclude 54_1/86."""

    manifests = {
        family: json.loads(
            (EPOCH_ROOT / "manifests" / f"{family}_manifest.json").read_text()
        )
        for family in ALL_FAMILIES
    }
    for family, manifest in manifests.items():
        assert manifest["epoch_count"] == EXPECTED_TRIAL_COUNT
        assert manifest["strict_clean_count"] == EXPECTED_STRICT_CLEAN_COUNT
        assert manifest["duration_warning_count"] == EXPECTED_DURATION_WARNING_COUNT
        assert len(manifest["shards"]) == 81
        assert sum(
            shard["row_count"]
            for shard in manifest["shards"]
            if shard["recording"] == "demi_49_data"
        ) == EXPECTED_FILE49_COUNT
        recordings = {shard["recording"] for shard in manifest["shards"]}
        assert "demi_49_data" in recordings
        assert not any("54_1" in item for item in recordings)
        assert not any(item in {"demi_86", "demi_86_data"} for item in recordings)
    sample_path = EPOCH_ROOT / manifests["response_onset"]["shards"][0]["artifact"]["relative_path"]
    epochs = mne.read_epochs(sample_path, preload=False, verbose="ERROR")
    observed_primary = [name for name in epochs.ch_names if name in PRIMARY_SCALP_CHANNELS]
    assert observed_primary == list(PRIMARY_SCALP_CHANNELS)
    assert not set(EXCLUDED_PRIMARY_CHANNELS) & set(observed_primary)
