"""Focused tests for the accepted DEMI stage-17 EEG feature contract."""

from __future__ import annotations

from pathlib import Path
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from feature_construction import (  # noqa: E402
    BAND_LIMITS,
    EXPECTED_CHANNEL_DEFINITIONS,
    EXPECTED_CHANNEL_ROWS,
    EXPECTED_DURATION_WARNING_COUNT,
    EXPECTED_FILE49_COUNT,
    EXPECTED_PARTICIPANT_COUNT,
    EXPECTED_ROI_DEFINITIONS,
    EXPECTED_ROI_ROWS,
    EXPECTED_SHORT_0P5_COUNT,
    EXPECTED_SHORT_1P0_COUNT,
    EXPECTED_STRICT_CLEAN_COUNT,
    EXPECTED_TRIAL_COUNT,
    FIXED_ROI_CHANNELS,
    LEFT_TASK_CONTRALATERAL,
    LEFT_TASK_IPSILATERAL,
    PRIMARY_CHANNELS,
    RIGHT_TASK_CONTRALATERAL,
    RIGHT_TASK_IPSILATERAL,
    WINDOW_LIMITS,
    apply_task_hand_mapping,
    atomic_combine_parquet,
    atomic_write_parquet,
    compare_source_snapshots,
    definition_table,
    derive_roi_rows,
    feature_definitions,
    inclusive_indices,
    load_and_validate_config,
    require_ignored_output_root,
    roi_definition_table,
    roi_equality_validation,
    source_snapshot,
    standard_1005_channel_table,
    summarize_transformed_cells,
)


CONFIG_PATH = EEG_DIR / "feature_config_v1.yaml"
PRODUCTION_ROOT = REPO_ROOT / "_Data/eeg/features_v1"
PRODUCTION_MANIFEST = PRODUCTION_ROOT / "feature_run_manifest.json"


def load_config() -> dict:
    """Return the validated tracked feature configuration."""

    return load_and_validate_config(CONFIG_PATH)[0]


def test_exact_band_membership_and_no_exploratory_or_beta_split() -> None:
    """Bands are inclusive 4--30 Hz equal-frequency definitions only."""

    config = load_config()
    frequencies = np.arange(4.0, 41.0)
    expected = {
        "theta": np.arange(4.0, 9.0),
        "alpha": np.arange(9.0, 13.0),
        "beta": np.arange(13.0, 31.0),
    }
    for band, members in expected.items():
        minimum, maximum, count = BAND_LIMITS[band]
        indices = inclusive_indices(
            frequencies, minimum, maximum, expected_count=count, label=band
        )
        assert np.array_equal(frequencies[indices], members)
    assert set(config["bands"]) == {"theta", "alpha", "beta"}
    assert config["diagnostic_frequency_surface"]["feature_rows_created"] is False
    assert max(item.maximum_hz for item in feature_definitions(config)) == 30.0


def test_exact_window_endpoint_and_sample_membership() -> None:
    """Every fixed 100-Hz window includes both declared endpoints."""

    config = load_config()
    times = np.arange(-50, 151, dtype=float) / 100.0
    frequencies = np.arange(4.0, 41.0)
    table = definition_table(config, frequencies, times)
    for name, (_, minimum, maximum, count, bands) in WINDOW_LIMITS.items():
        members = inclusive_indices(
            times, minimum, maximum, expected_count=count, label=name
        )
        assert times[members[0]] == minimum
        assert times[members[-1]] == maximum
        rows = table.loc[table["feature_definition"].eq(name)]
        assert len(rows) == len(bands)
        assert rows["time_sample_count"].eq(count).all()


def test_transformed_cell_arithmetic_means_and_log_power_order() -> None:
    """dB is averaged as stored; log power is transformed cell-by-cell first."""

    db = np.arange(1, 1 + 2 * 2 * 3 * 4, dtype=float).reshape(2, 2, 3, 4)
    raw = np.power(10.0, db / 10.0)
    frequency_indices = np.array([0, 1, 2])
    time_indices = np.array([1, 2, 3])
    observed_db, observed_log = summarize_transformed_cells(
        db, raw, frequency_indices, time_indices
    )
    expected = db[:, :, :, 1:4].mean(axis=(2, 3))
    assert np.allclose(observed_db, expected, atol=1e-12, rtol=0.0)
    assert np.allclose(observed_log, expected, atol=1e-12, rtol=0.0)
    # This deliberately differs, proving raw power is not averaged before log.
    wrong_order = 10.0 * np.log10(raw[:, :, :, 1:4].mean(axis=(2, 3)))
    assert not np.allclose(observed_log, wrong_order)


def test_task_hand_mapping_preserves_ambidextrous_code_and_owner_override() -> None:
    """Participant 70 keeps `a` but receives the right-task-hand mapping."""

    trials = pd.DataFrame(
        {
            "behavioural_id": [1, 28, 70],
            "handedness": ["r", "l", "a"],
        }
    )
    mapped = apply_task_hand_mapping(trials, load_config())
    assert mapped["original_handedness"].tolist() == ["r", "l", "a"]
    assert mapped["analysis_hand"].tolist() == ["right", "left", "right"]
    p70 = mapped.loc[mapped["behavioural_id"].eq(70)].iloc[0]
    assert p70["motor_laterality_mapping"] == "right_hand_task"
    assert p70["mapping_source"] == "owner_recollection"
    assert bool(p70["mapping_override"])


def synthetic_channel_rows() -> pd.DataFrame:
    """Return all channel definitions for one right, left, and P70 task hand."""

    config = load_config()
    rows = []
    trial_identity = [
        (0, "right-key", "right"),
        (1, "left-key", "left"),
        (2, "p70-key", "right"),
    ]
    for definition in feature_definitions(config):
        for order, key, hand in trial_identity:
            for channel_index, channel in enumerate(PRIMARY_CHANNELS):
                rows.append(
                    {
                        "canonical_order_index": order,
                        "canonical_event_key": key,
                        "feature_definition": definition.name,
                        "band": definition.band,
                        "analysis_hand": hand,
                        "physical_channel": channel,
                        "channel_index": channel_index,
                        "montage": "standard_1005",
                        "coordinate_frame": "head",
                        "x_m": 0.0,
                        "y_m": 0.0,
                        "z_m": 0.0,
                        "value_db": order * 100.0 + channel_index,
                        "value_log_power": order * 1000.0 + 2.0 * channel_index,
                    }
                )
    return pd.DataFrame(rows)


def test_exact_roi_membership_task_laterality_and_equal_weight_derivation() -> None:
    """Fixed and task-hand-normalized ROIs are exact unweighted channel means."""

    config = load_config()
    assert FIXED_ROI_CHANNELS["frontal_theta"] == ("F3", "FZ", "F4", "FCZ")
    assert FIXED_ROI_CHANNELS["posterior_alpha"] == ("P3", "PZ", "P4", "O1", "OZ", "O2")
    assert RIGHT_TASK_CONTRALATERAL == ("FC3", "C3", "CP3")
    assert RIGHT_TASK_IPSILATERAL == ("FC4", "C4", "CP4")
    assert LEFT_TASK_CONTRALATERAL == RIGHT_TASK_IPSILATERAL
    assert LEFT_TASK_IPSILATERAL == RIGHT_TASK_CONTRALATERAL

    channel_rows = synthetic_channel_rows()
    roi_rows = derive_roi_rows(channel_rows, config)
    assert len(roi_rows) == 3 * EXPECTED_ROI_DEFINITIONS
    equality = roi_equality_validation(channel_rows, roi_rows)
    assert equality["maximum_absolute_db_difference"] == 0.0
    assert equality["maximum_absolute_log_power_difference"] == 0.0

    right = roi_rows.query(
        "canonical_event_key == 'right-key' and roi_name == 'contralateral_motor_beta'"
    )
    left = roi_rows.query(
        "canonical_event_key == 'left-key' and roi_name == 'contralateral_motor_beta'"
    )
    p70 = roi_rows.query(
        "canonical_event_key == 'p70-key' and roi_name == 'contralateral_motor_beta'"
    )
    assert set(right["physical_source_channels"]) == {"FC3;C3;CP3"}
    assert set(left["physical_source_channels"]) == {"FC4;C4;CP4"}
    assert set(p70["physical_source_channels"]) == {"FC3;C3;CP3"}
    assert roi_rows["roi_weighting"].eq("equal_weight_arithmetic_mean").all()


def test_motor_roi_derivation_supports_a_single_task_hand_recording() -> None:
    """A recording containing only one task hand still yields complete motor ROIs."""

    channel_rows = synthetic_channel_rows()
    right_only = channel_rows.loc[channel_rows["analysis_hand"].eq("right")].copy()
    roi_rows = derive_roi_rows(right_only, load_config())
    assert len(roi_rows) == 2 * EXPECTED_ROI_DEFINITIONS
    assert set(
        roi_rows.loc[
            roi_rows["roi_name"].eq("contralateral_motor_beta"),
            "physical_source_channels",
        ]
    ) == {"FC3;C3;CP3"}


def test_contract_counts_scope_and_nonoperation_invariants() -> None:
    """Configuration pins accepted rows, participants, scopes, and non-operations."""

    config = load_config()
    assert EXPECTED_TRIAL_COUNT == 8_798
    assert EXPECTED_PARTICIPANT_COUNT == 81
    assert EXPECTED_STRICT_CLEAN_COUNT == 8_789
    assert EXPECTED_DURATION_WARNING_COUNT == 9
    assert EXPECTED_SHORT_0P5_COUNT == 47
    assert EXPECTED_SHORT_1P0_COUNT == 681
    assert EXPECTED_FILE49_COUNT == 117
    assert EXPECTED_CHANNEL_DEFINITIONS == 11
    assert EXPECTED_ROI_DEFINITIONS == 13
    assert EXPECTED_CHANNEL_ROWS == 8_798 * 30 * 11
    assert EXPECTED_ROI_ROWS == 8_798 * 13
    assert config["analysis_scopes"]["primary_assigned_condition"] == "offset_block <= 5"
    assert "imagery" in config["analysis_scopes"]["imagery_final_overt_bridge"]
    assert config["safety"]["participant_averaging"] is False
    assert config["safety"]["change_trial_eligibility"] is False
    assert config["safety"]["center_or_standardize_predictors"] is False
    assert config["safety"]["fit_models"] is False
    assert len(roi_definition_table(config)) == 13


def test_standard_1005_coordinate_contract() -> None:
    """All 30 physical channels have finite standard_1005 coordinates."""

    channels = standard_1005_channel_table()
    assert channels["physical_channel"].tolist() == list(PRIMARY_CHANNELS)
    assert np.isfinite(channels[["x_m", "y_m", "z_m"]]).all().all()


def test_atomic_parquet_publication_and_streamed_aggregation(tmp_path: Path) -> None:
    """Parquet shards and aggregate publish atomically with exact row counts."""

    first = pd.DataFrame({"key": [1, 2], "value": [1.0, 2.0]})
    second = pd.DataFrame({"key": [3], "value": [3.0]})
    first_path = tmp_path / "first.parquet"
    second_path = tmp_path / "second.parquet"
    atomic_write_parquet(first_path, first)
    atomic_write_parquet(second_path, second)
    output = tmp_path / "aggregate.parquet"
    descriptor = atomic_combine_parquet([first_path, second_path], output)
    assert descriptor["row_count"] == 3
    assert pq.read_metadata(output).num_rows == 3
    assert not list(tmp_path.glob(".*.tmp-*"))


def test_source_snapshot_detects_mutation(tmp_path: Path) -> None:
    """Source immutability evidence fails on any size/mtime change."""

    source = tmp_path / "source.npy"
    np.save(source, np.arange(3))
    before = source_snapshot([source], root=tmp_path)
    compare_source_snapshots(before, source_snapshot([source], root=tmp_path))
    source.touch()
    with pytest.raises(RuntimeError, match="source_tfr_size_or_mtime_changed"):
        compare_source_snapshots(before, source_snapshot([source], root=tmp_path))


def test_output_namespace_is_ignored_and_broad_targets_fail() -> None:
    """Production output is ignored and repository-root targets are forbidden."""

    require_ignored_output_root(REPO_ROOT, REPO_ROOT / "_Data/eeg/features_v1")
    with pytest.raises(RuntimeError, match="unsafe_feature_output_root"):
        require_ignored_output_root(REPO_ROOT, REPO_ROOT)


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="ignored production output absent")
def test_completed_production_surface_matches_all_acceptance_counts() -> None:
    """Current ignored production output has exact lineage, scope, and feature counts."""

    manifest = json.loads(PRODUCTION_MANIFEST.read_text(encoding="utf-8"))
    validation = json.loads(
        (PRODUCTION_ROOT / "validation/feature_validation.json").read_text(encoding="utf-8")
    )
    lineage = validation["lineage"]
    assert manifest["status"] == "complete"
    assert len(manifest["recordings"]) == 81
    assert pq.read_metadata(PRODUCTION_ROOT / "channel_features.parquet").num_rows == 2_903_340
    assert pq.read_metadata(PRODUCTION_ROOT / "roi_features.parquet").num_rows == 114_374
    assert lineage["accepted_eeg_trials"] == 8_798
    assert lineage["matched_frozen_behaviour_rows"] == 8_798
    assert lineage["unmatched_accepted_keys"] == 0
    assert lineage["duplicate_accepted_keys"] == 0
    assert lineage["duplicate_frozen_behaviour_keys"] == 0
    assert lineage["strict_clean_trials"] == 8_789
    assert lineage["duration_warning_trials"] == 9
    assert lineage["duration_lt_0p5_trials"] == 47
    assert lineage["duration_lt_1p0_trials"] == 681
    assert lineage["file49_trials"] == 117
    assert lineage["file54_1_present"] is False
    assert lineage["id86_present"] is False
    assert lineage["participant_4_trials"] == 48
    assert lineage["participant_60_trials"] == 34
    assert lineage["participant_70_original_handedness"] == ["a"]
    assert lineage["participant_70_analysis_hand"] == ["right"]
    assert lineage["participant_70_mapping_source"] == ["owner_recollection"]
    assert lineage["physical_missing_mt_clip"] == 0
    assert lineage["imagery_bad_figure_overlap"] == 0
    assert validation["all_channel_values_finite"]
    assert validation["all_roi_values_finite"]
    assert validation["all_roi_values_reproduce_from_channel_rows"]
    assert validation["no_31_40_hz_feature_input"]
    assert validation["no_low_high_beta_split"]
    assert validation["no_participant_averaging_authoritative_output"]
    assert validation["no_eligibility_change"]
    assert validation["source_tfr_unchanged"]

    scopes = pd.read_parquet(PRODUCTION_ROOT / "views/analysis_scope_indices.parquet")
    assert len(scopes) == 8_798
    assert int(scopes["scope_primary_blocks_1_5"].sum()) == 7_307
    assert int(scopes["scope_imagery_final_overt_bridge"].sum()) == 738
    participants = pd.read_parquet(PRODUCTION_ROOT / "summaries/participant_trial_counts.parquet")
    assert {4, 60, 70}.issubset(set(participants["behavioural_id"]))


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="ignored production output absent")
def test_completed_production_reuse_does_not_rewrite_accepted_artifacts() -> None:
    """Validation-only rerun reopens/hashes without changing accepted output files."""

    accepted = [
        PRODUCTION_MANIFEST,
        PRODUCTION_ROOT / "channel_features.parquet",
        PRODUCTION_ROOT / "roi_features.parquet",
        PRODUCTION_ROOT / "validation/feature_validation.json",
    ]
    before = [(path.stat().st_size, path.stat().st_mtime_ns) for path in accepted]
    result = subprocess.run(
        [
            str(REPO_ROOT / ".venv/bin/python"),
            str(EEG_DIR / "17_construct_trial_level_features.py"),
            "--verify-current",
        ],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    after = [(path.stat().st_size, path.stat().st_mtime_ns) for path in accepted]
    assert "PASS features_v1 unchanged reuse" in result.stdout
    assert before == after


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="ignored production output absent")
def test_completed_source_tfr_immutability_evidence_is_exact() -> None:
    """All 405 accepted TFR arrays retain identical before/after size/mtime evidence."""

    evidence = json.loads(
        (PRODUCTION_ROOT / "source_tfr_immutability.json").read_text(encoding="utf-8")
    )
    assert evidence["status"] == "unchanged"
    assert evidence["array_count"] == 405
    assert evidence["total_bytes"] == 31_445_863_440
    assert evidence["before"] == evidence["after"]
