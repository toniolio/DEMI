"""Focused tests for the accepted DEMI EEG epoch-construction contract."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from epoch_construction import (  # noqa: E402
    EXPECTED_DURATION_WARNING_COUNT,
    EXPECTED_EPOCH_COUNT,
    EXPECTED_FILE49_COUNT,
    EXPECTED_STRICT_CLEAN_COUNT,
    FAMILIES,
    FAMILY_BY_NAME,
    assess_cached_epoch_shard,
    atomic_write_json,
    assert_identical_family_keys,
    build_family_metadata,
    compare_source_snapshots,
    create_epochs,
    derive_anchor,
    expected_epoch_sample_count,
    forbidden_product_scan,
    reopen_and_validate_epochs,
    seconds_to_sample,
    select_epoch_rows,
    source_snapshot,
    write_epochs_atomic,
)


LOCAL_LEDGER = (
    REPO_ROOT
    / "_Data"
    / "eeg"
    / "event_epoch_eligibility_v1"
    / "event_epoch_eligibility_ledger.parquet"
)


def event(role: str, onset: float, sample: int) -> dict:
    """Return one minimal serialized accepted event-provenance record."""

    codes = {"stim_on": 28, "red_on": 30, "trace_start": 44, "trace_end": 46}
    return {
        "onset_seconds": onset,
        "sample": sample,
        "source_type": "annotation",
        "raw_event_name": role,
        "raw_value": str(codes[role]),
        "raw_code": codes[role],
        "raw_annotation_description": str(codes[role]),
        "proposed_event_name": role,
        "proposed_sequence_action": "unchanged",
        "proposed_sequence_reason": "raw_label_preserved",
    }


def synthetic_row(*, physical: bool = True, key: str = "key-1") -> dict:
    """Return one future-ready row carrying all metadata needed by helpers."""

    provenance = [
        event("stim_on", 1.0, 1000),
        event("red_on", 2.0, 2000),
        event("trace_start", 3.0, 3000),
        event("trace_end", 4.0, 4000),
    ]
    return {
        "canonical_order_index": 0,
        "canonical_event_key": key,
        "ledger_version": "event_epoch_eligibility_v1",
        "event_policy_version": "v1.0",
        "behavioural_participant_id": 1,
        "eeg_recording_id": 1,
        "source_filename": "demi_01 Data.edf",
        "file_role": "single",
        "source_type": "annotation",
        "selected_source": "annotation",
        "task_filename": "p1.example.txt",
        "audit_trial_count": 1,
        "raw_trial_sequence": 3,
        "join_trial_count": 1,
        "offset_trial_count": 1,
        "physical": physical,
        "trace_start_onset_seconds": 3.0,
        "trace_end_onset_seconds": 4.0,
        "red_on_onset_seconds": 2.0,
        "real_start": 0.25 if physical else 0.0,
        "real_end": 1.0,
        "trial_event_provenance_json": json.dumps(provenance),
        "primary_ledger_eligibility": True,
        "primary_ledger_eligibility_reason_code": "accepted_primary_candidate",
        "strict_clean_only_sensitivity_eligibility": True,
        "strict_clean_only_sensitivity_eligibility_reason_code": "strict_clean",
        "duration_warning_only": False,
        "event_policy_reason_code": "accepted_matched_event",
        "accepted_policy_ids": "ordinary",
        "accepted_specific_policy_id": "ordinary",
        "accepted_policy_action": "include",
        "continuous_v2_run_id": "run",
        "continuous_v2_manifest_path": "continuous/manifest.json",
        "continuous_v2_manifest_id": "manifest",
        "continuous_derivative_kind": "post_ica",
        "continuous_derivative_path": "continuous/post_ica_raw.fif",
        "continuous_source_sha256": "source",
        "continuous_qc_warning_codes": "",
        "interpolated_channel_count": 0,
        "interpolated_channel_proportion": 0.0,
        "interpolation_proportion_denominator": 30,
        "ica_terminal_route": "ordinary_post_ica_complete",
        "post_ica_derivative_availability": True,
        "post_ica_derivative_availability_reason_code": "available",
        "future_epoch_construction_readiness": True,
        "future_epoch_construction_readiness_reason_code": "ready",
    }


def synthetic_raw() -> mne.io.RawArray:
    """Return a small native-1000-Hz Raw with no annotations or projections."""

    rng = np.random.default_rng(20260716)
    data = rng.normal(scale=1e-6, size=(2, 10_000))
    info = mne.create_info(["C3", "C4"], sfreq=1000.0, ch_types="eeg")
    return mne.io.RawArray(data, info, verbose="ERROR")


def test_accepted_anchor_formulas_for_overt_and_imagery() -> None:
    """All four response formula branches and red_on use accepted fields."""

    overt = synthetic_row(physical=True)
    imagery = synthetic_row(physical=False)

    overt_onset = derive_anchor(overt, "response_onset")
    imagery_onset = derive_anchor(imagery, "response_onset")
    overt_end = derive_anchor(overt, "response_end")
    imagery_end = derive_anchor(imagery, "response_end")
    red_on = derive_anchor(overt, "red_on")

    assert overt_onset["anchor_time_seconds"] == 3.25
    assert overt_onset["anchor_source"] == "trace_start_plus_real_start"
    assert imagery_onset["anchor_time_seconds"] == 3.0
    assert imagery_onset["anchor_source"] == "raw_trace_start"
    assert overt_end["anchor_time_seconds"] == 4.0
    assert overt_end["anchor_source"] == "trace_start_plus_real_end"
    assert imagery_end["anchor_time_seconds"] == 4.0
    assert imagery_end["anchor_source"] == "raw_trace_end"
    assert red_on["anchor_time_seconds"] == 2.0
    assert red_on["anchor_source"] == "raw_red_on"


def test_deterministic_sample_conversion_and_half_sample_fail_closed() -> None:
    """Nearest-sample conversion is stable and exact half ties are rejected."""

    assert seconds_to_sample(3.2504, 1000.0) == 3250
    assert seconds_to_sample(3.2506, 1000.0) == 3251
    assert seconds_to_sample(3.2504, 1000.0, first_samp=50) == 3300
    with pytest.raises(ValueError, match="ambiguous half-sample"):
        seconds_to_sample(3.2505, 1000.0)


def test_inclusive_endpoint_sample_counts() -> None:
    """MNE includes both exact -1.5/+endpoint samples at native 1000 Hz."""

    assert expected_epoch_sample_count(-1.5, 2.5, 1000.0) == 4001
    assert expected_epoch_sample_count(-1.5, 0.8, 1000.0) == 2301


@pytest.mark.parametrize("family", FAMILIES, ids=lambda item: item.name)
def test_mne_epochs_are_preloaded_native_unbaselined_and_unrejected(
    family,
    tmp_path: Path,
) -> None:
    """Saved and reopened family shards preserve the signal-neutral contract."""

    raw = synthetic_raw()
    metadata = build_family_metadata(
        pd.DataFrame([synthetic_row()]),
        family,
        sfreq=raw.info["sfreq"],
        first_samp=raw.first_samp,
        last_samp=raw.last_samp,
    )
    epochs = create_epochs(raw, metadata, family)
    assert epochs.preload
    assert epochs.info["sfreq"] == 1000.0
    assert epochs.baseline is None
    assert epochs.reject is None
    assert epochs.flat is None
    assert epochs.detrend is None
    assert len(epochs.times) == expected_epoch_sample_count(
        family.tmin, family.tmax, 1000.0
    )
    assert epochs.times[0] == family.tmin
    assert epochs.times[-1] == family.tmax
    assert epochs.metadata["canonical_event_key"].tolist() == ["key-1"]

    path = tmp_path / f"{family.name}-epo.fif"
    written = write_epochs_atomic(epochs, path)
    validation = reopen_and_validate_epochs(
        path, metadata, family, expected_sha256=written["sha256"]
    )
    assert validation["valid"]
    assert validation["preloaded"]
    assert validation["drop_count"] == 0
    assert validation["baseline"] is None


def test_out_of_bounds_epoch_stops_construction() -> None:
    """An accepted row cannot be silently dropped at a recording boundary."""

    row = synthetic_row()
    row["red_on_onset_seconds"] = 0.5
    provenance = json.loads(row["trial_event_provenance_json"])
    provenance[1]["onset_seconds"] = 0.5
    provenance[1]["sample"] = 500
    row["trial_event_provenance_json"] = json.dumps(provenance)
    with pytest.raises(RuntimeError, match="out-of-bounds"):
        build_family_metadata(
            pd.DataFrame([row]),
            FAMILY_BY_NAME["red_on"],
            sfreq=1000.0,
            first_samp=0,
            last_samp=9999,
        )


def test_identical_ordered_keys_across_families() -> None:
    """Family metadata must carry exactly the same canonical row ordering."""

    raw = synthetic_raw()
    frames = {
        family.name: build_family_metadata(
            pd.DataFrame([synthetic_row()]),
            family,
            sfreq=1000.0,
            first_samp=raw.first_samp,
            last_samp=raw.last_samp,
        )
        for family in FAMILIES
    }
    digest = assert_identical_family_keys(frames)
    assert len(digest) == 64
    frames["response_end"].loc[0, "canonical_event_key"] = "different"
    with pytest.raises(RuntimeError, match="ordered canonical keys differ"):
        assert_identical_family_keys(frames)


def test_source_snapshot_detects_mutation(tmp_path: Path) -> None:
    """Source immutability evidence fails when size or mtime changes."""

    source = tmp_path / "post_ica_raw.fif"
    source.write_bytes(b"unchanged")
    before = source_snapshot([source])
    assert compare_source_snapshots(before, source_snapshot([source]))["unchanged"]
    source.write_bytes(b"changed-size")
    with pytest.raises(RuntimeError, match="changed"):
        compare_source_snapshots(before, source_snapshot([source]))


def test_terminal_shard_cache_is_deterministic_and_fail_closed(tmp_path: Path) -> None:
    """Only a complete matching hash-valid shard is reusable on rerun."""

    artifact = tmp_path / "response_onset-epo.fif"
    artifact.write_bytes(b"epoch artifact")
    from epoch_construction import sha256_file

    manifest = tmp_path / "manifest.json"
    atomic_write_json(
        manifest,
        {
            "status": "complete",
            "fingerprint": "current",
            "artifact": {
                "size_bytes": artifact.stat().st_size,
                "sha256": sha256_file(artifact),
            },
        },
    )
    assert assess_cached_epoch_shard(manifest, artifact, "current")["action"] == "reuse"
    assert assess_cached_epoch_shard(manifest, artifact, "changed") == {
        "action": "rebuild",
        "reason": "provenance_fingerprint_changed",
    }
    artifact.write_bytes(b"changed")
    assert assess_cached_epoch_shard(manifest, artifact, "current")["action"] == "rebuild"


def test_forbidden_product_scan_catches_analysis_scope_expansion(tmp_path: Path) -> None:
    """The epoch namespace rejects TFR/CSD/AutoReject/pre-ICA products."""

    (tmp_path / "response_onset-epo.fif").touch()
    assert forbidden_product_scan(tmp_path) == []
    (tmp_path / "response_onset-tfr.h5").touch()
    (tmp_path / "id86_pre_ica-epo.fif").touch()
    assert forbidden_product_scan(tmp_path) == [
        "id86_pre_ica-epo.fif",
        "response_onset-tfr.h5",
    ]


@pytest.fixture(scope="module")
def local_selected_rows() -> pd.DataFrame:
    """Load and validate the generated accepted ledger when locally present."""

    if not LOCAL_LEDGER.is_file():
        pytest.skip("local accepted event/epoch ledger is not available")
    return select_epoch_rows(pd.read_parquet(LOCAL_LEDGER))


def test_local_surface_exact_counts_and_special_routes(local_selected_rows: pd.DataFrame) -> None:
    """Local authority gives 8,798/8,789/nine, file 49, no 54_1, and no 86."""

    selected = local_selected_rows
    assert len(selected) == EXPECTED_EPOCH_COUNT
    assert int(selected["strict_clean_only_sensitivity_eligibility"].sum()) == (
        EXPECTED_STRICT_CLEAN_COUNT
    )
    assert int(selected["duration_warning_only"].sum()) == (
        EXPECTED_DURATION_WARNING_COUNT
    )
    assert len(selected[selected["eeg_recording_id"].eq(49)]) == EXPECTED_FILE49_COUNT
    assert not selected["source_filename"].str.contains("demi_54_1", regex=False).any()
    assert not selected["eeg_recording_id"].eq(86).any()


def test_local_all_three_anchors_available_and_family_keys_identical(
    local_selected_rows: pd.DataFrame,
) -> None:
    """Every selected local row yields all anchors in canonical order."""

    frames: dict[str, pd.DataFrame] = {}
    for family in FAMILIES:
        shards = []
        for _, rows in local_selected_rows.groupby(
            "continuous_derivative_path", sort=False
        ):
            shards.append(
                build_family_metadata(
                    rows,
                    family,
                    sfreq=1000.0,
                    first_samp=0,
                    last_samp=10_000_000,
                )
            )
        frames[family.name] = pd.concat(shards, ignore_index=True)
        assert len(frames[family.name]) == EXPECTED_EPOCH_COUNT
        assert frames[family.name]["anchor_time_seconds"].notna().all()
    assert_identical_family_keys(frames)
