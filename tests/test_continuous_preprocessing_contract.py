"""Synthetic contract and scientific-stage tests for continuous preprocessing."""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import mne
import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from continuous_preprocessing.contracts import (  # noqa: E402
    ALL_RETAINED_CHANNELS,
    EEG_TARGET_CHANNELS,
    MASTOID_CHANNELS,
    PRIMARY_PYPREP_CRITERIA,
    SCALP_SOURCE_CHANNELS,
    load_config,
    validate_config,
)
from continuous_preprocessing.ica import (  # noqa: E402
    apply_ica_exclusions,
    decide_component_route,
    fit_ica,
    propose_eog_components,
)
from continuous_preprocessing.source import (  # noqa: E402
    apply_channel_types,
    apply_montage,
    validate_required_channels,
)
from continuous_preprocessing.stages import (  # noqa: E402
    ObjectiveStop,
    accepted_global_bad_union,
    apply_reference,
    calculate_reference_sources,
    interpolate_global_bads,
)


CONFIG_PATH = EEG_DIR / "continuous_preprocessing_config_v1.yaml"


def config() -> dict:
    """Return a newly parsed valid tracked production configuration."""

    return load_config(CONFIG_PATH)


def synthetic_typed_montaged_raw(n_times: int = 1000) -> mne.io.RawArray:
    """Create a deterministic full 37-channel continuous RawArray."""

    rng = np.random.default_rng(20260712)
    data = rng.normal(0.0, 1e-6, (len(ALL_RETAINED_CHANNELS), n_times))
    data[-1] = 0.0
    raw = mne.io.RawArray(
        data,
        mne.create_info(list(ALL_RETAINED_CHANNELS), 250.0, ch_types="eeg"),
        verbose="ERROR",
    )
    apply_channel_types(raw)
    apply_montage(raw)
    return raw


def detection_with(primary: list[str], psd: list[str] | None = None) -> dict:
    """Build minimal criterion evidence with calls assigned to correlation."""

    criteria = {key: [] for key in PRIMARY_PYPREP_CRITERIA}
    criteria["bad_by_correlation"] = primary
    criteria["bad_by_psd"] = psd or []
    return {"criteria": criteria}


def component_scores(n_components: int, proposed: list[int]) -> dict:
    """Build complete synthetic per-component score/proposal rows."""

    rows = []
    for index in range(n_components):
        rows.append(
            {
                "component_index": index,
                "method": "infomax",
                "random_seed": 20260712,
                "estimated_rank": n_components,
                "heo_score": 0.2,
                "veo_score": 0.1,
                "maximum_absolute_eog_score": 0.2,
                "eog_support_channels": ["HEO"] if index in proposed else [],
                "automatic_proposal": index in proposed,
                "final_action": "pending_decision",
                "reason": "pending_exception_or_continuation_decision",
            }
        )
    return {
        "candidates_by_eog": {"HEO": proposed, "VEO": []},
        "component_rows": rows,
    }


def test_accepted_detector_reference_and_target_channel_contract() -> None:
    """The 30-source and retained 32-target channel surfaces stay exact."""

    cfg = config()
    assert len(SCALP_SOURCE_CHANNELS) == 30
    assert len(EEG_TARGET_CHANNELS) == 32
    assert tuple(cfg["channels"]["detector_source_channels"]) == SCALP_SOURCE_CHANNELS
    assert tuple(cfg["channels"]["reference_source_channels"]) == SCALP_SOURCE_CHANNELS
    assert tuple(cfg["channels"]["eeg_target_channels"]) == EEG_TARGET_CHANNELS
    assert not set(MASTOID_CHANNELS) & set(SCALP_SOURCE_CHANNELS)
    assert set(MASTOID_CHANNELS).issubset(EEG_TARGET_CHANNELS)


def test_required_channels_fail_on_missing_duplicate_unexpected_or_reorder() -> None:
    """Source validation fails closed on every required-channel conflict."""

    valid = list(ALL_RETAINED_CHANNELS)
    assert validate_required_channels(valid)["channel_order_matches_contract"]
    with pytest.raises(ValueError, match="missing"):
        validate_required_channels(valid[:-1])
    with pytest.raises(ValueError, match="duplicated"):
        validate_required_channels(valid[:-1] + [valid[0]])
    with pytest.raises(ValueError, match="unexpected"):
        validate_required_channels(valid[:-1] + ["EXTRA"])
    reordered = valid.copy()
    reordered[0], reordered[1] = reordered[1], reordered[0]
    with pytest.raises(ValueError, match="order"):
        validate_required_channels(reordered)


def test_published_criterion_union_excludes_psd_only_findings() -> None:
    """Modern PSD calls remain report-only and never enter the primary union."""

    result = accepted_global_bad_union(
        detection_with(["F7"], psd=["FP1", "F7"]), config()
    )
    assert result["accepted_global_bads"] == ["F7"]
    assert result["psd_report_only_findings"] == ["F7", "FP1"]
    assert result["psd_exclusive_report_only_findings"] == ["FP1"]
    assert result["psd_entered_authoritative_union"] is False


def test_25_percent_stop_uses_30_channel_denominator() -> None:
    """Seven accepted bads pass while eight trigger the >=25% objective stop."""

    seven = list(SCALP_SOURCE_CHANNELS[:7])
    result = accepted_global_bad_union(detection_with(seven), config())
    assert result["accepted_global_bad_proportion"] == pytest.approx(7 / 30)
    with pytest.raises(ObjectiveStop, match="8/30"):
        accepted_global_bad_union(
            detection_with(list(SCALP_SOURCE_CHANNELS[:8])), config()
        )


def test_reference_sources_exclude_global_bads_and_mastoids() -> None:
    """Only non-bad scalp sources estimate the reference; all EEG remain targets."""

    result = calculate_reference_sources(["F7", "FP1"], config())
    assert "F7" not in result["reference_sources"]
    assert "FP1" not in result["reference_sources"]
    assert "M1" not in result["reference_sources"]
    assert "M2" not in result["reference_sources"]
    assert result["reference_source_count"] == 28
    assert result["reference_targets"] == list(EEG_TARGET_CHANNELS)


def test_reference_is_applied_to_all_32_eeg_targets_including_mastoids() -> None:
    """M1/M2 receive the exact scalp reference without joining its source pool."""

    raw = synthetic_typed_montaged_raw()
    reference = calculate_reference_sources([], config())
    before = raw.get_data(picks=list(EEG_TARGET_CHANNELS))
    expected_reference = before[:17].sum(axis=0)
    # Use the exact accepted source order rather than relying on EEG order,
    # which contains M1/M2 interleaved with the scalp channels.
    expected_reference = raw.get_data(picks=list(SCALP_SOURCE_CHANNELS)).mean(axis=0)
    m1_before = raw.get_data(picks=["M1"])[0]
    evidence = apply_reference(raw, reference)
    m1_after = raw.get_data(picks=["M1"])[0]
    assert evidence["all_32_eeg_targets_received_reference"]
    assert np.allclose(m1_after, m1_before - expected_reference)


def test_interpolation_planning_requires_scalp_candidates_and_resets_bads() -> None:
    """M1/M2 cannot be candidates and a valid scalp interpolation resets bads."""

    raw = synthetic_typed_montaged_raw(2500)
    reference = calculate_reference_sources(["F7"], config())
    apply_reference(raw, reference)
    with pytest.raises(ValueError, match="M1/M2"):
        interpolate_global_bads(raw, ["M1"], config())
    result = interpolate_global_bads(raw, ["F7"], config())
    assert result["interpolation_count"] == 1
    assert result["interpolation_proportion"] == pytest.approx(1 / 30)
    assert result["operational_bad_list_after"] == []
    assert raw.info["bads"] == []


def test_line_filter_and_ica_copy_contracts_are_fixed() -> None:
    """The run toggle, analysis FIR, ICA copy, method, seed, and rank contract persist."""

    cfg = config()
    assert cfg["line_noise"]["enabled"] is False
    enabled = copy.deepcopy(cfg)
    enabled.pop("_provenance")
    enabled["line_noise"]["enabled"] = True
    validate_config(enabled)
    expected_filter = {
        "high_pass_hz": 0.5,
        "low_pass_hz": 45.0,
        "method": "fir",
        "phase": "zero",
        "fir_design": "firwin",
    }
    assert all(cfg["analysis_filter"][key] == value for key, value in expected_filter.items())
    assert cfg["ica"]["high_pass_hz"] == 1.0
    assert cfg["ica"]["low_pass_hz"] == 45.0
    assert cfg["ica"]["method"] == "infomax"
    assert cfg["ica"]["extended"] is True
    assert cfg["ica"]["random_seed"] == 20260712
    assert cfg["ica"]["n_components"] == "estimated_eeg_rank"


def test_fit_ica_passes_rank_method_seed_and_extended_contract(monkeypatch: pytest.MonkeyPatch) -> None:
    """The fitted component count comes from rank, never a historical fixed count."""

    captured: dict = {}

    class FakeICA:
        def __init__(self, **kwargs):
            captured["init"] = kwargs
            self.n_components_ = kwargs["n_components"]
            self.n_iter_ = 5
            self.current_fit = "raw"
            self.ch_names = []

        def fit(self, raw, **kwargs):
            captured["fit"] = kwargs
            self.ch_names = list(kwargs["picks"])
            return self

    monkeypatch.setattr(mne.preprocessing, "ICA", FakeICA)
    raw = synthetic_typed_montaged_raw()
    fitted, evidence = fit_ica(raw, {"estimated_eeg_rank": 27}, config())
    assert fitted.n_components_ == 27
    assert captured["init"]["method"] == "infomax"
    assert captured["init"]["fit_params"] == {"extended": True}
    assert captured["init"]["random_state"] == 20260712
    assert captured["fit"]["picks"] == list(EEG_TARGET_CHANNELS)
    assert evidence["fit_decim"] == config()["ica"]["fit_decim"]


def test_ordinary_zero_more_than_two_and_id86_component_routes() -> None:
    """Zero continues, >2 stops, and ID 86 always reaches its review boundary."""

    zero_scores = component_scores(5, [])
    zero_proposal = propose_eog_components(zero_scores, config())
    zero = decide_component_route(1, zero_scores, zero_proposal, config())
    assert zero["route"] == "continue_automatic"
    assert zero["final_exclusions"] == []
    assert zero["reason"] == "ordinary_zero_component_path"

    many_scores = component_scores(5, [0, 1, 2])
    many_proposal = propose_eog_components(many_scores, config())
    many = decide_component_route(1, many_scores, many_proposal, config())
    assert many["route"] == "stop_for_exceptional_review"
    assert many["final_exclusions"] == []

    id86_scores = component_scores(5, [0])
    id86_proposal = propose_eog_components(id86_scores, config())
    id86 = decide_component_route(86, id86_scores, id86_proposal, config())
    assert id86["route"] == "stop_for_component_review"
    assert id86["automatic_application_authorized"] is False
    assert id86["final_exclusions"] == []


def test_zero_component_application_does_not_reconstruct_or_change_samples() -> None:
    """The ordinary zero-component path records exclusions without ICA.apply."""

    class FakeICA:
        exclude: list[int] = []
        applied = False

        def apply(self, raw, exclude, verbose):
            self.applied = True

    raw = synthetic_typed_montaged_raw()
    before = raw.get_data().copy()
    fake = FakeICA()
    route = {
        "automatic_application_authorized": True,
        "final_exclusions": [],
        "reason": "ordinary_zero_component_path",
    }
    result = apply_ica_exclusions(fake, raw, route)
    assert result["zero_component_branch_left_analysis_samples_unchanged"]
    assert fake.applied is False
    assert np.array_equal(raw.get_data(), before)


def test_production_namespace_contains_no_epoch_construction_api() -> None:
    """Static production code has no MNE Epochs constructor or AutoReject call."""

    production_files = list((EEG_DIR / "continuous_preprocessing").glob("*.py")) + [
        EEG_DIR / "13_run_continuous_preprocessing_validation.py"
    ]
    text = "\n".join(path.read_text(encoding="utf-8") for path in production_files)
    assert "mne.Epochs(" not in text
    assert "EpochsArray(" not in text
    assert "AutoReject(" not in text
    assert "compute_current_source_density(" not in text
