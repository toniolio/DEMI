"""Focused public tests for event/epoch eligibility-ledger contracts.

These tests use synthetic rows and manifests.  They do not read private policy
records, EDF/FIF signals, or ignored generated evidence.  The numbered ledger
driver adds a separate optional local integration surface after implementation.
"""

from __future__ import annotations

import sys
import importlib.util
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from event_epoch_eligibility import (  # noqa: E402
    assert_no_signal_outputs,
    canonical_event_key,
    parse_continuous_recording_manifest,
    stable_frame_hash,
    status_fields,
    validate_ledger_contract,
)


LEDGER_SCRIPT_PATH = EEG_DIR / "14_build_event_epoch_eligibility_ledger.py"
LOCAL_POLICY_PATH = REPO_ROOT / "_Private/inventory/eeg_event_policy_v1.0.csv"
LOCAL_JOIN_PATH = (
    REPO_ROOT / "_Data/eeg/event_evidence_v1/proposed_offset_join_audit.csv"
)


def load_ledger_driver():
    """Import numbered script 14 by path for optional local integration tests."""

    spec = importlib.util.spec_from_file_location(
        "event_epoch_eligibility_driver_test_module", LEDGER_SCRIPT_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def local_ledger_result():
    """Build the local accepted-authority ledger once when ignored inputs exist."""

    if not LOCAL_POLICY_PATH.is_file() or not LOCAL_JOIN_PATH.is_file():
        pytest.skip("ignored accepted policy/event evidence is not available")
    driver = load_ledger_driver()
    return driver.build_ledger(REPO_ROOT)


def test_canonical_keys_distinguish_raw_and_offset_only_rows() -> None:
    """Native source/sample identity and frozen offset identity stay distinct."""

    raw = canonical_event_key(
        {
            "source_filename": "demi_13 Data.edf",
            "eeg_recording_id": 13,
            "behavioural_participant_id": 11,
            "raw_trial_sequence": 3,
            "anchor_sample": 123_456,
        }
    )
    offset = canonical_event_key(
        {
            "source_filename": "",
            "eeg_recording_id": None,
            "behavioural_participant_id": 6,
            "offset_trial_count": 1,
        }
    )
    assert raw != offset
    assert raw == canonical_event_key(
        {
            "source_filename": "demi_13 Data.edf",
            "eeg_recording_id": 13,
            "behavioural_participant_id": 11,
            "raw_trial_sequence": 3,
            "anchor_sample": 123_456,
        }
    )


def test_status_fields_always_explain_false_state() -> None:
    """Unavailable statuses always carry machine and readable reasons."""

    fields = status_fields(
        "post_ica_derivative_availability",
        False,
        "post_ica_available",
        "Post-ICA derivative exists.",
        "accepted_id86_stop_no_post_ica",
        "The accepted ID-86 route retains pre-ICA evidence only.",
    )
    assert fields["post_ica_derivative_availability"] is False
    assert fields["post_ica_derivative_availability_status"] == "unavailable"
    assert fields["post_ica_derivative_availability_reason_code"]
    assert fields["post_ica_derivative_availability_reason"]


@pytest.mark.parametrize(
    ("artifact_kind", "status", "stop_code", "expected_kind", "expected_route"),
    [
        (
            "post_ica_continuous",
            "complete",
            "",
            "post_ica",
            "ordinary_post_ica_complete",
        ),
        (
            "pre_ica_continuous",
            "stopped",
            "predeclared_id86_component_review_boundary",
            "retained_pre_ica",
            "accepted_historical_zero_component_stop_pre_ica_retained",
        ),
    ],
)
def test_continuous_manifest_parser_preserves_post_vs_pre_ica_route(
    tmp_path: Path,
    artifact_kind: str,
    status: str,
    stop_code: str,
    expected_kind: str,
    expected_route: str,
) -> None:
    """Ordinary completion and accepted ID-86 stop cannot collapse together."""

    recording_dir = tmp_path / "recordings" / "demi_86_data"
    manifest_path = recording_dir / "manifest.json"
    artifact_name = "post_ica_raw.fif" if expected_kind == "post_ica" else "pre_ica_raw.fif"
    manifest = {
        "manifest_id": "manifest-id",
        "status": status,
        "completion_class": status,
        "source": {"source_filename": "demi_86 Data.edf", "source_sha256": "abc"},
        "artifacts": [{"kind": artifact_kind, "relative_path": artifact_name}],
        "qc_warnings": [],
        "interpolation": {
            "interpolation_count": 4,
            "interpolation_proportion": 4 / 30,
            "proportion_denominator": 30,
        },
        "ica": {
            "proposal": {"proposal_count": 0},
            "application": {
                "exclusion_count": 0,
                "ica_application_performed": False,
            },
        },
        "stop_or_failure": {
            "code": stop_code,
            "detail": "Retained evidence.",
        },
    }
    parsed = parse_continuous_recording_manifest(
        manifest, manifest_path, tmp_path
    )
    assert parsed["continuous_derivative_kind"] == expected_kind
    assert parsed["ica_terminal_route"] == expected_route
    assert parsed["continuous_derivative_path"].endswith(artifact_name)


def test_signal_output_guard_rejects_epoch_fif(tmp_path: Path) -> None:
    """The ledger namespace must never contain epochs or signal derivatives."""

    (tmp_path / "unexpected-epo.fif").write_bytes(b"not a real FIF")
    with pytest.raises(RuntimeError, match="forbidden signal outputs"):
        assert_no_signal_outputs(tmp_path)


def test_signal_output_guard_accepts_ledger_tables(tmp_path: Path) -> None:
    """CSV, Parquet, JSON, and Markdown review artifacts remain allowed."""

    pd.DataFrame([{"a": 1}]).to_csv(tmp_path / "ledger.csv", index=False)
    (tmp_path / "summary.json").write_text("{}", encoding="utf-8")
    (tmp_path / "summary.md").write_text("# Summary\n", encoding="utf-8")
    assert_no_signal_outputs(tmp_path)


def test_local_accepted_counts_and_warning_only_difference(local_ledger_result) -> None:
    """Local accepted authority reproduces 8,905, 8,896, and exactly nine."""

    ledger = local_ledger_result["ledger"]
    counts = validate_ledger_contract(ledger)
    assert counts == {
        "primary": 8_905,
        "strict_clean": 8_896,
        "warning_only": 9,
        "possible_discrepancy": 72,
    }
    warning_only = ledger[
        ledger["primary_ledger_eligibility"]
        & ~ledger["strict_clean_only_sensitivity_eligibility"]
    ]
    assert len(warning_only) == 9
    assert warning_only["duration_warning_only"].all()


def test_local_identity_and_physical_trigger_contracts(local_ledger_result) -> None:
    """Behaviour 11 uses EEG 13; EEG 11 cannot join; EEG 94 stays physical."""

    ledger = local_ledger_result["ledger"]
    primary = ledger[ledger["primary_ledger_eligibility"]]
    mapped11 = primary[primary["behavioural_participant_id"].eq(11)]
    assert len(mapped11) == 99
    assert set(mapped11["eeg_recording_id"].astype(int)) == {13}
    assert not ledger["eeg_recording_id"].eq(11).any()
    source11 = local_ledger_result["source_status"].query("eeg_recording_id == 11")
    assert len(source11) == 1
    assert not source11.iloc[0]["behavioural_join_available"]
    assert source11.iloc[0]["behavioural_join_reason_code"] == (
        "eeg11_explicitly_unmapped_from_behavioural_join"
    )
    eeg94 = primary[primary["eeg_recording_id"].eq(94)]
    assert len(eeg94) == 104
    assert set(eeg94["source_type"]) == {"physical_trigger"}


def test_local_split_id5_and_discrepancy_policy_contracts(local_ledger_result) -> None:
    """Split IDs stay out; ID 5 stays pre-marker; all 72 remain unavailable."""

    ledger = local_ledger_result["ledger"]
    primary = ledger[ledger["primary_ledger_eligibility"]]
    assert not primary["eeg_recording_id"].isin({54, 56, 65}).any()
    id5 = primary[primary["behavioural_participant_id"].eq(5)]
    assert len(id5) == 114
    assert set(id5["id5_segment"]) == {"before_file_start"}
    discrepancy = ledger[
        ledger["event_policy_reason_code"].eq(
            "possible_event_sequence_discrepancy_unavailable"
        )
    ]
    assert len(discrepancy) == 72
    assert not discrepancy["primary_ledger_eligibility"].any()


def test_local_continuous_warning_and_id86_contracts(local_ledger_result) -> None:
    """49/54_1 warnings propagate without exclusion; ID 86 remains pre-ICA."""

    ledger = local_ledger_result["ledger"]
    for filename in ("demi_49 Data.edf", "demi_54_1 Data.edf"):
        rows = ledger[ledger["source_filename"].eq(filename)]
        assert not rows.empty
        assert rows["continuous_qc_warning_codes"].str.contains(
            "global_bad_proportion_above_25_percent"
        ).all()
        assert set(rows["interpolated_channel_count"].astype(int)) == {8}
    file49_primary = ledger[
        ledger["source_filename"].eq("demi_49 Data.edf")
        & ledger["primary_ledger_eligibility"]
    ]
    assert file49_primary["future_epoch_construction_readiness"].all()
    id86 = ledger[ledger["eeg_recording_id"].eq(86)]
    assert set(id86["continuous_derivative_kind"]) == {"retained_pre_ica"}
    assert id86["continuous_derivative_availability"].all()
    assert not id86["post_ica_derivative_availability"].any()
    assert not id86["future_epoch_construction_readiness"].any()
    assert set(id86["ica_exclusion_count"].astype(int)) == {0}


def test_local_keys_reasons_and_deterministic_rebuild(local_ledger_result) -> None:
    """Keys are unique, false statuses are explained, and a rebuild is stable."""

    ledger = local_ledger_result["ledger"]
    assert not ledger["canonical_event_key"].duplicated().any()
    for dimension in (
        "accepted_event_policy_candidate",
        "primary_ledger_eligibility",
        "strict_clean_only_sensitivity_eligibility",
        "continuous_derivative_availability",
        "post_ica_derivative_availability",
        "future_epoch_construction_readiness",
    ):
        unavailable = ledger[~ledger[dimension]]
        assert unavailable[f"{dimension}_reason_code"].fillna("").str.strip().ne("").all()
        assert unavailable[f"{dimension}_reason"].fillna("").str.strip().ne("").all()
    driver = load_ledger_driver()
    rebuilt = driver.build_ledger(REPO_ROOT)["ledger"]
    assert stable_frame_hash(rebuilt) == stable_frame_hash(ledger)


def test_local_output_namespace_contains_no_epoch_or_signal_derivative(
    local_ledger_result,
) -> None:
    """Generated ledger output remains tabular/review-only."""

    del local_ledger_result
    assert_no_signal_outputs(REPO_ROOT / "_Data/eeg/event_epoch_eligibility_v1")
