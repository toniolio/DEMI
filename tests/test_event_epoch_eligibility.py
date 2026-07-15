"""Focused public tests for event/epoch eligibility-ledger contracts.

These tests use synthetic rows and manifests.  They do not read private policy
records, EDF/FIF signals, or ignored generated evidence.  The numbered ledger
driver adds a separate optional local integration surface after implementation.
"""

from __future__ import annotations

import sys
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
    status_fields,
)


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
