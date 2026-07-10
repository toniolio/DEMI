"""Focused tests for DEMI participant identity and event-source selection.

The tests use small synthetic event streams and the tracked minimal identity
contract. They do not open private EDF, task, TraceLab, masterlist, or RDS
files. One optional local regression test checks the frozen old-compatible
offset CSV when that ignored file is available.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_ANALYSIS_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_ANALYSIS_DIR))

from event_source_contract import (  # noqa: E402
    IdentityMapping,
    compare_event_streams,
    load_identity_contract,
    require_identity_mapping,
    selected_event_rows,
    select_event_source,
    validate_identity_coverage,
)


IDENTITY_CONTRACT_PATH = EEG_ANALYSIS_DIR / "eeg_behavior_identity_contract.csv"
CORRECTED_EVIDENCE_SCRIPT_PATH = (
    EEG_ANALYSIS_DIR / "06_build_corrected_event_evidence.py"
)
OLD_COMPATIBLE_OFFSETS_PATH = (
    REPO_ROOT
    / "_Data"
    / "behavior"
    / "event_offsets"
    / "event_offsets_old_compatible.csv"
)

# This digest protects the exact local old-compatible target without requiring
# ignored behavioural data in public CI. The test is skipped when the local
# target is unavailable.
OLD_COMPATIBLE_OFFSETS_SHA256 = (
    "2a8b1e2bedfdce2c5853cee44624cb4beddc126e149f03137f1d94e60ac9cded"
)


def load_corrected_evidence_module():
    """Import script 06 by path because its filename starts with a number."""

    spec = importlib.util.spec_from_file_location(
        "corrected_event_evidence_test_module",
        CORRECTED_EVIDENCE_SCRIPT_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


corrected_evidence = load_corrected_evidence_module()


def event(
    sample: int,
    raw_code: int,
    normalized_candidate_code: int | None = None,
    source_type: str = "annotation",
) -> dict[str, object]:
    """Build one small event row for pure source-contract tests."""

    if normalized_candidate_code is None:
        normalized_candidate_code = raw_code
    return {
        "sample": sample,
        "raw_code": raw_code,
        "normalized_candidate_code": normalized_candidate_code,
        "source_type": source_type,
    }


def coherent_stream(
    *,
    source_type: str,
    raw_codes: tuple[int, ...] = (28, 30, 44, 46, 60),
    candidate_codes: tuple[int, ...] | None = None,
    start_sample: int = 100,
) -> list[dict[str, object]]:
    """Build two balanced synthetic DEMI trials."""

    if candidate_codes is None:
        candidate_codes = raw_codes
    rows: list[dict[str, object]] = []
    for trial in range(2):
        for position, (raw_code, candidate_code) in enumerate(
            zip(raw_codes, candidate_codes)
        ):
            rows.append(
                event(
                    start_sample + trial * 100 + position * 10,
                    raw_code,
                    candidate_code,
                    source_type,
                )
            )
    return rows


def selected_inventory_event(
    *,
    eeg_recording_id: int = 13,
    behavioural_participant_id: int | None = 11,
) -> dict[str, object]:
    """Return one selected-source row in the version-1 inventory schema."""

    return {
        "eeg_recording_id": eeg_recording_id,
        "eeg_recording_id_padded": f"{eeg_recording_id:03d}",
        "behavioural_participant_id": behavioural_participant_id,
        "identity_mapping_status": (
            "explicit_remap"
            if eeg_recording_id == 13
            else "explicitly_unmapped"
        ),
        "identity_contract_version": "v1",
        "source_filename": f"demi_{eeg_recording_id:02d} Data.edf",
        "file_path": f"_Data/eeg/raw/demi_{eeg_recording_id:02d} Data.edf",
        "file_role": "single",
        "split_part": None,
        "session_date": "2018-06-05",
        "sampling_frequency_hz": 1000.0,
        "source_type": "annotation",
        "event_row_order": 1,
        "sample": 100,
        "onset_seconds": 0.1,
        "duration_seconds": 0.0,
        "raw_value": "28",
        "previous_raw_value": "",
        "raw_code": 28,
        "previous_code": None,
        "normalized_candidate_code": 28,
        "normalization_candidate_rule": "annotation_label_preserved",
        "normalization_candidate_applied": False,
        "candidate_event_name": "stim_on",
        "task_event_candidate": True,
        "source_extraction_parameters_json": "{}",
        "selected_for_event_evidence": True,
        "selected_source": "annotation",
        "selection_status": "annotation_selected",
        "selection_reason": "coherent annotations preferred by contract",
        "selection_unresolved_issue": "",
        "normalization_applied": False,
        "selected_normalized_code": 28,
        "selected_event_name": "stim_on",
    }


def source_file_inventory(
    *,
    eeg_recording_id: int = 13,
    behavioural_participant_id: int | None = 11,
) -> dict[str, object]:
    """Return one compact file row used by corrected-evidence tests."""

    return {
        "source_filename": f"demi_{eeg_recording_id:02d} Data.edf",
        "recording_duration_seconds": 10.0,
        "source_discrepancy_class": "coherent_exact_agreement",
        "selected_source": "annotation",
        "selection_status": "annotation_selected",
        "selection_reason": "coherent annotations preferred by contract",
        "selection_unresolved_issue": "",
        "annotation_event_count": 10,
        "physical_event_count": 10,
        "eeg_recording_id": eeg_recording_id,
        "behavioural_participant_id": behavioural_participant_id,
    }


def test_public_identity_contract_encodes_reviewed_11_to_13_mapping() -> None:
    """Behaviour 11 must be reached through EEG 13, never EEG 11."""

    mappings = load_identity_contract(IDENTITY_CONTRACT_PATH)

    eeg13 = require_identity_mapping(mappings, 13)
    eeg11 = require_identity_mapping(mappings, 11)

    assert eeg13.behavioural_participant_id == 11
    assert eeg13.mapping_status == "explicit_remap"
    assert eeg11.behavioural_participant_id is None
    assert eeg11.mapping_status == "explicitly_unmapped"
    assert [
        mapping.eeg_recording_id
        for mapping in mappings.values()
        if mapping.behavioural_participant_id == 11
    ] == [13]


def test_missing_identity_never_falls_back_to_same_number() -> None:
    """An absent mapping must fail instead of silently using numeric equality."""

    mappings = {
        13: IdentityMapping(13, 11, "explicit_remap", "v1"),
    }
    with pytest.raises(
        RuntimeError, match="absent from the explicit identity contract"
    ):
        require_identity_mapping(mappings, 11)
    with pytest.raises(RuntimeError, match="lacks discovered EEG recording"):
        validate_identity_coverage(mappings, [11, 13])


def test_corrected_alignment_uses_behaviour_mapping_not_eeg_number() -> None:
    """Selected EEG 13 events must enter the join as behavioural ID 11."""

    selected = pd.DataFrame([selected_inventory_event()])
    files = pd.DataFrame([source_file_inventory()])

    aligned = corrected_evidence.build_selected_alignment_events(
        selected, files
    )

    assert aligned.loc[0, "eeg_recording_id"] == 13
    assert aligned.loc[0, "participant_id"] == 11
    assert aligned.loc[0, "behavioural_participant_id"] == 11
    assert aligned.loc[0, "identity_mapping_status"] == "explicit_remap"


def test_corrected_alignment_rejects_selected_event_without_mapping() -> None:
    """An explicitly unmapped EEG event cannot enter behavioural alignment."""

    selected = pd.DataFrame(
        [
            selected_inventory_event(
                eeg_recording_id=11,
                behavioural_participant_id=None,
            )
        ]
    )
    files = pd.DataFrame(
        [
            source_file_inventory(
                eeg_recording_id=11,
                behavioural_participant_id=None,
            )
        ]
    )
    with pytest.raises(RuntimeError, match="without a behavioural mapping"):
        corrected_evidence.build_selected_alignment_events(selected, files)


def test_coherent_annotation_and_physical_streams_select_annotations() -> None:
    """Exact coherent streams should retain annotations as the preferred source."""

    annotations = coherent_stream(source_type="annotation")
    physical = coherent_stream(source_type="physical_trigger")

    comparison = compare_event_streams(annotations, physical)
    selection = select_event_source(comparison)

    assert comparison.discrepancy_class == "coherent_exact_agreement"
    assert comparison.ordered_candidate_streams_match
    assert selection.selected_source == "annotation"
    assert selection.selection_status == "annotation_selected"
    assert not selection.normalization_applied


def test_annotation_absence_with_coherent_physical_stream_selects_fallback() -> None:
    """A coherent raw physical stream may be selected when annotations are absent."""

    comparison = compare_event_streams(
        [],
        coherent_stream(source_type="physical_trigger"),
    )
    selection = select_event_source(comparison)

    assert comparison.discrepancy_class == "physical_fallback_available"
    assert selection.selected_source == "physical_trigger"
    assert selection.selection_status == "physical_fallback_selected"
    assert selection.unresolved_issue == ""


def test_annotation_physical_off_by_one_codes_do_not_replace_annotations() -> None:
    """Candidate code agreement is recorded without applying normalization."""

    annotations = coherent_stream(source_type="annotation")
    physical = coherent_stream(
        source_type="physical_trigger",
        raw_codes=(27, 30, 43, 46, 60),
        candidate_codes=(28, 30, 44, 46, 60),
    )

    comparison = compare_event_streams(annotations, physical)
    selection = select_event_source(comparison)
    rows = selected_event_rows(annotations, physical, selection)

    assert comparison.discrepancy_class == "coherent_source_code_difference"
    assert comparison.raw_code_difference_count == 4
    assert selection.selected_source == "annotation"
    assert all(
        row["raw_code"] == row["selected_normalized_code"] for row in rows
    )
    assert all(row["normalization_applied"] is False for row in rows)


def test_candidate_only_physical_stream_remains_unresolved() -> None:
    """Physical fallback cannot silently depend on off-by-one normalization."""

    physical = coherent_stream(
        source_type="physical_trigger",
        raw_codes=(27, 30, 43, 46, 60),
        candidate_codes=(28, 30, 44, 46, 60),
    )
    comparison = compare_event_streams([], physical)
    selection = select_event_source(comparison)

    assert (
        comparison.discrepancy_class
        == "physical_candidate_normalization_review"
    )
    assert selection.selected_source is None
    assert (
        selection.selection_status
        == "unresolved_candidate_normalization_required"
    )
    assert (
        selection.unresolved_issue
        == "physical_code_normalization_requires_review"
    )


def test_sparse_unrecoverable_sources_remain_unresolved() -> None:
    """Markers or isolated non-task codes must not become a task event source."""

    annotations = [event(10, 127, source_type="annotation")]
    physical = [event(10, 127, source_type="physical_trigger")]

    comparison = compare_event_streams(annotations, physical)
    selection = select_event_source(comparison)

    assert comparison.discrepancy_class == "no_coherent_task_stream"
    assert selection.selected_source is None
    assert selection.selection_status == "no_coherent_source"
    assert selection.unresolved_issue == "no_coherent_task_event_source"


def test_coherent_but_timing_discrepant_sources_keep_unresolved_flag() -> None:
    """Preferred annotations do not erase a physical timing discrepancy."""

    annotations = coherent_stream(source_type="annotation")
    physical = coherent_stream(
        source_type="physical_trigger", start_sample=107
    )

    comparison = compare_event_streams(annotations, physical)
    selection = select_event_source(comparison)

    assert comparison.discrepancy_class == "coherent_stream_discrepancy"
    assert selection.selected_source == "annotation"
    assert (
        selection.selection_status
        == "annotation_selected_with_source_discrepancy"
    )
    assert selection.unresolved_issue == "annotation_physical_stream_discrepancy"


def test_selected_rows_preserve_raw_and_candidate_provenance() -> None:
    """Selected output must retain raw codes separately from candidate codes."""

    annotations: list[dict[str, object]] = []
    physical = coherent_stream(source_type="physical_trigger")
    comparison = compare_event_streams(annotations, physical)
    selection = select_event_source(comparison)
    rows = selected_event_rows(annotations, physical, selection)

    assert rows
    assert all("raw_code" in row for row in rows)
    assert all("normalized_candidate_code" in row for row in rows)
    assert all(row["selected_source"] == "physical_trigger" for row in rows)
    assert all(row["selection_reason"] for row in rows)
    assert all(row["normalization_applied"] is False for row in rows)


def test_selected_annotation_marker_name_is_preserved() -> None:
    """A nonnumeric selected annotation marker must retain its source label."""

    annotations = coherent_stream(source_type="annotation")
    annotations.append(
        {
            "sample": 500,
            "raw_code": None,
            "normalized_candidate_code": None,
            "candidate_event_name": "file start",
            "source_type": "annotation",
        }
    )
    physical = coherent_stream(source_type="physical_trigger")
    selection = select_event_source(
        compare_event_streams(annotations, physical)
    )
    rows = selected_event_rows(annotations, physical, selection)

    marker = [row for row in rows if row.get("candidate_event_name") == "file start"]
    assert len(marker) == 1
    assert marker[0]["selected_event_name"] == "file start"
    assert marker[0]["selected_normalized_code"] is None


def test_unaffected_join_regression_excludes_only_changed_identity_surface() -> None:
    """Changed IDs may differ while the unaffected regression remains exact."""

    common = {
        "audit_trial_count": 1,
        "source_filename": "demi_01 Data.edf",
        "join_status": "raw_annotation_and_offset",
        "old_trace_epoch_duration_mismatch": False,
        "old_stimulus_duration_mismatch": False,
        "missing_expected_event_names": "",
        "extra_expected_event_names": "",
        "alignment_problem_codes": "",
        "clean_timing_row": True,
    }
    old = pd.DataFrame(
        [
            {"audit_participant_id": 1, **common},
            {
                "audit_participant_id": 11,
                **{**common, "source_filename": "demi_11 Data.edf"},
            },
        ]
    )
    corrected = pd.DataFrame(
        [
            {"audit_participant_id": 1, **common},
            {
                "audit_participant_id": 11,
                **{**common, "source_filename": "demi_13 Data.edf"},
            },
        ]
    )

    assert corrected_evidence.normalized_regression_surface(old).equals(
        corrected_evidence.normalized_regression_surface(corrected)
    )


def test_join_metrics_uses_canonical_clean_timing_column() -> None:
    """Corrected aggregate metrics must count clean_timing_row values."""

    join = pd.DataFrame(
        {
            "join_status": [
                "raw_annotation_and_offset",
                "raw_annotation_without_offset",
            ],
            "clean_timing_row": [True, False],
            "alignment_problem_codes": ["", "raw_annotation_without_offset"],
        }
    )
    metrics = corrected_evidence.join_metrics(join)

    assert metrics["strict_clean_timing"] == 1
    assert metrics["direct_alignment_problem_rows"] == 1


@pytest.mark.skipif(
    not OLD_COMPATIBLE_OFFSETS_PATH.exists(),
    reason="ignored local old-compatible offset target is unavailable",
)
def test_old_compatible_offset_target_is_unchanged() -> None:
    """The remediation must not modify the frozen behavioural offset target."""

    digest = hashlib.sha256(
        OLD_COMPATIBLE_OFFSETS_PATH.read_bytes()
    ).hexdigest()
    assert digest == OLD_COMPATIBLE_OFFSETS_SHA256

    with OLD_COMPATIBLE_OFFSETS_PATH.open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 9444
