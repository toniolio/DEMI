"""Focused tests for DEMI participant identity and event-source selection.

The tests use small synthetic event streams and the tracked minimal identity
contract. They do not open private EDF, task, TraceLab, masterlist, or RDS
files. One optional local regression test checks the frozen old-compatible
offset CSV when that ignored file is available.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

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
