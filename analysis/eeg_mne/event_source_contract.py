"""Pure identity and event-source contracts for the DEMI EEG reanalysis.

This module contains the small, testable rules that sit between raw EEG files
and later event-alignment audits. It exists because DEMI uses several ID
namespaces (task/database, behavioural, figure, and EEG recording), and because
an EDF can expose task events through annotations, the physical Trigger
channel, or both. Neither same-number identity nor annotation availability is
a safe implicit assumption.

The module deliberately contains no MNE file I/O. The dual-source inventory
script supplies one-row-per-event dictionaries, and these helpers:

- load and validate the tracked minimal EEG-to-behaviour identity contract;
- require explicit coverage for every discovered EEG recording ID;
- assess whether an event stream contains a coherent DEMI task sequence;
- compare annotation and physical-Trigger streams;
- select a conservative evidence source without silently applying candidate
  trigger-code normalization.

The tracked identity CSV contains only the minimal reproducible mapping needed
by public code. Detailed historical notes, task/database mappings, confidence,
and unresolved questions remain in the ignored private crosswalk.

This module does not preprocess EEG, mutate annotations or Trigger samples,
repair events, make participant/trial inclusion decisions, or construct
epochs. A selected source is an event-evidence choice only.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


IDENTITY_CONTRACT_VERSION = "v1"

# The five dominant event classes define the repeated DEMI trial sequence.
# Vividness (62) occurs only at block-level intervals and is therefore optional
# when assessing whether a recording contains a coherent task stream.
CORE_TASK_EVENT_CODES = (28, 30, 44, 46, 60)
OPTIONAL_TASK_EVENT_CODES = (62,)
TASK_EVENT_CODES = CORE_TASK_EVENT_CODES + OPTIONAL_TASK_EVENT_CODES

EVENT_CODE_TO_NAME = {
    28: "stim_on",
    30: "red_on",
    44: "trace_start",
    46: "trace_end",
    60: "accuracy_submit",
    62: "vividness_submit",
}

# These differences were observed between physical Trigger values and the EDF
# annotation labels at the same samples. They are candidates only. The
# inventory records them so the two sources can be compared, but the selection
# contract never silently applies them to a physical fallback stream.
PHYSICAL_CODE_NORMALIZATION_CANDIDATES = {
    27: (28, "candidate_27_to_28_from_same_sample_annotation_comparison"),
    43: (44, "candidate_43_to_44_from_same_sample_annotation_comparison"),
    61: (62, "candidate_61_to_62_from_same_sample_annotation_comparison"),
}

VALID_MAPPING_STATUSES = {
    "reviewed_same_identifier",
    "explicit_remap",
    "explicitly_unmapped",
}


@dataclass(frozen=True)
class IdentityMapping:
    """One explicit EEG-recording to behavioural-participant mapping.

    Attributes:
        eeg_recording_id: Integer ID encoded in the raw EEG filename.
        behavioural_participant_id: Behavioural ID used by the old-compatible
            offset table, or None when the EEG recording is deliberately
            not mapped.
        mapping_status: Reviewed contract status.
        contract_version: Version string shared by all rows in the CSV.
    """

    eeg_recording_id: int
    behavioural_participant_id: int | None
    mapping_status: str
    contract_version: str


@dataclass(frozen=True)
class StreamAssessment:
    """Deterministic coherence summary for one event stream.

    Attributes:
        status: Coherent, empty, missing-code, or imbalanced status.
        coherent: Whether all dominant task codes are present with reasonably
            balanced counts.
        task_event_count: Count of dominant plus optional task events.
        core_code_counts: Count for each dominant task code.
        optional_code_counts: Count for optional block-level codes.
        minimum_core_count: Smallest dominant-code count.
        maximum_core_count: Largest dominant-code count.
    """

    status: str
    coherent: bool
    task_event_count: int
    core_code_counts: dict[int, int]
    optional_code_counts: dict[int, int]
    minimum_core_count: int
    maximum_core_count: int


@dataclass(frozen=True)
class StreamComparison:
    """Annotation-versus-physical comparison for one EDF file.

    Attributes:
        discrepancy_class: Descriptive comparison class.
        annotation: Assessment of annotation raw codes.
        physical_raw: Assessment of physical masked codes without candidate
            normalization.
        physical_candidate: Assessment using candidate codes for comparison
            only.
        compared_event_count: Ordered event pairs compared.
        same_sample_and_candidate_code_count: Pairs agreeing in sample and
            candidate semantic code.
        raw_code_difference_count: Compared pairs whose raw codes differ.
        maximum_absolute_sample_difference: Largest ordered-pair sample
            difference, or None when no pairs exist.
        ordered_candidate_streams_match: Whether task-event counts, samples,
            and candidate codes agree within the sample tolerance.
    """

    discrepancy_class: str
    annotation: StreamAssessment
    physical_raw: StreamAssessment
    physical_candidate: StreamAssessment
    compared_event_count: int
    same_sample_and_candidate_code_count: int
    raw_code_difference_count: int
    maximum_absolute_sample_difference: int | None
    ordered_candidate_streams_match: bool


@dataclass(frozen=True)
class SourceSelection:
    """Conservative event-source selection for one EDF file.

    Attributes:
        selected_source: Annotation, physical_trigger, or None.
        selection_status: Stable descriptive status.
        selection_reason: Human-readable deterministic reason.
        unresolved_issue: Empty when no source-level issue remains, otherwise
            a stable unresolved classification.
        normalization_applied: Always False in contract version 1.
    """

    selected_source: str | None
    selection_status: str
    selection_reason: str
    unresolved_issue: str
    normalization_applied: bool


def _text(value: Any) -> str:
    """Return stripped text, treating None as empty.

    Args:
        value: Arbitrary scalar value.

    Returns:
        Stripped text or an empty string.

    Side effects:
        None.
    """

    if value is None:
        return ""
    return str(value).strip()


def integer_or_none(value: Any) -> int | None:
    """Convert a scalar to an integer when possible.

    Args:
        value: Integer-like value, blank text, or None.

    Returns:
        Integer value or None for blank/missing input.

    Side effects:
        None. Raises ValueError for nonblank non-integer input.
    """

    text = _text(value)
    if not text:
        return None
    return int(float(text))


def load_identity_contract(path: Path) -> dict[int, IdentityMapping]:
    """Read and fail-closed validate the public identity contract CSV.

    Args:
        path: CSV containing one row per EEG recording ID.

    Returns:
        Dictionary keyed by EEG recording ID.

    Side effects:
        Reads path. Raises RuntimeError for missing files, schema errors,
        duplicate IDs, invalid statuses, inconsistent versions, or invalid
        status/value combinations.
    """

    if not path.exists():
        raise RuntimeError(f"identity contract is missing: {path}")

    required = {
        "eeg_recording_id",
        "behavioural_participant_id",
        "mapping_status",
        "contract_version",
    }
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                "identity contract is missing required column(s): "
                + ", ".join(sorted(missing))
            )

        mappings: dict[int, IdentityMapping] = {}
        behavioural_to_eeg: dict[int, int] = {}
        for row_number, row in enumerate(reader, start=2):
            try:
                eeg_id = integer_or_none(row["eeg_recording_id"])
                behavioural_id = integer_or_none(row["behavioural_participant_id"])
            except ValueError as error:
                raise RuntimeError(
                    f"identity contract row {row_number} has a non-integer ID"
                ) from error

            if eeg_id is None:
                raise RuntimeError(
                    f"identity contract row {row_number} has no EEG recording ID"
                )
            if eeg_id in mappings:
                raise RuntimeError(
                    f"identity contract repeats EEG recording ID {eeg_id}"
                )

            status = _text(row["mapping_status"])
            version = _text(row["contract_version"])
            if status not in VALID_MAPPING_STATUSES:
                raise RuntimeError(
                    f"identity contract row {row_number} has invalid mapping_status {status!r}"
                )
            if version != IDENTITY_CONTRACT_VERSION:
                raise RuntimeError(
                    f"identity contract row {row_number} has unsupported version {version!r}"
                )

            if status == "reviewed_same_identifier" and behavioural_id != eeg_id:
                raise RuntimeError(
                    f"same-identifier row for EEG {eeg_id} must map to the same behavioural ID"
                )
            if status == "explicit_remap" and (
                behavioural_id is None or behavioural_id == eeg_id
            ):
                raise RuntimeError(
                    f"explicit-remap row for EEG {eeg_id} must map to a different behavioural ID"
                )
            if status == "explicitly_unmapped" and behavioural_id is not None:
                raise RuntimeError(
                    f"explicitly-unmapped row for EEG {eeg_id} must have a blank behavioural ID"
                )

            if behavioural_id is not None:
                prior_eeg = behavioural_to_eeg.get(behavioural_id)
                if prior_eeg is not None:
                    raise RuntimeError(
                        "identity contract maps behavioural participant "
                        f"{behavioural_id} from both EEG {prior_eeg} and EEG {eeg_id}"
                    )
                behavioural_to_eeg[behavioural_id] = eeg_id

            mappings[eeg_id] = IdentityMapping(
                eeg_recording_id=eeg_id,
                behavioural_participant_id=behavioural_id,
                mapping_status=status,
                contract_version=version,
            )

    if not mappings:
        raise RuntimeError("identity contract contains no mappings")
    return mappings


def require_identity_mapping(
    mappings: Mapping[int, IdentityMapping],
    eeg_recording_id: int,
) -> IdentityMapping:
    """Return the explicit mapping for an EEG recording ID.

    Args:
        mappings: Mapping dictionary from load_identity_contract.
        eeg_recording_id: EEG recording ID to look up.

    Returns:
        The explicit IdentityMapping.

    Side effects:
        None. Raises RuntimeError when the ID is absent. There is intentionally
        no same-number fallback.
    """

    try:
        return mappings[int(eeg_recording_id)]
    except KeyError as error:
        raise RuntimeError(
            "EEG recording ID "
            f"{int(eeg_recording_id)} is absent from the explicit identity contract"
        ) from error


def validate_identity_coverage(
    mappings: Mapping[int, IdentityMapping],
    discovered_eeg_ids: Sequence[int],
) -> None:
    """Require explicit mapping rows for every discovered EEG recording ID.

    Args:
        mappings: Mapping dictionary from load_identity_contract.
        discovered_eeg_ids: IDs parsed from raw EDF filenames.

    Returns:
        None.

    Side effects:
        None. Raises RuntimeError when any discovered ID is absent.
    """

    missing = sorted(
        set(int(value) for value in discovered_eeg_ids).difference(mappings)
    )
    if missing:
        raise RuntimeError(
            "identity contract lacks discovered EEG recording ID(s): "
            + ", ".join(str(value) for value in missing)
        )


def physical_code_candidate(raw_code: int | None) -> tuple[int | None, str]:
    """Return a comparison-only normalized candidate for a physical code.

    Args:
        raw_code: Masked physical Trigger code.

    Returns:
        Candidate code and rule. Known one-off source differences receive a
        candidate and explicit rule; other codes are preserved. Missing input
        returns None and an empty rule.

    Side effects:
        None. Returning a candidate does not authorize its application.
    """

    if raw_code is None:
        return None, ""
    if raw_code in PHYSICAL_CODE_NORMALIZATION_CANDIDATES:
        return PHYSICAL_CODE_NORMALIZATION_CANDIDATES[raw_code]
    return int(raw_code), "raw_code_preserved"


def event_name_for_code(code: int | None) -> str:
    """Map a normalized candidate code to a DEMI semantic event name.

    Args:
        code: Numeric candidate code or None.

    Returns:
        Semantic event name, or an empty string when unmapped.

    Side effects:
        None.
    """

    if code is None:
        return ""
    return EVENT_CODE_TO_NAME.get(int(code), "")


def _codes_from_events(
    events: Sequence[Mapping[str, Any]],
    code_field: str,
) -> list[int]:
    """Extract integer task-event codes from event rows.

    Args:
        events: Event dictionaries.
        code_field: Field containing raw or comparison-candidate codes.

    Returns:
        Codes restricted to TASK_EVENT_CODES.

    Side effects:
        None.
    """

    codes: list[int] = []
    for event in events:
        code = integer_or_none(event.get(code_field))
        if code in TASK_EVENT_CODES:
            codes.append(int(code))
    return codes


def assess_stream(
    events: Sequence[Mapping[str, Any]],
    *,
    code_field: str,
    maximum_core_count_ratio: float = 1.5,
) -> StreamAssessment:
    """Assess whether rows contain a coherent repeated DEMI task stream.

    Args:
        events: Event dictionaries for one file and one source.
        code_field: Field containing the codes to assess.
        maximum_core_count_ratio: Largest permitted maximum/minimum count ratio
            across the five dominant event classes.

    Returns:
        StreamAssessment with descriptive status and counts.

    Side effects:
        None.
    """

    task_codes = _codes_from_events(events, code_field)
    counts = Counter(task_codes)
    core_counts = {
        code: int(counts.get(code, 0)) for code in CORE_TASK_EVENT_CODES
    }
    optional_counts = {
        code: int(counts.get(code, 0)) for code in OPTIONAL_TASK_EVENT_CODES
    }
    minimum = min(core_counts.values()) if core_counts else 0
    maximum = max(core_counts.values()) if core_counts else 0

    if not task_codes:
        status = "no_task_events"
        coherent = False
    elif minimum == 0:
        status = "missing_core_codes"
        coherent = False
    elif maximum / minimum > maximum_core_count_ratio:
        status = "imbalanced_core_counts"
        coherent = False
    else:
        status = "coherent"
        coherent = True

    return StreamAssessment(
        status=status,
        coherent=coherent,
        task_event_count=len(task_codes),
        core_code_counts=core_counts,
        optional_code_counts=optional_counts,
        minimum_core_count=minimum,
        maximum_core_count=maximum,
    )


def _ordered_task_events(
    events: Sequence[Mapping[str, Any]],
    code_field: str,
) -> list[tuple[int, int, int]]:
    """Return ordered sample, raw-code, and candidate-code task tuples.

    Args:
        events: Event dictionaries.
        code_field: Field used to decide whether a row is a task event.

    Returns:
        Stable sample-ordered tuples. Raw code is -1 when unavailable.

    Side effects:
        None.
    """

    ordered: list[tuple[int, int, int]] = []
    for event in events:
        candidate = integer_or_none(event.get(code_field))
        if candidate not in TASK_EVENT_CODES:
            continue
        sample = integer_or_none(event.get("sample"))
        if sample is None:
            raise ValueError("task event is missing sample provenance")
        raw_code = integer_or_none(event.get("raw_code"))
        ordered.append(
            (sample, raw_code if raw_code is not None else -1, candidate)
        )
    return sorted(ordered, key=lambda value: value[0])


def compare_event_streams(
    annotation_events: Sequence[Mapping[str, Any]],
    physical_events: Sequence[Mapping[str, Any]],
    *,
    sample_tolerance: int = 1,
) -> StreamComparison:
    """Compare annotation and physical Trigger streams descriptively.

    Args:
        annotation_events: Annotation-source event rows.
        physical_events: Physical-Trigger event rows.
        sample_tolerance: Maximum absolute sample difference considered the
            same acquisition event.

    Returns:
        StreamComparison retaining raw-code differences while using physical
        candidate codes only for comparison.

    Side effects:
        None.
    """

    annotation = assess_stream(annotation_events, code_field="raw_code")
    physical_raw = assess_stream(physical_events, code_field="raw_code")
    physical_candidate = assess_stream(
        physical_events, code_field="normalized_candidate_code"
    )

    annotation_ordered = _ordered_task_events(annotation_events, "raw_code")
    physical_ordered = _ordered_task_events(
        physical_events, "normalized_candidate_code"
    )
    compared = min(len(annotation_ordered), len(physical_ordered))
    same = 0
    raw_differences = 0
    sample_differences: list[int] = []
    for annotation_event, physical_event in zip(
        annotation_ordered, physical_ordered
    ):
        sample_difference = abs(annotation_event[0] - physical_event[0])
        sample_differences.append(sample_difference)
        if annotation_event[1] != physical_event[1]:
            raw_differences += 1
        if (
            sample_difference <= sample_tolerance
            and annotation_event[2] == physical_event[2]
        ):
            same += 1

    streams_match = (
        len(annotation_ordered) == len(physical_ordered)
        and compared > 0
        and same == compared
    )

    if annotation.coherent and physical_candidate.coherent:
        if streams_match:
            discrepancy = (
                "coherent_exact_agreement"
                if raw_differences == 0
                else "coherent_source_code_difference"
            )
        else:
            discrepancy = "coherent_stream_discrepancy"
    elif annotation.coherent:
        discrepancy = "annotation_coherent_physical_unusable"
    elif physical_raw.coherent:
        discrepancy = "physical_fallback_available"
    elif physical_candidate.coherent:
        discrepancy = "physical_candidate_normalization_review"
    else:
        discrepancy = "no_coherent_task_stream"

    return StreamComparison(
        discrepancy_class=discrepancy,
        annotation=annotation,
        physical_raw=physical_raw,
        physical_candidate=physical_candidate,
        compared_event_count=compared,
        same_sample_and_candidate_code_count=same,
        raw_code_difference_count=raw_differences,
        maximum_absolute_sample_difference=(
            max(sample_differences) if sample_differences else None
        ),
        ordered_candidate_streams_match=streams_match,
    )


def select_event_source(comparison: StreamComparison) -> SourceSelection:
    """Select an event-evidence source under the conservative version-1 rule.

    Args:
        comparison: Annotation-versus-physical comparison.

    Returns:
        SourceSelection. Candidate physical code normalization is never applied
        silently.

    Side effects:
        None.
    """

    if comparison.annotation.coherent:
        if comparison.discrepancy_class == "coherent_stream_discrepancy":
            return SourceSelection(
                selected_source="annotation",
                selection_status="annotation_selected_with_source_discrepancy",
                selection_reason=(
                    "annotations are coherent and remain the preferred source; "
                    "the physical discrepancy is retained for review"
                ),
                unresolved_issue="annotation_physical_stream_discrepancy",
                normalization_applied=False,
            )
        return SourceSelection(
            selected_source="annotation",
            selection_status="annotation_selected",
            selection_reason="coherent annotations preferred by contract",
            unresolved_issue="",
            normalization_applied=False,
        )

    if comparison.physical_raw.coherent:
        annotation_reason = (
            "absent"
            if comparison.annotation.task_event_count == 0
            else "insufficient"
        )
        return SourceSelection(
            selected_source="physical_trigger",
            selection_status="physical_fallback_selected",
            selection_reason=(
                f"annotation task stream is {annotation_reason}; "
                "physical raw codes form a coherent task stream"
            ),
            unresolved_issue="",
            normalization_applied=False,
        )

    if comparison.physical_candidate.coherent:
        return SourceSelection(
            selected_source=None,
            selection_status="unresolved_candidate_normalization_required",
            selection_reason=(
                "physical stream is coherent only after comparison-only "
                "candidate normalization"
            ),
            unresolved_issue="physical_code_normalization_requires_review",
            normalization_applied=False,
        )

    return SourceSelection(
        selected_source=None,
        selection_status="no_coherent_source",
        selection_reason="neither source contains a coherent DEMI task stream",
        unresolved_issue="no_coherent_task_event_source",
        normalization_applied=False,
    )


def selected_event_rows(
    annotation_events: Sequence[Mapping[str, Any]],
    physical_events: Sequence[Mapping[str, Any]],
    selection: SourceSelection,
) -> list[dict[str, Any]]:
    """Return copied rows from the selected source with selection provenance.

    Args:
        annotation_events: Annotation-source rows.
        physical_events: Physical-source rows.
        selection: Result from select_event_source.

    Returns:
        Copied selected rows. Raw and candidate codes remain separate, and
        normalization_applied is recorded on every row.

    Side effects:
        None.
    """

    if selection.selected_source == "annotation":
        source_rows = annotation_events
    elif selection.selected_source == "physical_trigger":
        source_rows = physical_events
    else:
        return []

    selected: list[dict[str, Any]] = []
    for row in source_rows:
        copied = dict(row)
        copied["selected_source"] = selection.selected_source
        copied["selection_status"] = selection.selection_status
        copied["selection_reason"] = selection.selection_reason
        copied["selection_unresolved_issue"] = selection.unresolved_issue
        copied["normalization_applied"] = selection.normalization_applied
        # Contract v1 only selects physical streams whose raw codes are already
        # coherent. For annotations, the annotation label is already the
        # authoritative selected representation. No candidate rewrite occurs.
        copied["selected_normalized_code"] = integer_or_none(
            copied.get("raw_code")
        )
        copied["selected_event_name"] = (
            event_name_for_code(copied["selected_normalized_code"])
            or _text(copied.get("candidate_event_name"))
        )
        selected.append(copied)
    return selected
