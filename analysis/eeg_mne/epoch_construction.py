"""Contracts and tested helpers for DEMI EEG epoch construction.

This module implements the accepted, signal-neutral epoch-definition contract
for the active MNE reanalysis.  It turns already-authoritative ledger timing
into deterministic MNE event samples and metadata.  The numbered stage-15
driver owns local file discovery, bounded FIF reads, resumable writes, and run
aggregation.

Inputs:
    The canonical ``event_epoch_eligibility_v1`` ledger, one linked
    ``continuous_v2`` post-ICA Raw object at a time, and the serialized raw
    trial-event provenance already carried by each ledger row.

Outputs:
    Validated future-ready row selections, derived response-onset,
    response-end, and ``red_on`` anchors, MNE-compatible event arrays and
    metadata, saved-artifact validation evidence, and deterministic hashes.

This module explicitly does not inspect FIF annotations, repair events,
resample, baseline voltages, reject amplitudes or annotations, run AutoReject,
apply CSD, compute time-frequency power, exclude participants, or create an
ID-86 sensitivity route.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import uuid

import mne
import numpy as np
import pandas as pd

from event_epoch_eligibility import bool_value, stable_frame_hash, text_value


EPOCH_NAMESPACE = "epochs_v1"
EPOCH_SCHEMA_VERSION = 1
LEDGER_VERSION = "event_epoch_eligibility_v1"
EXPECTED_EPOCH_COUNT = 8_798
EXPECTED_STRICT_CLEAN_COUNT = 8_789
EXPECTED_DURATION_WARNING_COUNT = 9
EXPECTED_FILE49_COUNT = 117
EXPECTED_SFREQ = 1_000.0
HALF_SAMPLE_AMBIGUITY_TOLERANCE = 1e-9


@dataclass(frozen=True)
class EpochFamily:
    """One accepted independently traceable epoch family.

    Args:
        name: Stable machine-readable family name.
        event_id: Positive MNE event code used only inside this family.
        tmin: Inclusive epoch start in seconds relative to the anchor.
        tmax: Inclusive epoch end in seconds relative to the anchor.
        scientific_role: Concise interpretation of the saved family.
    """

    name: str
    event_id: int
    tmin: float
    tmax: float
    scientific_role: str


FAMILIES: tuple[EpochFamily, ...] = (
    EpochFamily(
        name="response_onset",
        event_id=1,
        tmin=-1.5,
        tmax=2.5,
        scientific_role="Response-onset-locked broad task support.",
    ),
    EpochFamily(
        name="response_end",
        event_id=2,
        tmin=-1.5,
        tmax=2.5,
        scientific_role=(
            "Response-end-locked broad task support preserving pre-movement, "
            "movement, and post-movement periods."
        ),
    ),
    EpochFamily(
        name="red_on",
        event_id=3,
        tmin=-1.5,
        tmax=0.8,
        scientific_role=(
            "red_on-locked support for a later spectral reference; no spectral "
            "normalization is applied here."
        ),
    ),
)

FAMILY_BY_NAME = {family.name: family for family in FAMILIES}

REQUIRED_LEDGER_COLUMNS: tuple[str, ...] = (
    "ledger_version",
    "event_policy_version",
    "canonical_event_key",
    "behavioural_participant_id",
    "eeg_recording_id",
    "source_filename",
    "file_role",
    "source_type",
    "selected_source",
    "task_filename",
    "audit_trial_count",
    "raw_trial_sequence",
    "join_trial_count",
    "offset_trial_count",
    "physical",
    "trace_start_onset_seconds",
    "trace_end_onset_seconds",
    "red_on_onset_seconds",
    "real_start",
    "real_end",
    "trial_event_provenance_json",
    "primary_ledger_eligibility",
    "primary_ledger_eligibility_reason_code",
    "strict_clean_only_sensitivity_eligibility",
    "strict_clean_only_sensitivity_eligibility_reason_code",
    "duration_warning_only",
    "event_policy_reason_code",
    "accepted_policy_ids",
    "accepted_specific_policy_id",
    "accepted_policy_action",
    "continuous_v2_run_id",
    "continuous_v2_manifest_path",
    "continuous_v2_manifest_id",
    "continuous_derivative_kind",
    "continuous_derivative_path",
    "continuous_source_sha256",
    "continuous_qc_warning_codes",
    "interpolated_channel_count",
    "interpolated_channel_proportion",
    "interpolation_proportion_denominator",
    "ica_terminal_route",
    "post_ica_derivative_availability",
    "post_ica_derivative_availability_reason_code",
    "future_epoch_construction_readiness",
    "future_epoch_construction_readiness_reason_code",
)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest for one file.

    Args:
        path: File to hash without modification.
        chunk_size: Number of bytes read per chunk.

    Returns:
        Lowercase hexadecimal SHA-256 digest.

    Side effects:
        Reads ``path`` sequentially.
    """

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize a JSON-compatible value deterministically and reject NaN."""

    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def content_fingerprint(value: Any) -> str:
    """Return a deterministic SHA-256 digest for JSON-compatible content."""

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def atomic_write_json(path: Path, value: Any) -> None:
    """Write deterministic JSON through a sibling temporary file.

    Args:
        path: Final JSON path.
        value: JSON-compatible content.

    Side effects:
        Creates a temporary file, fsyncs it, and atomically replaces ``path``.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}.json")
    payload = canonical_json_bytes(value) + b"\n"
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_write_text(path: Path, value: str) -> None:
    """Write UTF-8 text through a sibling temporary file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}{path.suffix}")
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_write_parquet(path: Path, frame: pd.DataFrame) -> None:
    """Write a parquet table through a sibling temporary file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}.parquet")
    frame.to_parquet(temporary, index=False)
    with temporary.open("rb") as handle:
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Write a CSV review table through a sibling temporary file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}.csv")
    frame.to_csv(temporary, index=False, lineterminator="\n")
    with temporary.open("rb") as handle:
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def expected_epoch_sample_count(tmin: float, tmax: float, sfreq: float) -> int:
    """Return the inclusive MNE endpoint sample count.

    Both requested endpoints must fall exactly on the sampling grid.  MNE
    includes both endpoints, so the count is ``(tmax - tmin) * sfreq + 1``.
    At 1000 Hz this yields 4,001 samples for -1.5 to +2.5 seconds and 2,301
    samples for -1.5 to +0.8 seconds.
    """

    if not math.isfinite(sfreq) or sfreq <= 0:
        raise ValueError(f"invalid sampling frequency: {sfreq}")
    start = tmin * sfreq
    stop = tmax * sfreq
    for label, value in (("tmin", start), ("tmax", stop)):
        if not np.isclose(value, np.rint(value), atol=1e-9, rtol=0.0):
            raise ValueError(f"{label}={value / sfreq} is not on the {sfreq}-Hz sample grid")
    return int(np.rint(stop) - np.rint(start) + 1)


def seconds_to_sample(
    seconds: float,
    sfreq: float,
    *,
    first_samp: int = 0,
    ambiguity_tolerance: float = HALF_SAMPLE_AMBIGUITY_TOLERANCE,
) -> int:
    """Convert recording-relative seconds to an absolute MNE event sample.

    MNE's ``time_as_index(..., use_rounding=True)`` uses nearest-sample
    rounding.  This helper reproduces that NumPy rounding rule explicitly and
    fails closed if a timestamp is an exact half-sample tie within the stated
    numerical tolerance.  The returned event sample includes ``first_samp``.

    Args:
        seconds: Anchor time in seconds from recording start.
        sfreq: Sampling frequency in Hz.
        first_samp: Raw object's absolute first sample.
        ambiguity_tolerance: Maximum distance from a half-sample tie that is
            treated as ambiguous.

    Returns:
        Absolute integer MNE event sample.

    Raises:
        ValueError: For nonfinite/negative times, invalid sampling, or an
            ambiguous half-sample tie.
    """

    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError(f"anchor seconds must be finite and nonnegative: {seconds}")
    if not math.isfinite(sfreq) or sfreq <= 0:
        raise ValueError(f"invalid sampling frequency: {sfreq}")
    scaled = seconds * sfreq
    fraction = scaled - math.floor(scaled)
    if abs(fraction - 0.5) <= ambiguity_tolerance:
        raise ValueError(
            f"anchor {seconds:.15g} s is an ambiguous half-sample at {sfreq:.15g} Hz"
        )
    return int(np.rint(scaled)) + int(first_samp)


def _required_numeric(row: Mapping[str, Any], column: str) -> float:
    """Return one required finite numeric ledger scalar."""

    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"ledger row lacks numeric {column}") from error
    if not math.isfinite(value):
        raise ValueError(f"ledger row has nonfinite {column}")
    return value


def parse_trial_event_provenance(value: Any) -> list[dict[str, Any]]:
    """Parse and minimally validate one accepted trial-event provenance list."""

    try:
        parsed = json.loads(text_value(value))
    except json.JSONDecodeError as error:
        raise ValueError("trial_event_provenance_json is invalid JSON") from error
    if not isinstance(parsed, list) or not parsed:
        raise ValueError("trial_event_provenance_json must be a nonempty list")
    if not all(isinstance(item, dict) for item in parsed):
        raise ValueError("trial_event_provenance_json contains a non-object entry")
    return parsed


def raw_event_for_role(row: Mapping[str, Any], role: str) -> dict[str, Any]:
    """Return the unique accepted raw-provenance entry for an event role.

    The lookup uses ``proposed_event_name`` because the accepted event policy
    authorizes the provenance-preserving proposed cleanup layer.  It never
    consults a FIF annotation.
    """

    matches = [
        item
        for item in parse_trial_event_provenance(row["trial_event_provenance_json"])
        if text_value(item.get("proposed_event_name")) == role
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one {role} provenance event, observed {len(matches)}")
    event = matches[0]
    _required_numeric(event, "onset_seconds")
    _required_numeric(event, "sample")
    return event


def derive_anchor(row: Mapping[str, Any], family_name: str) -> dict[str, Any]:
    """Derive one accepted scientific anchor from canonical ledger fields.

    Args:
        row: One future-ready canonical ledger row.
        family_name: ``response_onset``, ``response_end``, or ``red_on``.

    Returns:
        Anchor time, formula, condition semantics, base raw-event role and
        provenance, and applied timing offset.  Sample conversion is separate
        because it must use the reopened continuous Raw sampling contract.
    """

    if family_name not in FAMILY_BY_NAME:
        raise KeyError(f"unknown epoch family: {family_name}")
    physical = bool_value(row.get("physical"))
    condition_semantics = "overt_movement" if physical else "motor_imagery"
    trace_start = _required_numeric(row, "trace_start_onset_seconds")

    if family_name == "response_onset":
        raw_role = "trace_start"
        if physical:
            offset = _required_numeric(row, "real_start")
            anchor_time = trace_start + offset
            source = "trace_start_plus_real_start"
            formula = "trace_start_onset_seconds + real_start"
        else:
            offset = 0.0
            anchor_time = trace_start
            source = "raw_trace_start"
            formula = "trace_start_onset_seconds"
    elif family_name == "response_end":
        if physical:
            raw_role = "trace_start"
            offset = _required_numeric(row, "real_end")
            anchor_time = trace_start + offset
            source = "trace_start_plus_real_end"
            formula = "trace_start_onset_seconds + real_end"
        else:
            raw_role = "trace_end"
            offset = 0.0
            anchor_time = _required_numeric(row, "trace_end_onset_seconds")
            source = "raw_trace_end"
            formula = "trace_end_onset_seconds"
    else:
        raw_role = "red_on"
        offset = 0.0
        anchor_time = _required_numeric(row, "red_on_onset_seconds")
        source = "raw_red_on"
        formula = "red_on_onset_seconds"

    raw_event = raw_event_for_role(row, raw_role)
    raw_onset = _required_numeric(raw_event, "onset_seconds")
    raw_sample = int(_required_numeric(raw_event, "sample"))
    ledger_raw_onset = _required_numeric(row, f"{raw_role}_onset_seconds")
    if not np.isclose(raw_onset, ledger_raw_onset, atol=1e-12, rtol=0.0):
        raise ValueError(
            f"{raw_role} ledger/provenance onset mismatch: {ledger_raw_onset} != {raw_onset}"
        )

    return {
        "derived_anchor_type": family_name,
        "condition_semantics": condition_semantics,
        "anchor_source": source,
        "anchor_formula": formula,
        "anchor_time_seconds": float(anchor_time),
        "applied_timing_offset_seconds": float(offset),
        "raw_anchor_event_role": raw_role,
        "raw_anchor_onset_seconds": raw_onset,
        "raw_anchor_sample": raw_sample,
        "raw_anchor_source_type": text_value(raw_event.get("source_type")),
        "raw_anchor_event_name": text_value(raw_event.get("raw_event_name")),
        "raw_anchor_value": text_value(raw_event.get("raw_value")),
        "raw_anchor_code": _required_numeric(raw_event, "raw_code"),
        "raw_anchor_annotation_description": text_value(
            raw_event.get("raw_annotation_description")
        ),
        "raw_anchor_proposed_event_name": text_value(
            raw_event.get("proposed_event_name")
        ),
        "raw_anchor_proposed_sequence_action": text_value(
            raw_event.get("proposed_sequence_action")
        ),
        "raw_anchor_proposed_sequence_reason": text_value(
            raw_event.get("proposed_sequence_reason")
        ),
    }


def select_epoch_rows(ledger: pd.DataFrame) -> pd.DataFrame:
    """Select and validate the exact accepted ordinary epoch surface.

    The returned order is the input ledger order.  A zero-based
    ``canonical_order_index`` is added so every sharded artifact can be placed
    back into the global family order without inference.
    """

    missing = [column for column in REQUIRED_LEDGER_COLUMNS if column not in ledger]
    if missing:
        raise RuntimeError(f"canonical ledger lacks required epoch columns: {missing}")
    if set(ledger["ledger_version"].dropna().astype(str)) != {LEDGER_VERSION}:
        raise RuntimeError("ledger version differs from event_epoch_eligibility_v1")
    if ledger["canonical_event_key"].duplicated().any():
        raise RuntimeError("canonical ledger keys are not unique")

    selected = ledger.loc[
        ledger["future_epoch_construction_readiness"].map(bool_value)
    ].copy()
    selected.insert(0, "canonical_order_index", np.arange(len(selected), dtype=np.int64))
    strict = selected["strict_clean_only_sensitivity_eligibility"].map(bool_value)
    warnings = selected["duration_warning_only"].map(bool_value)

    observed = {
        "selected": len(selected),
        "strict_clean": int(strict.sum()),
        "duration_warning": int(warnings.sum()),
    }
    expected = {
        "selected": EXPECTED_EPOCH_COUNT,
        "strict_clean": EXPECTED_STRICT_CLEAN_COUNT,
        "duration_warning": EXPECTED_DURATION_WARNING_COUNT,
    }
    if observed != expected:
        raise RuntimeError(f"accepted epoch surface differs: observed={observed}, expected={expected}")
    if not (warnings == ~strict).all():
        raise RuntimeError("strict-clean difference is not exactly the duration-warning rows")
    if not selected["primary_ledger_eligibility"].map(bool_value).all():
        raise RuntimeError("a future-ready row lacks primary eligibility")
    if not selected["post_ica_derivative_availability"].map(bool_value).all():
        raise RuntimeError("a future-ready row lacks post-ICA availability")
    if set(selected["continuous_derivative_kind"].astype(str)) != {"post_ica"}:
        raise RuntimeError("ordinary epoch rows do not all use post-ICA derivatives")
    if selected["eeg_recording_id"].astype(int).eq(86).any():
        raise RuntimeError("ID 86 entered the ordinary epoch surface")
    if selected["continuous_derivative_path"].fillna("").astype(str).str.strip().eq("").any():
        raise RuntimeError("a future-ready row lacks a derivative path")

    file49 = selected[selected["eeg_recording_id"].astype(int).eq(49)]
    if len(file49) != EXPECTED_FILE49_COUNT:
        raise RuntimeError(f"file 49 contributes {len(file49)}, expected {EXPECTED_FILE49_COUNT}")
    if not file49["continuous_qc_warning_codes"].str.contains(
        "global_bad_proportion_above_25_percent", na=False
    ).all():
        raise RuntimeError("file 49 interpolation warning is not visible on every selected row")
    if selected["source_filename"].astype(str).str.contains("demi_54_1", regex=False).any():
        raise RuntimeError("file 54_1 entered the accepted epoch surface")

    required_numeric = (
        "behavioural_participant_id",
        "eeg_recording_id",
        "audit_trial_count",
        "raw_trial_sequence",
        "join_trial_count",
        "offset_trial_count",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
        "red_on_onset_seconds",
        "real_start",
        "real_end",
    )
    for column in required_numeric:
        values = pd.to_numeric(selected[column], errors="coerce")
        if values.isna().any() or not np.isfinite(values.to_numpy(dtype=float)).all():
            raise RuntimeError(f"selected ledger rows have missing/nonfinite {column}")
    if not (
        selected["audit_trial_count"].astype(int)
        == selected["offset_trial_count"].astype(int)
    ).all():
        raise RuntimeError("accepted audit and offset trial identifiers differ")
    return selected.reset_index(drop=True)


def derive_offset_identifiers(offset_trial_count: Any) -> tuple[int, int, int]:
    """Recover accepted session, block, and within-block trial identifiers.

    The old-compatible offset contract numbers trials sequentially within the
    single task session as ``trial + (block - 1) * 20``.  Active selected rows
    all belong to session 1.  This function reverses that accepted mapping.
    """

    value = int(float(offset_trial_count))
    if value < 1:
        raise ValueError(f"offset_trial_count must be positive: {value}")
    return 1, (value - 1) // 20 + 1, (value - 1) % 20 + 1


def build_family_metadata(
    rows: pd.DataFrame,
    family: EpochFamily,
    *,
    sfreq: float,
    first_samp: int,
    last_samp: int,
) -> pd.DataFrame:
    """Build fully traceable MNE metadata for one recording/family shard.

    Every row is retained in input order.  Raw event sample/onset consistency,
    derived-anchor bounds, sample ties, and duplicate MNE event samples fail
    closed before an Epochs object is created.
    """

    records: list[dict[str, Any]] = []
    start_offset = int(np.rint(family.tmin * sfreq))
    stop_offset = int(np.rint(family.tmax * sfreq))

    for row in rows.to_dict(orient="records"):
        anchor = derive_anchor(row, family.name)
        raw_sample = seconds_to_sample(
            anchor["raw_anchor_onset_seconds"], sfreq, first_samp=first_samp
        )
        if raw_sample != anchor["raw_anchor_sample"]:
            raise RuntimeError(
                f"raw provenance sample mismatch for {row['canonical_event_key']} "
                f"{family.name}: converted={raw_sample}, ledger={anchor['raw_anchor_sample']}"
            )
        if anchor["raw_anchor_source_type"] != text_value(row["source_type"]):
            raise RuntimeError(
                f"raw provenance source mismatch for {row['canonical_event_key']}"
            )
        anchor_sample = seconds_to_sample(
            anchor["anchor_time_seconds"], sfreq, first_samp=first_samp
        )
        epoch_first = anchor_sample + start_offset
        epoch_last = anchor_sample + stop_offset
        if epoch_first < first_samp or epoch_last > last_samp:
            raise RuntimeError(
                f"out-of-bounds {family.name} epoch for {row['canonical_event_key']}: "
                f"[{epoch_first}, {epoch_last}] outside [{first_samp}, {last_samp}]"
            )
        session, block, trial = derive_offset_identifiers(row["offset_trial_count"])
        record = {
            "canonical_order_index": int(row["canonical_order_index"]),
            "canonical_event_key": text_value(row["canonical_event_key"]),
            "ledger_version": text_value(row["ledger_version"]),
            "event_policy_version": text_value(row["event_policy_version"]),
            "behavioural_id": int(row["behavioural_participant_id"]),
            "eeg_source_id": int(row["eeg_recording_id"]),
            "source_recording_filename": text_value(row["source_filename"]),
            "source_file_role": text_value(row["file_role"]),
            "task_filename": text_value(row["task_filename"]),
            "offset_session": session,
            "offset_block": block,
            "offset_trial": trial,
            "audit_trial_count": int(row["audit_trial_count"]),
            "raw_trial_sequence": int(row["raw_trial_sequence"]),
            "join_trial_count": int(row["join_trial_count"]),
            "offset_trial_count": int(row["offset_trial_count"]),
            "physical": bool_value(row["physical"]),
            "condition_semantics": anchor["condition_semantics"],
            "event_source_type": text_value(row["source_type"]),
            "selected_event_source": text_value(row["selected_source"]),
            "trial_event_provenance_json": text_value(row["trial_event_provenance_json"]),
            **anchor,
            "anchor_sample": anchor_sample,
            "epoch_first_sample": epoch_first,
            "epoch_last_sample": epoch_last,
            "epoch_tmin_seconds": family.tmin,
            "epoch_tmax_seconds": family.tmax,
            "epoch_endpoint_inclusive": True,
            "epoch_expected_sample_count": expected_epoch_sample_count(
                family.tmin, family.tmax, sfreq
            ),
            "primary_eligibility": bool_value(row["primary_ledger_eligibility"]),
            "primary_eligibility_reason_code": text_value(
                row["primary_ledger_eligibility_reason_code"]
            ),
            "strict_clean_eligibility": bool_value(
                row["strict_clean_only_sensitivity_eligibility"]
            ),
            "strict_clean_eligibility_reason_code": text_value(
                row["strict_clean_only_sensitivity_eligibility_reason_code"]
            ),
            "duration_warning_flag": bool_value(row["duration_warning_only"]),
            "event_policy_reason_code": text_value(row["event_policy_reason_code"]),
            "accepted_policy_ids": text_value(row["accepted_policy_ids"]),
            "accepted_specific_policy_id": text_value(
                row["accepted_specific_policy_id"]
            ),
            "accepted_policy_action": text_value(row["accepted_policy_action"]),
            "continuous_qc_warning": text_value(row["continuous_qc_warning_codes"]),
            "interpolation_count": int(float(row["interpolated_channel_count"])),
            "interpolation_proportion": float(row["interpolated_channel_proportion"]),
            "interpolation_denominator": int(
                float(row["interpolation_proportion_denominator"])
            ),
            "continuous_v2_run_id": text_value(row["continuous_v2_run_id"]),
            "continuous_v2_manifest_path": text_value(
                row["continuous_v2_manifest_path"]
            ),
            "continuous_v2_manifest_id": text_value(row["continuous_v2_manifest_id"]),
            "post_ica_derivative_kind": text_value(row["continuous_derivative_kind"]),
            "post_ica_derivative_path": text_value(row["continuous_derivative_path"]),
            "continuous_source_sha256": text_value(row["continuous_source_sha256"]),
            "ica_terminal_route": text_value(row["ica_terminal_route"]),
            "post_ica_availability_reason_code": text_value(
                row["post_ica_derivative_availability_reason_code"]
            ),
            "future_readiness_reason_code": text_value(
                row["future_epoch_construction_readiness_reason_code"]
            ),
        }
        records.append(record)

    metadata = pd.DataFrame(records)
    if metadata["canonical_event_key"].duplicated().any():
        raise RuntimeError(f"duplicate canonical keys in {family.name} metadata")
    if metadata["anchor_sample"].duplicated().any():
        duplicated = metadata.loc[
            metadata["anchor_sample"].duplicated(keep=False),
            ["canonical_event_key", "anchor_sample"],
        ]
        raise RuntimeError(
            f"duplicate MNE event samples in {family.name}: {duplicated.to_dict(orient='records')}"
        )
    return metadata


def events_from_metadata(metadata: pd.DataFrame, event_id: int) -> np.ndarray:
    """Return an MNE event array in the exact metadata row order."""

    if event_id <= 0:
        raise ValueError("MNE event_id must be positive")
    samples = metadata["anchor_sample"].to_numpy(dtype=np.int64)
    return np.column_stack(
        [samples, np.zeros(len(samples), dtype=np.int64), np.full(len(samples), event_id, dtype=np.int64)]
    )


def create_epochs(
    raw: mne.io.BaseRaw,
    metadata: pd.DataFrame,
    family: EpochFamily,
) -> mne.Epochs:
    """Create one preloaded MNE Epochs shard without signal alteration.

    ``reject_by_annotation=False`` is deliberate: accepted ledger rows, not
    annotations carried for continuous provenance, own epoch eligibility.
    No reject/flat thresholds, baseline, detrending, projection, decimation,
    or resampling are applied.
    """

    if not np.isclose(raw.info["sfreq"], EXPECTED_SFREQ, atol=1e-9, rtol=0.0):
        raise RuntimeError(f"continuous derivative sfreq is {raw.info['sfreq']}, expected 1000 Hz")
    events = events_from_metadata(metadata, family.event_id)
    epochs = mne.Epochs(
        raw,
        events,
        event_id={family.name: family.event_id},
        tmin=family.tmin,
        tmax=family.tmax,
        baseline=None,
        picks=None,
        preload=True,
        reject=None,
        flat=None,
        proj=False,
        detrend=None,
        decim=1,
        reject_by_annotation=False,
        event_repeated="error",
        on_missing="raise",
        metadata=metadata,
        verbose="ERROR",
    )
    validate_epochs_object(epochs, metadata, family)
    return epochs


def validate_epochs_object(
    epochs: mne.BaseEpochs,
    expected_metadata: pd.DataFrame,
    family: EpochFamily,
) -> dict[str, Any]:
    """Validate a constructed or reopened Epochs object against the contract."""

    expected_n_times = expected_epoch_sample_count(
        family.tmin, family.tmax, EXPECTED_SFREQ
    )
    problems: list[str] = []
    if len(epochs) != len(expected_metadata):
        problems.append(f"epoch_count={len(epochs)} expected={len(expected_metadata)}")
    if epochs.metadata is None:
        problems.append("metadata_missing")
        observed_keys: list[str] = []
    else:
        observed_keys = epochs.metadata["canonical_event_key"].astype(str).tolist()
        expected_keys = expected_metadata["canonical_event_key"].astype(str).tolist()
        if observed_keys != expected_keys:
            problems.append("canonical_key_order_mismatch")
    expected_samples = expected_metadata["anchor_sample"].to_numpy(dtype=np.int64)
    if not np.array_equal(epochs.events[:, 0], expected_samples):
        problems.append("event_sample_order_mismatch")
    if not np.isclose(epochs.info["sfreq"], EXPECTED_SFREQ, atol=1e-9, rtol=0.0):
        problems.append(f"sfreq={epochs.info['sfreq']}")
    if epochs.baseline is not None:
        problems.append(f"baseline={epochs.baseline}")
    if epochs.reject is not None or epochs.flat is not None:
        problems.append("reject_or_flat_present")
    if epochs.detrend is not None:
        problems.append(f"detrend={epochs.detrend}")
    if len(epochs.times) != expected_n_times:
        problems.append(f"n_times={len(epochs.times)} expected={expected_n_times}")
    if not np.isclose(epochs.times[0], family.tmin, atol=1e-12, rtol=0.0):
        problems.append(f"first_time={epochs.times[0]}")
    if not np.isclose(epochs.times[-1], family.tmax, atol=1e-12, rtol=0.0):
        problems.append(f"last_time={epochs.times[-1]}")
    if len(epochs.selection) != len(expected_metadata):
        problems.append("selection_indicates_drop")
    if any(entry for row in epochs.drop_log for entry in row):
        problems.append("drop_log_not_empty")
    data = epochs.get_data(copy=False)
    if data.shape[2] != expected_n_times:
        problems.append(f"data_n_times={data.shape[2]}")
    if not np.isfinite(data).all():
        problems.append("nonfinite_epoch_data")
    if problems:
        raise RuntimeError(f"{family.name} Epochs validation failed: {problems}")
    return {
        "valid": True,
        "family": family.name,
        "epoch_count": len(epochs),
        "strict_clean_count": int(expected_metadata["strict_clean_eligibility"].sum()),
        "duration_warning_count": int(expected_metadata["duration_warning_flag"].sum()),
        "sampling_frequency_hz": float(epochs.info["sfreq"]),
        "sample_count_per_epoch": expected_n_times,
        "tmin_seconds": float(epochs.times[0]),
        "tmax_seconds": float(epochs.times[-1]),
        "endpoint_inclusive": True,
        "baseline": None,
        "reject": None,
        "flat": None,
        "detrend": None,
        "decimation": 1,
        "channel_count": len(epochs.ch_names),
        "canonical_key_sha256": content_fingerprint(observed_keys),
        "metadata_content_sha256": stable_frame_hash(expected_metadata),
        "nonfinite_sample_count": 0,
        "drop_count": 0,
    }


def write_epochs_atomic(epochs: mne.BaseEpochs, path: Path) -> dict[str, Any]:
    """Save one bounded Epochs shard through a sibling temporary FIF.

    Epoch data use MNE's explicit single-precision FIF storage.  The driver
    shards by recording, so a shard is not expected to split; an unexpected
    split fails closed instead of publishing an incomplete artifact set.
    """

    if not path.name.endswith("-epo.fif"):
        raise ValueError("Epoch FIF names must end in -epo.fif")
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing epoch artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}-epo.fif")
    epochs.save(
        temporary,
        fmt="single",
        overwrite=False,
        split_size="2GB",
        split_naming="neuromag",
        verbose="ERROR",
    )
    split_candidates = list(temporary.parent.glob(f"{temporary.stem}-*.fif"))
    if split_candidates:
        raise RuntimeError(f"unexpected split Epochs shard: {split_candidates}")
    os.replace(temporary, path)
    return {
        "path": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "fif_format": "single",
    }


def reopen_and_validate_epochs(
    path: Path,
    expected_metadata: pd.DataFrame,
    family: EpochFamily,
    *,
    expected_sha256: str | None = None,
) -> dict[str, Any]:
    """Reopen one saved Epochs FIF with preload and validate content/hash."""

    observed_sha256 = sha256_file(path)
    if expected_sha256 and observed_sha256 != expected_sha256:
        raise RuntimeError(f"saved Epochs checksum mismatch: {path}")
    reopened = mne.read_epochs(path, preload=True, verbose="ERROR")
    validation = validate_epochs_object(reopened, expected_metadata, family)
    validation.update(
        {
            "path": path.as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": observed_sha256,
            "preloaded": bool(reopened.preload),
        }
    )
    return validation


def ordered_key_hash(keys: Iterable[str]) -> str:
    """Hash an ordered canonical-key sequence deterministically."""

    return content_fingerprint([str(key) for key in keys])


def assert_identical_family_keys(family_metadata: Mapping[str, pd.DataFrame]) -> str:
    """Require identical ordered canonical keys across all accepted families."""

    expected_names = [family.name for family in FAMILIES]
    if set(family_metadata) != set(expected_names):
        raise RuntimeError(
            f"family metadata set differs: observed={sorted(family_metadata)}, expected={expected_names}"
        )
    reference = family_metadata[expected_names[0]]["canonical_event_key"].astype(str).tolist()
    for name in expected_names[1:]:
        observed = family_metadata[name]["canonical_event_key"].astype(str).tolist()
        if observed != reference:
            raise RuntimeError(f"ordered canonical keys differ for {name}")
    return ordered_key_hash(reference)


def forbidden_product_scan(output_root: Path) -> list[str]:
    """Return unauthorized analysis products found in the epoch namespace."""

    forbidden_fragments = (
        "autoreject",
        "csd",
        "tfr",
        "time_frequency",
        "time-frequency",
        "band_power",
        "spectral_normalization",
        "id86",
        "pre_ica",
        "pre-ica",
    )
    found: list[str] = []
    if not output_root.exists():
        return found
    for path in output_root.rglob("*"):
        if path.is_file() and any(fragment in path.name.lower() for fragment in forbidden_fragments):
            found.append(path.relative_to(output_root).as_posix())
    return sorted(found)


def source_snapshot(paths: Sequence[Path]) -> list[dict[str, Any]]:
    """Capture non-mutating size and modification-time evidence for sources."""

    rows = []
    for path in sorted({item.resolve() for item in paths}):
        stat = path.stat()
        rows.append(
            {
                "path": path.as_posix(),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
    return rows


def compare_source_snapshots(
    before: Sequence[Mapping[str, Any]], after: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """Compare source path/size/mtime snapshots and fail on any change."""

    normalized_before = [dict(row) for row in before]
    normalized_after = [dict(row) for row in after]
    unchanged = normalized_before == normalized_after
    if not unchanged:
        raise RuntimeError("linked continuous source FIF size/mtime surface changed during epoch construction")
    return {
        "unchanged": True,
        "file_count": len(normalized_before),
        "total_bytes": sum(int(row["size_bytes"]) for row in normalized_before),
    }
