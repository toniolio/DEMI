"""Contracts and pure helpers for the DEMI event/epoch eligibility ledger.

This module supports the accepted-policy ledger stage without reading EEG
signals or constructing MNE epochs.  It keeps the high-risk distinctions
between event-policy candidacy, primary and strict-clean event surfaces,
continuous-derivative availability, post-ICA availability, and future epoch
construction readiness explicit and testable.

The numbered ledger driver owns local file discovery and output writing.  The
helpers here are deliberately small and deterministic so public tests can
exercise the contracts without access to private policy records, local EEG
files, or generated evidence.

Inputs:
    In-memory row dictionaries, pandas data frames, and parsed continuous-v2
    recording manifests supplied by the caller.

Outputs:
    Pure status dictionaries, canonical keys, parsed continuous provenance,
    validation results, and deterministic table-content hashes.

This module explicitly does not open EDF/FIF files, mutate events, repair
trials, decide participant or analytic inclusion, or create epochs or any
signal derivative.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


LEDGER_VERSION = "event_epoch_eligibility_v1"
POLICY_VERSION = "v1.0"
EXPECTED_PRIMARY_COUNT = 8_905
EXPECTED_STRICT_CLEAN_COUNT = 8_896
EXPECTED_WARNING_ONLY_COUNT = 9
EXPECTED_POSSIBLE_DISCREPANCY_COUNT = 72
SPLIT_EEG_IDS = frozenset({54, 56, 65})

FORBIDDEN_SIGNAL_OUTPUT_SUFFIXES = (
    ".edf",
    ".bdf",
    ".fif",
    ".set",
    ".vhdr",
    ".vmrk",
)
FORBIDDEN_OUTPUT_NAME_FRAGMENTS = (
    "autoreject",
    "csd",
    "time_frequency",
    "time-frequency",
)


def bool_value(value: Any) -> bool:
    """Return a strict, missing-safe Boolean interpretation.

    Args:
        value: Scalar value read from CSV, JSON, or a pandas row.

    Returns:
        ``True`` only for ordinary true-like values; missing and false-like
        strings return ``False``.

    Side effects:
        None.
    """

    if value is None or value is pd.NA:
        return False
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)


def nullable_int(value: Any) -> int | None:
    """Convert a scalar to an integer while preserving missing values.

    Args:
        value: Numeric, string, or missing scalar.

    Returns:
        Integer value, or ``None`` when the input is missing or blank.

    Side effects:
        None.
    """

    if value is None or value is pd.NA:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return int(float(value))


def text_value(value: Any) -> str:
    """Return stripped text while treating pandas missing values as blank.

    Args:
        value: Scalar from a CSV/JSON/pandas row.

    Returns:
        Stripped string, or an empty string for missing values.

    Side effects:
        None.
    """

    if value is None or value is pd.NA:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def status_fields(
    prefix: str,
    available: bool,
    true_code: str,
    true_reason: str,
    false_code: str,
    false_reason: str,
) -> dict[str, Any]:
    """Build one explicit Boolean/status/reason field family.

    Args:
        prefix: Output column prefix.
        available: Boolean result for the status.
        true_code: Machine-readable code when ``available`` is true.
        true_reason: Readable explanation when ``available`` is true.
        false_code: Machine-readable code when ``available`` is false.
        false_reason: Readable explanation when ``available`` is false.

    Returns:
        Four fields: the Boolean, status string, reason code, and readable
        reason.

    Side effects:
        None.
    """

    return {
        prefix: bool(available),
        f"{prefix}_status": "available" if available else "unavailable",
        f"{prefix}_reason_code": true_code if available else false_code,
        f"{prefix}_reason": true_reason if available else false_reason,
    }


def canonical_event_key(row: Mapping[str, Any]) -> str:
    """Build a stable canonical key for a corrected join row.

    Raw-present rows retain the native source filename, reconstructed raw trial
    sequence, and selected anchor sample.  Offset-only rows use the frozen
    behavioural participant/trial key and are explicitly labelled as lacking
    a raw event.  A short digest prevents punctuation differences from making
    downstream identifiers unwieldy while the readable basis remains in a
    separate ledger column.

    Args:
        row: Mapping containing join, identity, source, and anchor fields.

    Returns:
        Versioned canonical event key.

    Side effects:
        None.
    """

    source_filename = text_value(row.get("source_filename"))
    behavioural_id = nullable_int(row.get("behavioural_participant_id"))
    eeg_id = nullable_int(row.get("eeg_recording_id"))
    raw_trial = nullable_int(row.get("raw_trial_sequence"))
    offset_trial = nullable_int(row.get("offset_trial_count"))
    anchor_sample = nullable_int(row.get("anchor_sample"))

    if source_filename:
        basis = (
            f"raw|eeg={eeg_id}|beh={behavioural_id}|file={source_filename}|"
            f"raw_trial={raw_trial}|sample={anchor_sample}"
        )
    else:
        basis = (
            f"offset_only|eeg={eeg_id}|beh={behavioural_id}|"
            f"offset_trial={offset_trial}"
        )
    digest = hashlib.sha256(basis.encode("utf-8")).hexdigest()[:20]
    return f"{LEDGER_VERSION}:{digest}"


def parse_continuous_recording_manifest(
    manifest: Mapping[str, Any],
    manifest_path: Path,
    repo_root: Path,
) -> dict[str, Any]:
    """Extract ledger-safe facts from one continuous-v2 recording manifest.

    Args:
        manifest: Parsed per-recording continuous manifest.
        manifest_path: Path from which the manifest was read.
        repo_root: Repository root used to render stable relative paths.

    Returns:
        Flat dictionary containing terminal status, canonical derivative kind
        and path, QC warning codes, interpolation facts, and ICA route.

    Side effects:
        None.  The referenced FIF and ICA artifacts are not opened.
    """

    artifacts = {
        str(item.get("kind")): item
        for item in manifest.get("artifacts", [])
        if isinstance(item, Mapping)
    }
    recording_dir = manifest_path.parent
    post_artifact = artifacts.get("post_ica_continuous")
    pre_artifact = artifacts.get("pre_ica_continuous")
    canonical_artifact = post_artifact or pre_artifact
    derivative_kind = (
        "post_ica"
        if post_artifact is not None
        else "retained_pre_ica"
        if pre_artifact is not None
        else "none"
    )
    derivative_path = ""
    if canonical_artifact is not None:
        derivative_path = (recording_dir / str(canonical_artifact["relative_path"])).resolve().relative_to(
            repo_root.resolve()
        ).as_posix()

    qc_rows = [
        item for item in manifest.get("qc_warnings", []) if isinstance(item, Mapping)
    ]
    qc_codes = sorted({str(item.get("code")) for item in qc_rows if item.get("code")})
    interpolation = manifest.get("interpolation", {})
    ica = manifest.get("ica", {})
    proposal = ica.get("proposal", {})
    application = ica.get("application", {})
    stop = manifest.get("stop_or_failure", {})

    if derivative_kind == "post_ica":
        ica_route = "ordinary_post_ica_complete"
    elif (
        derivative_kind == "retained_pre_ica"
        and str(stop.get("code")) == "predeclared_id86_component_review_boundary"
    ):
        ica_route = "accepted_historical_zero_component_stop_pre_ica_retained"
    elif derivative_kind == "retained_pre_ica":
        ica_route = "pre_ica_review_stop"
    else:
        ica_route = "no_continuous_derivative"

    source = manifest.get("source", {})
    return {
        "source_filename": str(source.get("source_filename", "")),
        "continuous_v2_manifest_path": manifest_path.resolve().relative_to(
            repo_root.resolve()
        ).as_posix(),
        "continuous_v2_manifest_id": str(manifest.get("manifest_id", "")),
        "continuous_v2_terminal_status": str(manifest.get("status", "")),
        "continuous_v2_completion_class": str(manifest.get("completion_class", "")),
        "continuous_derivative_kind": derivative_kind,
        "continuous_derivative_path": derivative_path,
        "continuous_qc_warning_codes": ";".join(qc_codes),
        "continuous_qc_warnings_json": json.dumps(
            qc_rows, sort_keys=True, separators=(",", ":")
        ),
        "interpolated_channel_count": nullable_int(
            interpolation.get("interpolation_count")
        ),
        "interpolated_channel_proportion": interpolation.get(
            "interpolation_proportion"
        ),
        "interpolation_proportion_denominator": nullable_int(
            interpolation.get("proportion_denominator")
        ),
        "ica_terminal_route": ica_route,
        "ica_proposal_count": nullable_int(proposal.get("proposal_count")),
        "ica_exclusion_count": nullable_int(application.get("exclusion_count")),
        "ica_application_performed": bool_value(
            application.get("ica_application_performed")
        ),
        "continuous_stop_reason_code": str(stop.get("code", "")),
        "continuous_stop_reason": str(stop.get("detail", "")),
        "continuous_source_sha256": str(source.get("source_sha256", "")),
    }


def stable_frame_hash(frame: pd.DataFrame) -> str:
    """Hash table content deterministically after stable column normalization.

    Args:
        frame: Data frame whose current row order is part of the contract.

    Returns:
        SHA-256 digest of deterministic CSV serialization.

    Side effects:
        None.
    """

    normalized = frame.copy()
    normalized = normalized.reindex(sorted(normalized.columns), axis=1)
    payload = normalized.to_csv(index=False, lineterminator="\n", na_rep="").encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def validate_ledger_contract(frame: pd.DataFrame) -> dict[str, int]:
    """Validate the accepted v1.0 count, key, and reason-code contracts.

    Args:
        frame: Completed canonical eligibility ledger.

    Returns:
        Exact primary, strict-clean, warning-only, and discrepancy counts.

    Raises:
        RuntimeError: If a count, uniqueness, implication, or explicit-reason
            contract fails.

    Side effects:
        None.
    """

    required_booleans = (
        "accepted_event_policy_candidate",
        "primary_ledger_eligibility",
        "strict_clean_only_sensitivity_eligibility",
        "continuous_derivative_availability",
        "post_ica_derivative_availability",
        "future_epoch_construction_readiness",
    )
    for column in required_booleans:
        if column not in frame.columns:
            raise RuntimeError(f"ledger is missing required status column {column}")
        false_rows = ~frame[column].map(bool_value)
        for suffix in ("_reason_code", "_reason"):
            reason_column = f"{column}{suffix}"
            if reason_column not in frame.columns:
                raise RuntimeError(f"ledger is missing {reason_column}")
            missing_reason = false_rows & frame[reason_column].fillna("").astype(str).str.strip().eq("")
            if missing_reason.any():
                raise RuntimeError(
                    f"{int(missing_reason.sum())} false {column} rows lack {reason_column}"
                )

    if frame["canonical_event_key"].duplicated().any():
        raise RuntimeError("canonical_event_key contains duplicates")

    primary = frame["primary_ledger_eligibility"].map(bool_value)
    strict = frame["strict_clean_only_sensitivity_eligibility"].map(bool_value)
    warning_only = primary & ~strict
    discrepancy = frame["event_policy_reason_code"].eq(
        "possible_event_sequence_discrepancy_unavailable"
    )
    counts = {
        "primary": int(primary.sum()),
        "strict_clean": int(strict.sum()),
        "warning_only": int(warning_only.sum()),
        "possible_discrepancy": int(discrepancy.sum()),
    }
    expected = {
        "primary": EXPECTED_PRIMARY_COUNT,
        "strict_clean": EXPECTED_STRICT_CLEAN_COUNT,
        "warning_only": EXPECTED_WARNING_ONLY_COUNT,
        "possible_discrepancy": EXPECTED_POSSIBLE_DISCREPANCY_COUNT,
    }
    if counts != expected:
        raise RuntimeError(f"accepted ledger counts differ: observed={counts}, expected={expected}")
    if (strict & ~primary).any():
        raise RuntimeError("strict-clean eligibility must imply primary eligibility")
    if not frame.loc[warning_only, "duration_warning_only"].map(bool_value).all():
        raise RuntimeError("every primary-only row must be an explicit duration warning")
    return counts


def assert_no_signal_outputs(output_dir: Path) -> None:
    """Fail if the ledger namespace contains an epoch or signal derivative.

    Args:
        output_dir: Ledger output namespace to scan recursively.

    Returns:
        None.

    Raises:
        RuntimeError: If a forbidden signal suffix or stage-name fragment is
            present.

    Side effects:
        Reads directory entries only.
    """

    if not output_dir.exists():
        return
    forbidden: list[str] = []
    for path in output_dir.rglob("*"):
        if not path.is_file():
            continue
        lower_name = path.name.lower()
        if lower_name.endswith(FORBIDDEN_SIGNAL_OUTPUT_SUFFIXES) or any(
            fragment in lower_name for fragment in FORBIDDEN_OUTPUT_NAME_FRAGMENTS
        ):
            forbidden.append(path.as_posix())
    if forbidden:
        raise RuntimeError(f"ledger namespace contains forbidden signal outputs: {forbidden}")
