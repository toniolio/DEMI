"""Deterministic validation-cohort selection from existing technical evidence.

The selector uses the already generated detector-audit group labels and
30-scalp published-criterion counts to choose a compact architecture-validation
cohort. It includes two ordinary controls, both line-noise comparison files, a
high-count accepted bad-channel/interpolation case not otherwise selected, and
the predeclared ID-86 ICA review boundary.

Selection is based on technical QC/audit roles, never scientific outcomes. The
module does not decide participant inclusion, event/epoch eligibility, or
whether a recording belongs in an analysis. It reads CSV evidence and raw-file
names only; it writes nothing.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .storage import parse_recording_id
from .contracts import sha256_file


def _selection_groups(value: str) -> set[str]:
    """Parse one deterministic JSON list of audit selection groups."""

    parsed = json.loads(value)
    if not isinstance(parsed, list) or not all(isinstance(item, str) for item in parsed):
        raise ValueError(f"Invalid selection_groups_json value: {value!r}")
    return set(parsed)


def _add_reason(reasons: dict[str, list[str]], filename: str, reason: str) -> None:
    """Add one stable selection reason without duplication."""

    reasons.setdefault(filename, [])
    if reason not in reasons[filename]:
        reasons[filename].append(reason)


def select_validation_cohort(
    *,
    raw_root: Path,
    audit_runs_path: Path,
    sensitivity_path: Path,
    event_inventory_path: Path | None = None,
) -> dict[str, Any]:
    """Select and document the fixed six-recording validation cohort.

    Args:
        raw_root: Local immutable EDF directory.
        audit_runs_path: Existing global-method audit run table.
        sensitivity_path: Existing 30-versus-32 detector-pool evidence table.
        event_inventory_path: Optional factual event-source inventory used only
            to state availability; it does not gate continuous preprocessing.

    Returns:
        Exact ordered cohort members, reasons, evidence hashes, and selection
        rules suitable for a generated run-level manifest.

    Side effects:
        Reads local CSV evidence and checks raw-file existence. Writes nothing.
    """

    runs = pd.read_csv(audit_runs_path)
    required_run_columns = {"source_filename", "selection_groups_json", "method", "repeat"}
    if not required_run_columns.issubset(runs.columns):
        raise ValueError(f"Audit run table lacks columns: {required_run_columns - set(runs.columns)}")
    primary_runs = runs[
        runs["method"].eq("pyprep_noisy_channels") & runs["repeat"].eq(1)
    ].drop_duplicates("source_filename")
    groups_by_file = {
        row.source_filename: _selection_groups(row.selection_groups_json)
        for row in primary_runs.itertuples(index=False)
    }

    reasons: dict[str, list[str]] = {}
    ordered: list[str] = []

    ordinary = sorted(
        filename
        for filename, groups in groups_by_file.items()
        if "ordinary_non_candidate" in groups
    )
    if len(ordinary) < 2:
        raise ValueError("Audit evidence does not contain two ordinary non-candidate files.")
    for filename in ordinary[:2]:
        ordered.append(filename)
        _add_reason(reasons, filename, "ordinary_non_candidate_control")
    _add_reason(reasons, ordinary[0], "representative_accepted_global_bad_handling_evidence")
    _add_reason(reasons, ordinary[1], "ordinary_zero_primary_bad_audit_evidence")

    line_files = sorted(
        filename
        for filename, groups in groups_by_file.items()
        if "line_noise_comparison" in groups
    )
    if line_files != ["demi_28 Data.edf", "demi_29 Data.edf"]:
        raise ValueError(f"Expected deterministic line-noise files 28 and 29; found {line_files}")
    for filename in line_files:
        ordered.append(filename)
        _add_reason(reasons, filename, "high_60_hz_line_noise_comparison_evidence")
        _add_reason(reasons, filename, "multiple_suspicious_scalp_channel_audit_context")

    sensitivity = pd.read_csv(sensitivity_path)
    required_sensitivity = {
        "source_filename",
        "method",
        "sensitivity_30ch_policy_count",
        "sensitivity_30ch_policy_calls_json",
    }
    if not required_sensitivity.issubset(sensitivity.columns):
        raise ValueError(
            f"M1/M2 sensitivity table lacks columns: {required_sensitivity - set(sensitivity.columns)}"
        )
    candidates = sensitivity[
        sensitivity["method"].eq("pyprep_noisy_channels")
        & ~sensitivity["source_filename"].isin(ordered)
        & sensitivity["source_filename"].map(parse_recording_id).ne(86)
    ].copy()
    if candidates.empty:
        raise ValueError("No independent 30-channel bad/interpolation candidate is available.")
    candidates = candidates.sort_values(
        ["sensitivity_30ch_policy_count", "source_filename"],
        ascending=[False, True],
        kind="stable",
    )
    representative = str(candidates.iloc[0]["source_filename"])
    ordered.append(representative)
    _add_reason(reasons, representative, "highest_accepted_30_scalp_published_criterion_count_not_already_selected")
    _add_reason(reasons, representative, "representative_multi_channel_reference_and_interpolation_case")

    id86_matches = sorted(
        filename for filename in groups_by_file if parse_recording_id(filename) == 86
    )
    if len(id86_matches) != 1:
        raise ValueError(f"Expected exactly one ID-86 audit recording; found {id86_matches}")
    ordered.append(id86_matches[0])
    _add_reason(reasons, id86_matches[0], "predeclared_ica_component_review_boundary")
    _add_reason(reasons, id86_matches[0], "intentional_stop_route_validation")

    if not 5 <= len(ordered) <= 7 or len(ordered) != len(set(ordered)):
        raise ValueError(f"Validation cohort is not a unique compact 5--7 file set: {ordered}")

    sensitivity_by_file = sensitivity[
        sensitivity["method"].eq("pyprep_noisy_channels")
    ].set_index("source_filename")
    event_by_file: dict[str, dict[str, Any]] = {}
    if event_inventory_path is not None and event_inventory_path.is_file():
        event_rows = pd.read_csv(event_inventory_path).drop_duplicates("source_filename")
        event_by_file = {
            str(row.source_filename): {
                "event_source_availability": (
                    "selected_coherent_source_available"
                    if "selected" in str(row.selection_status).lower()
                    and pd.notna(row.selected_source)
                    else f"event_selection_status:{row.selection_status}"
                ),
                "selected_event_source": row.selected_source,
                "event_selection_status": row.selection_status,
                "event_selection_reason": row.selection_reason,
            }
            for row in event_rows.itertuples(index=False)
        }

    members: list[dict[str, Any]] = []
    for order, filename in enumerate(ordered, start=1):
        path = raw_root / filename
        if not path.is_file():
            raise FileNotFoundError(f"Selected validation EDF is missing: {path}")
        sensitivity_row = sensitivity_by_file.loc[filename]
        members.append(
            {
                "order": order,
                "source_filename": filename,
                "recording_id": parse_recording_id(filename),
                "source_path": path.as_posix(),
                "selection_reasons": reasons[filename],
                "audit_selection_groups": sorted(groups_by_file[filename]),
                "audit_30ch_published_criterion_count": int(
                    sensitivity_row["sensitivity_30ch_policy_count"]
                ),
                "audit_30ch_published_criterion_calls": json.loads(
                    sensitivity_row["sensitivity_30ch_policy_calls_json"]
                ),
                **event_by_file.get(
                    filename,
                    {"event_source_availability": "event_inventory_not_available"},
                ),
                "event_or_epoch_status_used_for_selection": False,
                "scientific_outcome_used_for_selection": False,
            }
        )

    evidence_paths = [audit_runs_path, sensitivity_path]
    if event_inventory_path is not None and event_inventory_path.is_file():
        evidence_paths.append(event_inventory_path)
    return {
        "schema_version": 1,
        "cohort_kind": "production_architecture_validation_not_full_dataset",
        "selection_algorithm": [
            "first_two_lexicographic_ordinary_non_candidate_audit_files",
            "both_fixed_line_noise_comparison_files_28_and_29",
            "highest_30_scalp_published_criterion_count_not_already_selected_with_filename_tiebreak",
            "sole_predeclared_id86_component_review_exception",
        ],
        "selection_uses_scientific_outcomes": False,
        "selection_changes_participant_inclusion": False,
        "selection_changes_event_or_epoch_eligibility": False,
        "cohort_size": len(members),
        "members": members,
        "evidence_files": [
            {
                "path": path.as_posix(),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in evidence_paths
        ],
    }
