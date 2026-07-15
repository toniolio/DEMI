"""Sequential continuous-preprocessing execution and aggregate run evidence.

The runner selects either the fixed technical validation cohort or the
explicitly authorized 95-file readable production surface, hashes code/config/
environment state, processes one EDF at a time, and flushes aggregate progress
after every recording. Each mode has a separate fixed ignored output root.

``--max-files`` supports bounded smoke/resume checks; valid completed/stopped
recording directories are then hash-validated and skipped. The runner never
constructs epochs, runs AutoReject, computes CSD, or changes event/participant
eligibility.
"""

from __future__ import annotations

import json
import shutil
import statistics
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from .cohort import select_production_surface, select_validation_cohort
from .contracts import canonical_json_bytes, load_config, sha256_bytes, sha256_file
from .pipeline import process_recording
from .repair import (
    REPAIR_KIND,
    compare_completed_manifest_to_historical_rule,
    is_historical_ica_repair_target,
    repair_recording_from_saved_ica,
)
from .storage import (
    atomic_write_json,
    atomic_write_text,
    run_provenance,
    validate_derivative_path,
    validate_saved_ica,
    validate_saved_raw,
    validate_output_root,
)


SURFACE_MODES = ("validation", "production", "production_v2")


def production_code_paths(repo_root: Path) -> list[Path]:
    """Return the exact code files that invalidate saved recording caches."""

    package = repo_root / "analysis" / "eeg_mne" / "continuous_preprocessing"
    paths = sorted(package.glob("*.py"))
    paths.extend(
        [
            repo_root / "analysis" / "eeg_mne" / "13_run_continuous_preprocessing_validation.py",
            repo_root / "requirements-eeg.txt",
        ]
    )
    return paths


def cohort_evidence_paths(repo_root: Path) -> dict[str, Path]:
    """Return local deterministic audit evidence paths used by selection."""

    audit_root = (
        repo_root
        / "_Data"
        / "eeg"
        / "mne_preprocessing"
        / "preprocessing_parameter_audit_v1"
    )
    return {
        "audit_runs": audit_root / "global_bad_method_runs.csv",
        "sensitivity": audit_root / "m1_m2_detector_sensitivity.csv",
        "event_inventory": repo_root
        / "_Data"
        / "eeg"
        / "event_source_inventory_v1"
        / "event_source_file_inventory.csv",
    }


def build_validation_cohort(repo_root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    """Build the exact cohort selection package without processing signals."""

    paths = cohort_evidence_paths(repo_root)
    raw_root = repo_root / config["paths"]["raw_root"]
    return select_validation_cohort(
        raw_root=raw_root,
        audit_runs_path=paths["audit_runs"],
        sensitivity_path=paths["sensitivity"],
        event_inventory_path=paths["event_inventory"],
    )


def build_production_surface(repo_root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    """Build the exact authorized readable-EDF production surface."""

    raw_root = (repo_root / config["paths"]["raw_root"]).resolve()
    manifest_path = (repo_root / config["production_surface"]["raw_manifest"]).resolve()
    return select_production_surface(
        raw_root=raw_root,
        raw_manifest_path=manifest_path,
        expected_readable_count=int(
            config["production_surface"]["expected_readable_edf_count"]
        ),
    )


def build_surface(
    repo_root: Path, config: Mapping[str, Any], surface_mode: str
) -> dict[str, Any]:
    """Return the selected validation or authorized production surface."""

    if surface_mode == "validation":
        cohort = build_validation_cohort(repo_root, config)
        return {
            **cohort,
            "surface_kind": cohort["cohort_kind"],
            "surface_size": cohort["cohort_size"],
        }
    if surface_mode in {"production", "production_v2"}:
        return build_production_surface(repo_root, config)
    raise ValueError(f"Unknown continuous-preprocessing surface mode: {surface_mode!r}")


def output_root_for_mode(
    repo_root: Path, config: Mapping[str, Any], surface_mode: str
) -> tuple[Path, str]:
    """Resolve and fail closed on the configured output namespace for a mode."""

    if surface_mode == "validation":
        configured = str(config["paths"]["output_root"])
    elif surface_mode == "production":
        configured = str(config["paths"]["production_output_root"])
    elif surface_mode == "production_v2":
        configured = str(config["paths"]["production_v2_output_root"])
    else:
        raise ValueError(f"Unknown continuous-preprocessing surface mode: {surface_mode!r}")
    output_root = (repo_root / configured).resolve()
    validate_output_root(repo_root, output_root, configured)
    return output_root, configured


def validate_production_preflight(
    repo_root: Path, config: Mapping[str, Any], output_root: Path
) -> dict[str, Any]:
    """Verify disk, ignore, and separate-root requirements before production.

    This check writes nothing. It deliberately runs before even a bounded
    production smoke so an unattended invocation cannot drift into an
    unignored, overlapping, or undersized destination.
    """

    if config["line_noise"]["enabled"] is not False:
        raise ValueError("The authorized production run requires line-noise removal off.")
    if config["safety"]["parallel_recordings"] is not False or config["safety"][
        "process_one_recording_at_a_time"
    ] is not True:
        raise ValueError("The authorized production run requires recording concurrency one.")
    validation_root = (repo_root / config["paths"]["output_root"]).resolve()
    production_v1_root = (repo_root / config["paths"]["production_output_root"]).resolve()
    production_root = output_root.resolve()
    protected_roots = {validation_root, production_v1_root}
    protected_roots.discard(production_root)
    overlap = any(
        production_root == other
        or production_root in other.parents
        or other in production_root.parents
        for other in protected_roots
    )
    if overlap:
        raise ValueError("Production and validation output roots overlap.")
    ignore = subprocess.run(
        ["git", "check-ignore", "--quiet", production_root.as_posix()],
        cwd=repo_root,
        check=False,
    )
    if ignore.returncode != 0:
        raise ValueError(f"Production output root is not ignored: {production_root}")
    disk = shutil.disk_usage(production_root.parent)
    minimum = int(config["production_surface"]["minimum_free_bytes"])
    if disk.free < minimum:
        raise RuntimeError(
            f"Production requires at least {minimum} free bytes; found {disk.free}."
        )
    return {
        "production_output_root": production_root.as_posix(),
        "validation_output_root": validation_root.as_posix(),
        "preserved_production_v1_output_root": production_v1_root.as_posix(),
        "roots_are_separate": True,
        "production_output_root_is_gitignored": True,
        "minimum_free_bytes": minimum,
        "free_bytes_before_invocation": disk.free,
        "free_space_requirement_passed": True,
        "recording_concurrency": 1,
    }


def snapshot_source_inventory(members: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Hash one immutable source surface with size and nanosecond mtime evidence.

    Args:
        members: Complete selected surface, not merely a bounded invocation.

    Returns:
        Stable per-file SHA-256, size, and timestamp rows.

    Side effects:
        Reads each EDF once and writes nothing.
    """

    rows: list[dict[str, Any]] = []
    for member in members:
        path = Path(member["source_path"])
        before = path.stat()
        digest = sha256_file(path)
        after = path.stat()
        if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
            raise RuntimeError(f"Raw EDF changed while source inventory was hashed: {path}")
        rows.append(
            {
                "order": int(member["order"]),
                "source_filename": str(member["source_filename"]),
                "source_path": path.as_posix(),
                "sha256": digest,
                "size_bytes": before.st_size,
                "mtime_ns": before.st_mtime_ns,
            }
        )
    return {
        "file_count": len(rows),
        "total_size_bytes": sum(row["size_bytes"] for row in rows),
        "rows": rows,
    }


def compare_source_inventories(
    before: Mapping[str, Any], after: Mapping[str, Any]
) -> dict[str, Any]:
    """Compare two complete raw inventories without weakening identity checks."""

    before_rows = {row["source_filename"]: row for row in before["rows"]}
    after_rows = {row["source_filename"]: row for row in after["rows"]}
    missing = sorted(set(before_rows) - set(after_rows))
    added = sorted(set(after_rows) - set(before_rows))
    changed: list[dict[str, Any]] = []
    for filename in sorted(set(before_rows) & set(after_rows)):
        first = before_rows[filename]
        second = after_rows[filename]
        differences = [
            field
            for field in ("sha256", "size_bytes", "mtime_ns", "source_path")
            if first[field] != second[field]
        ]
        if differences:
            changed.append({"source_filename": filename, "changed_fields": differences})
    unchanged = not missing and not added and not changed
    return {
        "unchanged": unchanged,
        "before_file_count": before["file_count"],
        "after_file_count": after["file_count"],
        "before_total_size_bytes": before["total_size_bytes"],
        "after_total_size_bytes": after["total_size_bytes"],
        "missing_files": missing,
        "added_files": added,
        "changed_files": changed,
    }


def _load_terminal_manifests(output_root: Path, members: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Read currently published per-recording manifests for surface members."""

    manifests: list[dict[str, Any]] = []
    for member in members:
        key = Path(member["source_filename"]).stem.lower()
        import re

        key = re.sub(r"[^a-z0-9]+", "_", key).strip("_")
        path = output_root / "recordings" / key / "manifest.json"
        if path.is_file():
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["_manifest_path"] = path.as_posix()
            manifest["_manifest_sha256"] = sha256_file(path)
            manifests.append(manifest)
    return manifests


def _stage_runtime_summary(manifests: Sequence[Mapping[str, Any]], output_root: Path) -> dict[str, Any]:
    """Aggregate completed-stage runtimes from terminal ledgers."""

    values: dict[str, list[float]] = defaultdict(list)
    for manifest in manifests:
        result_dir = Path(manifest["_manifest_path"]).parent
        ledger_path = result_dir / "stage_ledger.json"
        if not ledger_path.is_file():
            continue
        ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
        for stage in ledger["stages"]:
            if stage["status"] == "complete" and "elapsed_seconds" in stage:
                values[stage["stage"]].append(float(stage["elapsed_seconds"]))
    return {
        stage: {
            "recording_count": len(times),
            "total_seconds": sum(times),
            "median_seconds": statistics.median(times),
            "minimum_seconds": min(times),
            "maximum_seconds": max(times),
        }
        for stage, times in values.items()
    }


def _artifact_bytes(manifest: Mapping[str, Any], kinds: set[str]) -> int:
    """Sum recorded artifact bytes for selected artifact kinds."""

    return sum(
        int(row["size_bytes"])
        for row in manifest.get("artifacts", [])
        if row.get("kind") in kinds
    )


def aggregate_continuous_run(
    *,
    repo_root: Path,
    output_root: Path,
    surface: Mapping[str, Any],
    invocation_results: Sequence[Mapping[str, Any]],
    run_provenance_data: Mapping[str, Any],
    not_attempted: Sequence[str],
    run_started_at: str,
    run_elapsed_seconds: float,
    run_state: str,
    preflight: Mapping[str, Any] | None,
    source_inventory_before: Mapping[str, Any] | None,
    source_inventory_after: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Build run-level status, scientific-stage, runtime, and storage evidence."""

    manifests = _load_terminal_manifests(output_root, surface["members"])
    invocation_counts = Counter(result["status"] for result in invocation_results)
    terminal_counts = Counter(manifest["status"] for manifest in manifests)
    terminal_by_file = {
        manifest["source"]["source_filename"]: manifest["status"] for manifest in manifests
    }
    missing = [
        member["source_filename"]
        for member in surface["members"]
        if member["source_filename"] not in terminal_by_file
    ]

    recording_rows: list[dict[str, Any]] = []
    for manifest in manifests:
        union = manifest.get("detector", {}).get("accepted_union") or {}
        reference = manifest.get("reference", {}).get("estimation") or {}
        interpolation = manifest.get("interpolation") or {}
        ica = manifest.get("ica", {})
        rank = ica.get("rank") or {}
        proposal = ica.get("proposal") or {}
        application = ica.get("application") or {}
        stop = manifest.get("stop_or_failure") or {}
        recording_rows.append(
            {
                "source_filename": manifest["source"]["source_filename"],
                "status": manifest["status"],
                "completion_class": manifest.get(
                    "completion_class", manifest["status"]
                ),
                "qc_warning_codes": [
                    warning["code"] for warning in manifest.get("qc_warnings", [])
                ],
                "stop_or_failure_code": stop.get("code"),
                "stop_or_failure_stage": stop.get("stage"),
                "stop_or_failure_detail": stop.get("detail"),
                "input_bytes": manifest["source"].get("source_size_bytes"),
                "output_bytes": manifest.get("storage", {}).get("output_bytes_excluding_manifest"),
                "total_elapsed_seconds": manifest.get("runtime", {}).get("total_elapsed_seconds"),
                "process_peak_rss_bytes": manifest.get("runtime", {}).get("process_peak_rss_bytes"),
                "accepted_global_bads": union.get("accepted_global_bads", []),
                "accepted_global_bad_count": union.get("accepted_global_bad_count"),
                "psd_report_only_findings": union.get("psd_report_only_findings", []),
                "reference_source_count": reference.get("reference_source_count"),
                "interpolation_count": interpolation.get("interpolation_count"),
                "interpolation_proportion": interpolation.get("interpolation_proportion"),
                "ica_rank": rank.get("estimated_eeg_rank"),
                "ica_proposal_count": proposal.get("proposal_count"),
                "ica_final_exclusions": application.get("final_exclusions", []),
                "post_ica_derivative_written": any(
                    row.get("kind") == "post_ica_continuous"
                    for row in manifest.get("artifacts", [])
                ),
                "pre_ica_derivative_written": any(
                    row.get("kind") == "pre_ica_continuous"
                    for row in manifest.get("artifacts", [])
                ),
                "ica_object_written": any(
                    row.get("kind") == "ica_solution" for row in manifest.get("artifacts", [])
                ),
                "derivatives_reopened_and_valid": all(
                    value.get("valid", False)
                    for value in manifest.get("derivative_reopen_validation", {}).values()
                ),
                "raw_source_unchanged": manifest["source"]
                .get("immutability_validation", {})
                .get("unchanged", False),
                "manifest_path": manifest["_manifest_path"],
                "manifest_sha256": manifest["_manifest_sha256"],
            }
        )

    stage_runtime = _stage_runtime_summary(manifests, output_root)
    full_processing_times = [
        float(row["total_elapsed_seconds"])
        for row in recording_rows
        if row["total_elapsed_seconds"] is not None
    ]
    projection: dict[str, Any] = {"basis_recording_count": len(full_processing_times)}
    if full_processing_times:
        median_seconds = statistics.median(full_processing_times)
        mean_seconds = statistics.mean(full_processing_times)
        projection.update(
            {
                "median_seconds_per_recording": median_seconds,
                "mean_seconds_per_recording": mean_seconds,
                "projected_95_files_seconds_from_median": median_seconds * 95,
                "projected_95_files_hours_from_median": median_seconds * 95 / 3600,
                "projected_95_files_seconds_from_mean": mean_seconds * 95,
                "projected_95_files_hours_from_mean": mean_seconds * 95 / 3600,
            }
        )

    raw_files = sorted((repo_root / "_Data" / "eeg" / "raw").glob("*.edf"))
    total_raw_bytes = sum(path.stat().st_size for path in raw_files)
    surface_input_bytes = sum(
        int(manifest["source"]["source_size_bytes"])
        for manifest in manifests
        if manifest["source"].get("source_size_bytes") is not None
    )
    continuous_bytes = sum(
        _artifact_bytes(manifest, {"pre_ica_continuous", "post_ica_continuous"})
        for manifest in manifests
    )
    ica_sizes = [
        _artifact_bytes(manifest, {"ica_solution"})
        for manifest in manifests
        if _artifact_bytes(manifest, {"ica_solution"}) > 0
    ]
    evidence_sizes = [
        _artifact_bytes(manifest, {"detector_evidence", "ica_component_evidence", "stage_ledger"})
        + Path(manifest["_manifest_path"]).stat().st_size
        for manifest in manifests
    ]
    storage_projection: dict[str, Any] = {
        "raw_edf_file_count": len(raw_files),
        "total_raw_edf_bytes": total_raw_bytes,
        "surface_terminal_input_bytes": surface_input_bytes,
        "surface_terminal_continuous_fif_bytes": continuous_bytes,
        "fif_precision": "single",
        "retention_contract": "one canonical post-ICA continuous FIF plus ICA object for ordinary files; pre-ICA FIF plus unapplied ICA evidence for review stops",
    }
    if surface_input_bytes and continuous_bytes:
        ratio = continuous_bytes / surface_input_bytes
        projected_continuous = total_raw_bytes * ratio
        projected_ica = statistics.median(ica_sizes) * 95 if ica_sizes else 0
        projected_evidence = statistics.median(evidence_sizes) * 95 if evidence_sizes else 0
        canonical_total = projected_continuous + projected_ica + projected_evidence
        storage_projection.update(
            {
                "continuous_fif_to_source_edf_byte_ratio": ratio,
                "projected_95_file_canonical_continuous_bytes": int(projected_continuous),
                "projected_95_file_ica_bytes": int(projected_ica),
                "projected_95_file_manifest_and_evidence_bytes": int(projected_evidence),
                "projected_95_file_canonical_contract_bytes": int(canonical_total),
                "projected_95_file_both_pre_and_post_branches_bytes": int(
                    canonical_total + projected_continuous
                ),
                "projected_increment_for_redundant_second_continuous_branch_bytes": int(
                    projected_continuous
                ),
            }
        )

    peak_values = [
        int(row["process_peak_rss_bytes"])
        for row in recording_rows
        if row["process_peak_rss_bytes"] is not None
    ]
    peak = max(peak_values) if peak_values else None
    operational = {
        "default_recording_concurrency": 1,
        "concurrency_above_one_recommended": False,
        "expected_peak_memory_bytes": peak,
        "overnight_run_appropriate": (
            projection.get("projected_95_files_hours_from_median", float("inf")) <= 16
        ),
        "light_work_likely_while_running": peak is not None and peak <= 10 * 1024**3,
        "recommendation": "Run sequentially on AC power, prevent system sleep, keep ample free disk space, and avoid other memory-intensive work.",
    }

    recordings_dir = output_root / "recordings"
    history_dir = output_root / "history"
    incomplete = sorted(path.name for path in recordings_dir.glob(".*.incomplete")) if recordings_dir.exists() else []
    history = sorted(path.name for path in history_dir.iterdir()) if history_dir.exists() else []
    history_counts = Counter()
    for name in history:
        matched = next(
            (
                state
                for state in ("stale", "incomplete", "invalid", "retry", "force")
                if f"__{state}__" in name
            ),
            "other",
        )
        history_counts[matched] += 1
    inventory_comparison = (
        compare_source_inventories(source_inventory_before, source_inventory_after)
        if source_inventory_before is not None and source_inventory_after is not None
        else None
    )
    explicit_terminal_counts = {
        status: int(terminal_counts.get(status, 0))
        for status in ("complete", "stopped", "failed")
    }
    qc_warning_counts = Counter(
        code for row in recording_rows for code in row["qc_warning_codes"]
    )
    return {
        "schema_version": 1,
        "run_state": run_state,
        "run_started_at": run_started_at,
        "run_elapsed_seconds": run_elapsed_seconds,
        "surface": surface,
        "output_root": output_root.as_posix(),
        "preflight": dict(preflight or {}),
        "run_provenance": run_provenance_data,
        "invocation_results": list(invocation_results),
        "invocation_counts": {
            status: int(invocation_counts.get(status, 0))
            for status in ("complete", "stopped", "failed", "skipped")
        },
        "not_attempted_due_to_bounded_run": list(not_attempted),
        "current_surface_terminal_counts": {
            **explicit_terminal_counts,
            "missing_or_incomplete": len(missing),
        },
        "current_surface_qc_warning_counts": dict(sorted(qc_warning_counts.items())),
        "current_surface_missing_or_incomplete_files": missing,
        "recordings": recording_rows,
        "stage_runtime": stage_runtime,
        "runtime_projection": projection,
        "storage_projection": storage_projection,
        "operational_recommendation": operational,
        "validation": {
            "all_saved_derivatives_reopened_and_valid": bool(recording_rows)
            and all(row["derivatives_reopened_and_valid"] for row in recording_rows),
            "all_terminal_raw_sources_unchanged": bool(recording_rows)
            and all(row["raw_source_unchanged"] for row in recording_rows),
            "all_surface_members_have_terminal_manifests": not missing,
            "full_source_inventory_unchanged": (
                inventory_comparison["unchanged"]
                if inventory_comparison is not None
                else None
            ),
            "id86_has_no_post_ica_derivative": all(
                not row["post_ica_derivative_written"]
                for row in recording_rows
                if row["source_filename"].lower().startswith("demi_86 ")
            ),
            "epochs_written": False,
            "csd_written": False,
            "full_dataset_run_performed": (
                run_state == "finalized"
                and surface.get("surface_kind")
                == "authorized_all_inventoried_readable_edfs_v1"
                and surface.get("surface_size") == 95
                and not not_attempted
            ),
            "saved_signal_recording_count": sum(
                row["pre_ica_derivative_written"] or row["post_ica_derivative_written"]
                for row in recording_rows
            ),
        },
        "cache_state": {
            "incomplete_directories": incomplete,
            "preserved_history_directories": history,
            "preserved_history_counts": dict(sorted(history_counts.items())),
            "stale_incomplete_or_superseded_detection_enabled": True,
        },
        "source_inventory": {
            "before": source_inventory_before,
            "after": source_inventory_after,
            "comparison": inventory_comparison,
        },
    }


def markdown_summary(aggregate: Mapping[str, Any]) -> str:
    """Render a concise local Markdown summary from aggregate run evidence."""

    counts = aggregate["current_surface_terminal_counts"]
    runtime = aggregate["runtime_projection"]
    storage = aggregate["storage_projection"]
    validation = aggregate["validation"]
    surface = aggregate["surface"]
    production = surface["surface_kind"] == "authorized_all_inventoried_readable_edfs_v1"
    lines = [
        "# Continuous preprocessing run",
        "",
        f"Surface: {surface['surface_kind']} ({surface['surface_size']} EDF files).",
        f"Run state: {aggregate['run_state']}.",
        (
            "This is the authorized all-readable-EDF production surface."
            if production
            else "This is the separate saved architecture-validation cohort."
        ),
        "Continuous status does not imply event, epoch, or analytic eligibility.",
        "",
        "## Current terminal state",
        "",
        f"- Complete: {counts.get('complete', 0)}",
        f"- Stopped: {counts.get('stopped', 0)}",
        f"- Failed: {counts.get('failed', 0)}",
        f"- Missing or incomplete: {counts.get('missing_or_incomplete', 0)}",
        f"- QC warnings: {aggregate.get('current_surface_qc_warning_counts', {})}",
        "",
        "## Validation",
        "",
        f"- Saved derivatives reopened and valid: {validation['all_saved_derivatives_reopened_and_valid']}",
        f"- Raw checksums, sizes, and mtimes unchanged: {validation['all_terminal_raw_sources_unchanged']}",
        f"- Full source inventory unchanged: {validation['full_source_inventory_unchanged']}",
        f"- ID 86 has no post-ICA derivative: {validation['id86_has_no_post_ica_derivative']}",
        "- Epochs, AutoReject, and CSD were not run.",
        "",
        "## Projection",
        "",
        f"- Median-derived runtime for 95 files: {runtime.get('projected_95_files_hours_from_median', float('nan')):.2f} h",
        f"- Canonical full-run storage: {storage.get('projected_95_file_canonical_contract_bytes', 0) / 1024**3:.2f} GiB",
        f"- Both continuous branches for every file: {storage.get('projected_95_file_both_pre_and_post_branches_bytes', 0) / 1024**3:.2f} GiB",
        "- Recommended concurrency: one recording at a time.",
        "",
        "## Per recording",
        "",
        "| File | Status | Completion class | Global bads | Interpolated | Reference sources | ICA rank | Proposals |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in aggregate["recordings"]:
        lines.append(
            f"| {row['source_filename']} | {row['status']} | {row['completion_class']} | {row['accepted_global_bad_count']} | "
            f"{row['interpolation_count']} | {row['reference_source_count']} | "
            f"{row['ica_rank']} | {row['ica_proposal_count']} |"
        )
    return "\n".join(lines) + "\n"


def reopen_current_surface_derivatives(
    *,
    output_root: Path,
    raw_root: Path,
    surface: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Reopen and independently revalidate every current Raw/ICA FIF artifact.

    Args:
        output_root: Authorized validation or production derivative root.
        raw_root: Immutable source EDF root, used for path-safety checks.
        surface: Complete selected surface package.
        config: Accepted tracked configuration.

    Returns:
        Per-file and aggregate validation evidence. Missing manifests, stray
        current FIF files, or any MNE/hash/metadata mismatch make the result
        invalid.

    Side effects:
        Reopens and scans current saved FIF files read-only. Writes nothing.
    """

    rows: list[dict[str, Any]] = []
    raw_artifact_count = 0
    ica_artifact_count = 0
    for member in surface["members"]:
        filename = str(member["source_filename"])
        key = Path(filename).stem.lower()
        import re

        key = re.sub(r"[^a-z0-9]+", "_", key).strip("_")
        result_dir = output_root / "recordings" / key
        manifest_path = result_dir / "manifest.json"
        if not manifest_path.is_file():
            rows.append(
                {
                    "source_filename": filename,
                    "valid": False,
                    "errors": ["manifest_missing"],
                    "artifacts": [],
                }
            )
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        declared_fifs = {
            str(artifact["relative_path"])
            for artifact in manifest.get("artifacts", [])
            if str(artifact.get("relative_path", "")).lower().endswith(".fif")
        }
        observed_fifs = {
            path.relative_to(result_dir).as_posix()
            for path in result_dir.rglob("*.fif")
            if path.is_file()
        }
        file_errors: list[str] = []
        if declared_fifs != observed_fifs:
            file_errors.append("declared_and_observed_fif_inventory_mismatch")
        artifact_rows: list[dict[str, Any]] = []
        for artifact in manifest.get("artifacts", []):
            kind = artifact.get("kind")
            if kind not in {
                "pre_ica_continuous",
                "post_ica_continuous",
                "ica_solution",
            }:
                continue
            path = result_dir / artifact["relative_path"]
            try:
                validate_derivative_path(
                    path,
                    output_root,
                    raw_root,
                )
                if kind in {"pre_ica_continuous", "post_ica_continuous"}:
                    raw_artifact_count += 1
                    evidence = validate_saved_raw(
                        path,
                        expected_sample_count=int(manifest["source"]["sample_count"]),
                        expected_sfreq=float(manifest["source"]["sampling_frequency_hz"]),
                        expected_highpass=float(
                            config["analysis_filter"]["high_pass_hz"]
                        ),
                        expected_lowpass=float(config["analysis_filter"]["low_pass_hz"]),
                        expected_sha256=str(artifact["sha256"]),
                        chunk_seconds=float(
                            config["validation"]["finite_scan_chunk_seconds"]
                        ),
                    )
                else:
                    ica_artifact_count += 1
                    application = manifest.get("ica", {}).get("application") or {}
                    evidence = validate_saved_ica(
                        path,
                        expected_components=int(
                            manifest["ica"]["rank"]["estimated_eeg_rank"]
                        ),
                        expected_exclusions=list(
                            application.get("final_exclusions", [])
                        ),
                        expected_sha256=str(artifact["sha256"]),
                    )
            except Exception as error:  # Keep validating later artifacts/files.
                evidence = {
                    "valid": False,
                    "errors": [f"{type(error).__name__}: {error}"],
                    "reopened_with_mne": False,
                }
            artifact_rows.append(
                {
                    "kind": kind,
                    "relative_path": artifact["relative_path"],
                    "validation": evidence,
                }
            )
            if not evidence["valid"]:
                file_errors.append(f"invalid:{artifact['relative_path']}")
        if filename.lower().startswith("demi_86 ") and "post_ica_raw.fif" in observed_fifs:
            file_errors.append("id86_unauthorized_post_ica_derivative")
        rows.append(
            {
                "source_filename": filename,
                "terminal_status": manifest.get("status"),
                "completion_class": manifest.get(
                    "completion_class", manifest.get("status")
                ),
                "qc_warning_codes": [
                    warning["code"] for warning in manifest.get("qc_warnings", [])
                ],
                "valid": not file_errors,
                "errors": file_errors,
                "artifacts": artifact_rows,
            }
        )
    return {
        "schema_version": 1,
        "surface_kind": surface["surface_kind"],
        "surface_size": surface["surface_size"],
        "recording_manifest_count": sum("terminal_status" in row for row in rows),
        "continuous_raw_artifact_count": raw_artifact_count,
        "ica_artifact_count": ica_artifact_count,
        "all_current_fif_and_ica_artifacts_reopened_and_valid": bool(rows)
        and all(row["valid"] for row in rows),
        "recordings": rows,
    }


def execute_members(
    members: Sequence[Mapping[str, Any]],
    processor: Any,
    progress_callback: Any | None = None,
) -> list[dict[str, Any]]:
    """Execute recording jobs sequentially and continue after one failure.

    Args:
        members: Ordered cohort member mappings.
        processor: Callable accepting one member and returning a result mapping.
        progress_callback: Optional callable receiving ``(index, results)``
            after every attempt.

    Returns:
        One result per member in input order. Raised per-recording exceptions
        become failed result rows so the next member still runs.

    Side effects:
        Determined by the supplied callables. This helper itself writes nothing.
    """

    results: list[dict[str, Any]] = []
    for index, member in enumerate(members, start=1):
        try:
            result = dict(processor(member))
        except Exception as error:
            source_path = Path(member["source_path"])
            result = {
                "source_filename": member["source_filename"],
                "status": "failed",
                "reason": "unpublished_incomplete_result",
                "exception_type": type(error).__name__,
                "exception_message": str(error),
                "total_elapsed_seconds": None,
                "input_bytes": source_path.stat().st_size if source_path.is_file() else None,
                "output_bytes": None,
            }
        results.append(result)
        if progress_callback is not None:
            progress_callback(index, results)
    return results


def compare_completed_ica_routing(
    output_root: Path, members: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """Run the read-only historical-rule regression over original completes."""

    rows: list[dict[str, Any]] = []
    for manifest in _load_terminal_manifests(output_root, members):
        if manifest.get("status") != "complete" or manifest.get("repair"):
            continue
        rows.append(
            compare_completed_manifest_to_historical_rule(
                Path(manifest["_manifest_path"]).parent
            )
        )
    return {
        "completed_recordings_checked": len(rows),
        "set_decision_difference_count": sum(
            not row["same_selected_set"] for row in rows
        ),
        "ordering_only_difference_count": sum(
            row["same_selected_set"] and not row["same_recorded_order"] for row in rows
        ),
        "set_decision_differences": [
            row for row in rows if not row["same_selected_set"]
        ],
        "ordering_only_differences": [
            row for row in rows if row["same_selected_set"] and not row["same_recorded_order"]
        ],
        "rows": rows,
    }


def run_historical_ica_routing_repair(
    *, repo_root: Path, config_path: Path
) -> dict[str, Any]:
    """Repair exactly the superseded over-two ICA stops, sequentially.

    The target set is discovered from current stopped manifests on the first
    invocation and from the explicit repair marker on unchanged reruns. All 95
    raw EDFs are hashed before and after, but no EDF is opened through MNE.
    Progress and current aggregate state are flushed after every recording.
    """

    config = load_config(config_path)
    raw_root = (repo_root / config["paths"]["raw_root"]).resolve()
    output_root, _ = output_root_for_mode(repo_root, config, "production")
    preflight = validate_production_preflight(repo_root, config, output_root)
    surface = build_production_surface(repo_root, config)
    manifests = _load_terminal_manifests(output_root, surface["members"])
    by_name = {
        manifest["source"]["source_filename"]: manifest for manifest in manifests
    }
    targets = [
        member
        for member in surface["members"]
        if is_historical_ica_repair_target(by_name.get(member["source_filename"], {}))
    ]
    if len(targets) != 25:
        raise RuntimeError(f"Expected the bounded 25-recording repair surface; found {len(targets)}.")
    regression = compare_completed_ica_routing(output_root, surface["members"])
    if regression["set_decision_difference_count"]:
        raise RuntimeError(
            "Historical ICA rule changes completed exclusion sets; repair aborted before writes."
        )

    before = snapshot_source_inventory(surface["members"])
    execution_context = {
        "surface_mode": "production_historical_ica_artifact_repair",
        "repair_kind": REPAIR_KIND,
        "target_count": 25,
        "recording_concurrency": 1,
        "earlier_preprocessing_rerun": False,
        "ica_refit": False,
    }
    provenance = run_provenance(
        repo_root,
        config,
        production_code_paths(repo_root),
        execution_context=execution_context,
    )
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ") + "__" + provenance["code"]["sha256"][:10]
    run_dir = output_root / "repair_runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    atomic_write_json(run_dir / "regression.json", regression, output_root, raw_root)
    atomic_write_json(run_dir / "source_inventory_before.json", before, output_root, raw_root)

    started = time.perf_counter()

    def processor(member: Mapping[str, Any]) -> dict[str, Any]:
        return repair_recording_from_saved_ica(
            Path(member["source_path"]),
            output_root=output_root,
            raw_root=raw_root,
            config=config,
            run_provenance=provenance,
        )

    def flush(index: int, results: list[dict[str, Any]]) -> None:
        current = aggregate_continuous_run(
            repo_root=repo_root,
            output_root=output_root,
            surface=surface,
            invocation_results=results,
            run_provenance_data=provenance,
            not_attempted=[row["source_filename"] for row in targets[index:]],
            run_started_at=run_id,
            run_elapsed_seconds=time.perf_counter() - started,
            run_state="running",
            preflight=preflight,
            source_inventory_before=before,
            source_inventory_after=None,
        )
        atomic_write_json(
            run_dir / "repair_progress.json",
            {"attempted": index, "target_count": 25, "results": results},
            output_root,
            raw_root,
        )
        atomic_write_json(run_dir / "run_manifest.json", current, output_root, raw_root)
        atomic_write_text(run_dir / "run_summary.md", markdown_summary(current), output_root, raw_root)

    results = execute_members(targets, processor, flush)
    after = snapshot_source_inventory(surface["members"])
    immutability = compare_source_inventories(before, after)
    atomic_write_json(run_dir / "source_inventory_after.json", after, output_root, raw_root)
    atomic_write_json(run_dir / "source_inventory_comparison.json", immutability, output_root, raw_root)
    aggregate = aggregate_continuous_run(
        repo_root=repo_root,
        output_root=output_root,
        surface=surface,
        invocation_results=results,
        run_provenance_data=provenance,
        not_attempted=[],
        run_started_at=run_id,
        run_elapsed_seconds=time.perf_counter() - started,
        run_state="finalized",
        preflight=preflight,
        source_inventory_before=before,
        source_inventory_after=after,
    )
    aggregate.update(
        {
            "repair_kind": REPAIR_KIND,
            "repair_target_count": 25,
            "regression": regression,
            "raw_source_immutability": immutability,
            "run_directory": run_dir.as_posix(),
        }
    )
    atomic_write_json(run_dir / "run_manifest.json", aggregate, output_root, raw_root)
    atomic_write_text(run_dir / "run_summary.md", markdown_summary(aggregate), output_root, raw_root)
    if not immutability["unchanged"]:
        raise RuntimeError("Raw EDF source inventory changed during ICA repair.")
    return aggregate


def run_continuous_surface(
    *,
    repo_root: Path,
    config_path: Path,
    surface_mode: str = "validation",
    max_files: int | None = None,
    only_recording: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Execute a bounded or complete validation/production invocation.

    Args:
        repo_root: Repository root.
        config_path: Tracked configuration path.
        surface_mode: ``validation`` or explicitly authorized ``production``.
        max_files: Optional bounded prefix for smoke/interruption validation.
        only_recording: Optional exact surface filename for focused processing.
        force: Explicitly preserve and recompute existing selected results.

    Returns:
        Aggregate run manifest plus ``run_directory``.

    Side effects:
        Hashes the full selected raw surface before and after the invocation.
        Writes versioned run evidence and per-recording derivatives only under
        the mode's authorized ignored root.
    """

    if surface_mode not in SURFACE_MODES:
        raise ValueError(f"surface_mode must be one of {SURFACE_MODES}; found {surface_mode!r}")
    if force and only_recording is None:
        raise ValueError("Force recomputation requires one exact --recording filename.")
    if max_files is not None and only_recording is not None:
        raise ValueError("Use either --max-files or --recording, not both.")
    config = load_config(config_path)
    raw_root = (repo_root / config["paths"]["raw_root"]).resolve()
    output_root, configured_output_root = output_root_for_mode(
        repo_root, config, surface_mode
    )
    preflight = (
        validate_production_preflight(repo_root, config, output_root)
        if surface_mode in {"production", "production_v2"}
        else {
            "validation_output_root": output_root.as_posix(),
            "recording_concurrency": 1,
        }
    )
    surface = build_surface(repo_root, config, surface_mode)
    members = list(surface["members"])
    if only_recording is not None:
        members = [member for member in members if member["source_filename"] == only_recording]
        if len(members) != 1:
            raise ValueError(
                f"--recording must exactly match one {surface_mode} member: {only_recording!r}"
            )
    if max_files is not None:
        if max_files < 1:
            raise ValueError("--max-files must be positive.")
        members = members[:max_files]

    run_started_at = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    source_inventory_before = snapshot_source_inventory(surface["members"])
    execution_context = {
        "surface_mode": surface_mode,
        "surface_kind": surface["surface_kind"],
        "surface_size": surface["surface_size"],
        "surface_sha256": sha256_bytes(canonical_json_bytes(surface)),
        "configured_output_root": configured_output_root,
        "resolved_output_root": output_root.as_posix(),
        "line_noise_enabled": bool(config["line_noise"]["enabled"]),
        "recording_concurrency": 1,
    }
    run_provenance_data = run_provenance(
        repo_root,
        config,
        production_code_paths(repo_root),
        execution_context=execution_context,
    )
    run_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        + "__"
        + run_provenance_data["code"]["sha256"][:10]
    )
    run_dir = output_root / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    atomic_write_json(run_dir / "surface_selection.json", surface, output_root, raw_root)
    atomic_write_json(run_dir / "preflight.json", preflight, output_root, raw_root)
    atomic_write_json(
        run_dir / "source_inventory_before.json",
        source_inventory_before,
        output_root,
        raw_root,
    )

    attempted_names = {member["source_filename"] for member in members}
    invocation_position = {
        member["source_filename"]: index
        for index, member in enumerate(members, start=1)
    }
    not_attempted = [
        member["source_filename"]
        for member in surface["members"]
        if member["source_filename"] not in attempted_names
    ]

    def processor(member: Mapping[str, Any]) -> dict[str, Any]:
        """Process one selected member with the shared fixed run state."""

        index = invocation_position[member["source_filename"]]
        print(
            f"[{index}/{len(members)}] Starting {member['source_filename']}",
            flush=True,
        )

        return process_recording(
            Path(member["source_path"]),
            repo_root=repo_root,
            raw_root=raw_root,
            output_root=output_root,
            config=config,
            run_provenance=run_provenance_data,
            selection_context=member,
            force=force,
        )

    def write_aggregate(
        results: Sequence[Mapping[str, Any]], run_state: str,
        source_inventory_after: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Flush a current aggregate manifest and summary atomically."""

        aggregate = aggregate_continuous_run(
            repo_root=repo_root,
            output_root=output_root,
            surface=surface,
            invocation_results=results,
            run_provenance_data=run_provenance_data,
            not_attempted=not_attempted,
            run_started_at=run_started_at,
            run_elapsed_seconds=time.perf_counter() - started,
            run_state=run_state,
            preflight=preflight,
            source_inventory_before=source_inventory_before,
            source_inventory_after=source_inventory_after,
        )
        aggregate["run_id"] = run_id
        aggregate["run_directory"] = run_dir.as_posix()
        atomic_write_json(run_dir / "run_manifest.json", aggregate, output_root, raw_root)
        atomic_write_text(
            run_dir / "run_summary.md", markdown_summary(aggregate), output_root, raw_root
        )
        return aggregate

    def flush_progress(index: int, results: list[dict[str, Any]]) -> None:
        """Atomically flush invocation and aggregate state after every recording."""

        latest = results[-1]
        elapsed = latest.get("total_elapsed_seconds")
        elapsed_text = f"{elapsed:.1f}s" if isinstance(elapsed, (int, float)) else "n/a"
        print(
            f"[{index}/{len(members)}] {latest['source_filename']}: "
            f"{latest['status']} ({latest.get('reason', 'completed')}, {elapsed_text})",
            flush=True,
        )

        atomic_write_json(
            run_dir / "progress.json",
            {
                "run_id": run_id,
                "attempted": index,
                "planned_for_invocation": len(members),
                "results": results,
                "progress_flushed_after_each_recording": True,
            },
            output_root,
            raw_root,
        )
        write_aggregate(results, "running")

    write_aggregate([], "running")
    results = execute_members(members, processor, flush_progress)
    source_inventory_after = snapshot_source_inventory(surface["members"])
    comparison = compare_source_inventories(
        source_inventory_before, source_inventory_after
    )
    atomic_write_json(
        run_dir / "source_inventory_after.json",
        source_inventory_after,
        output_root,
        raw_root,
    )
    atomic_write_json(
        run_dir / "source_inventory_comparison.json", comparison, output_root, raw_root
    )
    aggregate = write_aggregate(results, "finalized", source_inventory_after)
    if not comparison["unchanged"]:
        raise RuntimeError("Raw EDF source inventory changed during the invocation.")
    return aggregate


def run_validation_cohort(
    *,
    repo_root: Path,
    config_path: Path,
    max_files: int | None = None,
    only_recording: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Backward-compatible wrapper for the separate validation surface."""

    return run_continuous_surface(
        repo_root=repo_root,
        config_path=config_path,
        surface_mode="validation",
        max_files=max_files,
        only_recording=only_recording,
        force=force,
    )


def verify_current_surface(
    *, repo_root: Path, config_path: Path, surface_mode: str
) -> dict[str, Any]:
    """Reopen all current surface FIF artifacts and save verification evidence."""

    config = load_config(config_path)
    raw_root = (repo_root / config["paths"]["raw_root"]).resolve()
    output_root, _ = output_root_for_mode(repo_root, config, surface_mode)
    surface = build_surface(repo_root, config, surface_mode)
    evidence = reopen_current_surface_derivatives(
        output_root=output_root, raw_root=raw_root, surface=surface, config=config
    )
    status_counts = Counter(
        row["terminal_status"]
        for row in evidence["recordings"]
        if "terminal_status" in row
    )
    evidence["current_surface_terminal_counts"] = {
        status: int(status_counts.get(status, 0))
        for status in ("complete", "stopped", "failed")
    }
    evidence["current_surface_terminal_counts"]["missing_or_incomplete"] = (
        int(surface["surface_size"]) - int(evidence["recording_manifest_count"])
    )
    evidence["current_surface_qc_warning_counts"] = dict(
        sorted(
            Counter(
                code
                for row in evidence["recordings"]
                for code in row.get("qc_warning_codes", [])
            ).items()
        )
    )
    verification_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    verification_path = output_root / "verifications" / f"{verification_id}.json"
    atomic_write_json(verification_path, evidence, output_root, raw_root)
    evidence["verification_path"] = verification_path.as_posix()
    return evidence
