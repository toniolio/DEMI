"""Sequential validation-cohort execution and aggregate run evidence.

The runner selects the fixed technical validation cohort, hashes code/config/
environment state, processes one EDF at a time, flushes progress after every
recording, and writes a versioned run directory with aggregate runtime,
storage, detector, reference, interpolation, ICA, derivative-validation, and
raw-immutability summaries.

This runner has no all-recording mode. ``--max-files`` supports a bounded first
pass that can be resumed later; valid completed/stopped recording directories
are then hash-validated and skipped. It writes only below the authorized
ignored validation root and never constructs epochs, runs AutoReject, computes
CSD, or changes event/participant eligibility.
"""

from __future__ import annotations

import json
import statistics
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from .cohort import select_validation_cohort
from .contracts import load_config, sha256_file
from .pipeline import process_recording
from .storage import (
    atomic_write_json,
    atomic_write_text,
    run_provenance,
    validate_output_root,
)


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


def _load_terminal_manifests(output_root: Path, members: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Read currently published per-recording manifests for cohort members."""

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


def aggregate_validation_run(
    *,
    repo_root: Path,
    output_root: Path,
    cohort: Mapping[str, Any],
    invocation_results: Sequence[Mapping[str, Any]],
    run_provenance_data: Mapping[str, Any],
    not_attempted: Sequence[str],
    run_started_at: str,
    run_elapsed_seconds: float,
) -> dict[str, Any]:
    """Build run-level status, scientific-stage, runtime, and storage evidence."""

    manifests = _load_terminal_manifests(output_root, cohort["members"])
    invocation_counts = Counter(result["status"] for result in invocation_results)
    terminal_counts = Counter(manifest["status"] for manifest in manifests)
    terminal_by_file = {
        manifest["source"]["source_filename"]: manifest["status"] for manifest in manifests
    }
    missing = [
        member["source_filename"]
        for member in cohort["members"]
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
        recording_rows.append(
            {
                "source_filename": manifest["source"]["source_filename"],
                "status": manifest["status"],
                "stop_or_failure_code": (manifest.get("stop_or_failure") or {}).get("code"),
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
    cohort_input_bytes = sum(
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
        "cohort_terminal_input_bytes": cohort_input_bytes,
        "cohort_terminal_continuous_fif_bytes": continuous_bytes,
        "fif_precision": "single",
        "retention_contract": "one canonical post-ICA continuous FIF plus ICA object for ordinary files; pre-ICA FIF plus unapplied ICA evidence for review stops",
    }
    if cohort_input_bytes and continuous_bytes:
        ratio = continuous_bytes / cohort_input_bytes
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
    return {
        "schema_version": 1,
        "run_started_at": run_started_at,
        "run_elapsed_seconds": run_elapsed_seconds,
        "cohort": cohort,
        "run_provenance": run_provenance_data,
        "invocation_results": list(invocation_results),
        "invocation_counts": dict(sorted(invocation_counts.items())),
        "not_attempted_due_to_bounded_run": list(not_attempted),
        "current_cohort_terminal_counts": {
            **dict(sorted(terminal_counts.items())),
            "missing_or_incomplete": len(missing),
        },
        "current_cohort_missing_or_incomplete_files": missing,
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
            "id86_has_no_post_ica_derivative": all(
                not row["post_ica_derivative_written"]
                for row in recording_rows
                if row["source_filename"].lower().startswith("demi_86 ")
            ),
            "epochs_written": False,
            "csd_written": False,
            "full_dataset_run_performed": False,
            "saved_signal_recording_count": sum(
                row["pre_ica_derivative_written"] or row["post_ica_derivative_written"]
                for row in recording_rows
            ),
        },
        "cache_state": {
            "incomplete_directories": incomplete,
            "preserved_history_directories": history,
            "stale_incomplete_or_superseded_detection_enabled": True,
        },
    }


def markdown_summary(aggregate: Mapping[str, Any]) -> str:
    """Render a concise local Markdown summary from aggregate run evidence."""

    counts = aggregate["current_cohort_terminal_counts"]
    runtime = aggregate["runtime_projection"]
    storage = aggregate["storage_projection"]
    validation = aggregate["validation"]
    lines = [
        "# Continuous preprocessing validation run",
        "",
        f"Cohort size: {aggregate['cohort']['cohort_size']}. This is a saved validation cohort, not the 95-file production run.",
        "",
        "## Current terminal state",
        "",
        f"- Complete: {counts.get('complete', 0)}",
        f"- Stopped: {counts.get('stopped', 0)}",
        f"- Failed: {counts.get('failed', 0)}",
        f"- Missing or incomplete: {counts.get('missing_or_incomplete', 0)}",
        "",
        "## Validation",
        "",
        f"- Saved derivatives reopened and valid: {validation['all_saved_derivatives_reopened_and_valid']}",
        f"- Raw checksums, sizes, and mtimes unchanged: {validation['all_terminal_raw_sources_unchanged']}",
        f"- ID 86 has no post-ICA derivative: {validation['id86_has_no_post_ica_derivative']}",
        "- Epochs, AutoReject, and CSD were not run.",
        "",
        "## Projection",
        "",
        f"- Median projected runtime for 95 files: {runtime.get('projected_95_files_hours_from_median', float('nan')):.2f} h",
        f"- Canonical full-run storage: {storage.get('projected_95_file_canonical_contract_bytes', 0) / 1024**3:.2f} GiB",
        f"- Both continuous branches for every file: {storage.get('projected_95_file_both_pre_and_post_branches_bytes', 0) / 1024**3:.2f} GiB",
        "- Recommended concurrency: one recording at a time.",
        "",
        "## Per recording",
        "",
        "| File | Status | Global bads | Interpolated | Reference sources | ICA rank | Proposals |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in aggregate["recordings"]:
        lines.append(
            f"| {row['source_filename']} | {row['status']} | {row['accepted_global_bad_count']} | "
            f"{row['interpolation_count']} | {row['reference_source_count']} | "
            f"{row['ica_rank']} | {row['ica_proposal_count']} |"
        )
    return "\n".join(lines) + "\n"


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


def run_validation_cohort(
    *,
    repo_root: Path,
    config_path: Path,
    max_files: int | None = None,
    only_recording: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Execute a bounded or complete validation-cohort invocation.

    Args:
        repo_root: Repository root.
        config_path: Tracked configuration path.
        max_files: Optional bounded prefix for interruption/resume validation.
        only_recording: Optional exact cohort filename for focused recomputation.
        force: Explicitly preserve and recompute existing selected results.

    Returns:
        Aggregate run manifest plus ``run_directory``.

    Side effects:
        Writes versioned run evidence and per-recording derivatives under the
        authorized ignored validation root.
    """

    config = load_config(config_path)
    raw_root = (repo_root / config["paths"]["raw_root"]).resolve()
    output_root = (repo_root / config["paths"]["output_root"]).resolve()
    validate_output_root(repo_root, output_root, config["paths"]["output_root"])
    cohort = build_validation_cohort(repo_root, config)
    members = list(cohort["members"])
    if only_recording is not None:
        members = [member for member in members if member["source_filename"] == only_recording]
        if len(members) != 1:
            raise ValueError(f"--recording must exactly match one validation member: {only_recording!r}")
    if max_files is not None:
        if max_files < 1:
            raise ValueError("--max-files must be positive.")
        members = members[:max_files]

    run_provenance_data = run_provenance(
        repo_root, config, production_code_paths(repo_root)
    )
    run_started_at = datetime.now(timezone.utc).isoformat()
    run_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        + "__"
        + run_provenance_data["code"]["sha256"][:10]
    )
    run_dir = output_root / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    atomic_write_json(run_dir / "cohort_selection.json", cohort, output_root, raw_root)

    started = time.perf_counter()
    attempted_names = {member["source_filename"] for member in members}
    not_attempted = [
        member["source_filename"]
        for member in cohort["members"]
        if member["source_filename"] not in attempted_names
    ]
    def processor(member: Mapping[str, Any]) -> dict[str, Any]:
        """Process one selected member with the shared fixed run state."""

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

    def flush_progress(index: int, results: list[dict[str, Any]]) -> None:
        """Atomically flush the invocation ledger after every recording."""

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

    results = execute_members(members, processor, flush_progress)

    aggregate = aggregate_validation_run(
        repo_root=repo_root,
        output_root=output_root,
        cohort=cohort,
        invocation_results=results,
        run_provenance_data=run_provenance_data,
        not_attempted=not_attempted,
        run_started_at=run_started_at,
        run_elapsed_seconds=time.perf_counter() - started,
    )
    aggregate["run_id"] = run_id
    aggregate["run_directory"] = run_dir.as_posix()
    atomic_write_json(run_dir / "run_manifest.json", aggregate, output_root, raw_root)
    atomic_write_text(
        run_dir / "run_summary.md", markdown_summary(aggregate), output_root, raw_root
    )
    return aggregate
