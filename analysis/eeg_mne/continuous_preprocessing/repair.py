"""Bounded artifact-only repair for the superseded ICA proposal guardrail.

This module repairs only production results stopped by the removed
``more_than_two_automatic_eog_component_proposals`` rule. It validates and
reuses the existing pre-ICA FIF, fitted ICA object, rank, scores, manifest, and
hash evidence; applies the recovered historical selection; and atomically
publishes a canonical post-ICA result while preserving the stopped directory
in local history.

It never reads an EDF through MNE, reruns detector/filter/reference/
interpolation stages, refits ICA, creates epochs, runs AutoReject, or computes
CSD. The source EDF is read only for hash/stat immutability checks.
"""

from __future__ import annotations

import gc
import json
import os
import time
from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping

import mne

from .contracts import PIPELINE_VERSION, sha256_file
from .ica import apply_ica_exclusions, decide_component_route, propose_eog_components
from .storage import (
    artifact_descriptor,
    assess_terminal_result,
    atomic_write_json,
    deterministic_manifest_id,
    history_stamp,
    parse_recording_id,
    peak_rss_bytes,
    recording_key,
    recording_provenance,
    utc_now,
    validate_derivative_path,
    validate_saved_ica,
    validate_saved_raw,
    write_ica_atomic,
    write_raw_fif_atomic,
)


SUPERSEDED_STOP_REASON = "more_than_two_automatic_eog_component_proposals"
REPAIR_KIND = "historical_ica_routing_from_saved_artifacts_v1"


def is_historical_ica_repair_target(manifest: Mapping[str, Any]) -> bool:
    """Return whether a current manifest belongs to the bounded repair set."""

    return (
        manifest.get("status") == "stopped"
        and (manifest.get("stop_or_failure") or {}).get("code")
        == SUPERSEDED_STOP_REASON
    ) or (
        manifest.get("status") == "complete"
        and manifest.get("repair", {}).get("repair_kind") == REPAIR_KIND
    )


def _artifact_by_kind(manifest: Mapping[str, Any], kind: str) -> Mapping[str, Any]:
    """Return the unique recorded artifact descriptor for ``kind``."""

    matches = [row for row in manifest.get("artifacts", []) if row.get("kind") == kind]
    if len(matches) != 1:
        raise ValueError(f"Expected one {kind} artifact; found {len(matches)}.")
    return matches[0]


def _descriptor(path: Path, result_dir: Path, kind: str) -> dict[str, Any]:
    """Return the standard artifact descriptor with an explicit kind."""

    return {"kind": kind, **artifact_descriptor(path, result_dir)}


def compare_completed_manifest_to_historical_rule(result_dir: Path) -> dict[str, Any]:
    """Compare one completed manifest's exclusions with the recovered rule.

    This read-only regression check distinguishes exclusion-set changes, which
    would change the signal, from ordering-only changes, which do not.
    """

    manifest = json.loads((result_dir / "manifest.json").read_text(encoding="utf-8"))
    evidence = json.loads((result_dir / "ica_components.json").read_text(encoding="utf-8"))
    proposal = propose_eog_components(evidence["scores"], {})
    historical = list(proposal["proposed_components"])
    current = list(manifest["ica"]["application"]["final_exclusions"])
    return {
        "source_filename": manifest["source"]["source_filename"],
        "current_exclusions": current,
        "historical_exclusions": historical,
        "same_selected_set": set(current) == set(historical),
        "same_recorded_order": current == historical,
        "signal_rewrite_required": set(current) != set(historical),
    }


def repair_recording_from_saved_ica(
    source_path: Path,
    *,
    output_root: Path,
    raw_root: Path,
    config: Mapping[str, Any],
    run_provenance: Mapping[str, Any],
) -> dict[str, Any]:
    """Repair or validly skip one exact superseded-guardrail result.

    Args:
        source_path: Immutable EDF path, used only for hash/stat verification.
        output_root: Authorized production ``continuous_v1`` root.
        raw_root: Immutable EDF directory.
        config: Current validated production configuration.
        run_provenance: Stable repair-run config/code/environment evidence.

    Returns:
        Per-recording repair, skip, or non-target summary.

    Side effects:
        Writes only below ``output_root``. A successful repair preserves the
        prior stopped directory under ``history`` and publishes a replacement
        recording directory atomically.
    """

    started = time.perf_counter()
    key = recording_key(source_path.name)
    final_dir = output_root / "recordings" / key
    work_dir = output_root / "recordings" / f".{key}.ica-repair-incomplete"
    history_dir = output_root / "history"
    for path in (final_dir, work_dir, history_dir):
        validate_derivative_path(path, output_root, raw_root)
    if not final_dir.is_dir():
        return {"source_filename": source_path.name, "status": "incomplete", "reason": "current_result_missing"}

    manifest_path = final_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_hash = sha256_file(source_path)
    expected_provenance = recording_provenance(run_provenance, source_hash)
    if manifest.get("repair", {}).get("repair_kind") == REPAIR_KIND:
        assessment = assess_terminal_result(final_dir, expected_provenance)
        if assessment["cache_status"] == "valid" and manifest.get("status") == "complete":
            return {
                "source_filename": source_path.name,
                "status": "skipped",
                "reason": "valid_repaired_result_unchanged",
                "result_dir": final_dir.as_posix(),
                "manifest_id": manifest["manifest_id"],
                "total_elapsed_seconds": time.perf_counter() - started,
            }
        raise RuntimeError(f"Existing repaired result is not a valid current cache: {assessment['reason']}")

    stop_code = (manifest.get("stop_or_failure") or {}).get("code")
    if manifest.get("status") != "stopped" or stop_code != SUPERSEDED_STOP_REASON:
        return {"source_filename": source_path.name, "status": "not_target", "reason": stop_code or manifest.get("status")}
    if parse_recording_id(source_path.name) == 86:
        raise RuntimeError("ID 86 may not enter the bounded guardrail repair path.")

    old_assessment = assess_terminal_result(
        final_dir, {"fingerprint": manifest["provenance"]["fingerprint"]}
    )
    if old_assessment["cache_status"] != "valid":
        raise RuntimeError(f"Stopped result failed artifact validation: {old_assessment['reason']}")
    source_stat = source_path.stat()
    if (
        source_hash != manifest["source"]["source_sha256"]
        or source_stat.st_size != manifest["source"]["source_size_bytes"]
        or source_stat.st_mtime_ns != manifest["source"]["source_mtime_ns"]
    ):
        raise RuntimeError("Raw EDF hash, size, or mtime differs from the stopped manifest.")

    pre_descriptor = _artifact_by_kind(manifest, "pre_ica_continuous")
    ica_descriptor = _artifact_by_kind(manifest, "ica_solution")
    pre_path = final_dir / pre_descriptor["relative_path"]
    ica_path = final_dir / ica_descriptor["relative_path"]
    rank = manifest["ica"]["rank"]
    pre_validation = validate_saved_raw(
        pre_path,
        expected_sample_count=int(manifest["source"]["sample_count"]),
        expected_sfreq=float(manifest["source"]["sampling_frequency_hz"]),
        expected_highpass=float(config["analysis_filter"]["high_pass_hz"]),
        expected_lowpass=float(config["analysis_filter"]["low_pass_hz"]),
        expected_sha256=pre_descriptor["sha256"],
        chunk_seconds=float(config["validation"]["finite_scan_chunk_seconds"]),
    )
    ica_validation = validate_saved_ica(
        ica_path,
        expected_components=int(rank["estimated_eeg_rank"]),
        expected_exclusions=[],
        expected_sha256=ica_descriptor["sha256"],
    )
    if not pre_validation["valid"] or not ica_validation["valid"]:
        raise RuntimeError("Saved pre-ICA or ICA artifact failed bounded repair validation.")

    component_path = final_dir / _artifact_by_kind(manifest, "ica_component_evidence")["relative_path"]
    component_evidence = json.loads(component_path.read_text(encoding="utf-8"))
    scores = component_evidence["scores"]
    proposal = propose_eog_components(scores, config)
    route = decide_component_route(parse_recording_id(source_path.name), scores, proposal, config)
    if not route["automatic_application_authorized"] or len(route["final_exclusions"]) <= 2:
        raise RuntimeError("Recovered rule did not authorize the expected over-two repair.")

    if work_dir.exists():
        history_dir.mkdir(parents=True, exist_ok=True)
        os.replace(work_dir, history_dir / f"{key}__incomplete_ica_repair__{history_stamp()}")
    work_dir.mkdir(parents=False, exist_ok=False)
    raw = None
    ica = None
    try:
        raw = mne.io.read_raw_fif(pre_path, preload=True, verbose="ERROR")
        ica = mne.preprocessing.read_ica(ica_path, verbose="ERROR")
        application = apply_ica_exclusions(ica, raw, route)

        detector_payload = json.loads(
            (final_dir / _artifact_by_kind(manifest, "detector_evidence")["relative_path"]).read_text(encoding="utf-8")
        )
        detector_out = work_dir / "detector_criteria.json"
        atomic_write_json(detector_out, detector_payload, output_root, raw_root)
        component_evidence.update(
            {"proposal": proposal, "route": route, "application": application}
        )
        component_out = work_dir / "ica_components.json"
        atomic_write_json(component_out, component_evidence, output_root, raw_root)

        raw_out = work_dir / "post_ica_raw.fif"
        raw_write = write_raw_fif_atomic(
            raw,
            raw_out,
            fif_format=config["derivatives"]["fif_format"],
            output_root=output_root,
            raw_root=raw_root,
        )
        raw_reopen = validate_saved_raw(
            raw_out,
            expected_sample_count=int(raw.n_times),
            expected_sfreq=float(raw.info["sfreq"]),
            expected_highpass=float(config["analysis_filter"]["high_pass_hz"]),
            expected_lowpass=float(config["analysis_filter"]["low_pass_hz"]),
            expected_sha256=raw_write["sha256"],
            chunk_seconds=float(config["validation"]["finite_scan_chunk_seconds"]),
        )
        ica_out = work_dir / "solution-ica.fif"
        ica_write = write_ica_atomic(ica, ica_out, output_root=output_root, raw_root=raw_root)
        ica_reopen = validate_saved_ica(
            ica_out,
            expected_components=int(rank["estimated_eeg_rank"]),
            expected_exclusions=route["final_exclusions"],
            expected_sha256=ica_write["sha256"],
        )
        if not raw_reopen["valid"] or not ica_reopen["valid"]:
            raise RuntimeError("New repaired derivative failed reopening validation.")

        old_ledger = json.loads((final_dir / "stage_ledger.json").read_text(encoding="utf-8"))
        old_ledger["overall_status"] = "complete"
        old_ledger["pipeline_version"] = PIPELINE_VERSION
        old_ledger["repair"] = {
            "repair_kind": REPAIR_KIND,
            "earlier_preprocessing_stages_rerun": False,
            "ica_refit": False,
            "reused_stages": [row["stage"] for row in old_ledger["stages"][:17]],
            "corrected_stages": [
                "automatic_component_proposal",
                "exception_or_continuation_decision",
                "derivative_write",
                "manifest_finalization",
            ],
            "superseded_terminal_code": SUPERSEDED_STOP_REASON,
        }
        ledger_out = work_dir / "stage_ledger.json"
        atomic_write_json(ledger_out, old_ledger, output_root, raw_root)

        history_name = f"{key}__superseded_ica_rule__{history_stamp()}"
        history_destination = history_dir / history_name
        repaired = deepcopy(manifest)
        repaired.update(
            {
                "created_at": utc_now(),
                "status": "complete",
                "provenance": expected_provenance,
                "runtime": {
                    "total_elapsed_seconds": time.perf_counter() - started,
                    "process_peak_rss_bytes": peak_rss_bytes(),
                    "repair_only": True,
                },
                "ica": {
                    **repaired["ica"],
                    "proposal": proposal,
                    "route": route,
                    "application": application,
                },
                "derivative_contract": {
                    **repaired["derivative_contract"],
                    "ordinary_signal_derivative": "post_ica_only",
                    "pre_ica_review_derivative_written": False,
                    "post_ica_derivative_written": True,
                },
                "derivative_reopen_validation": {
                    "post_ica_raw.fif": raw_reopen,
                    "solution-ica.fif": ica_reopen,
                },
                "repair": {
                    "repair_kind": REPAIR_KIND,
                    "artifact_only": True,
                    "edf_read_with_mne": False,
                    "pyprep_rerun": False,
                    "filter_reference_interpolation_rerun": False,
                    "ica_refit": False,
                    "validated_inputs": {
                        "pre_ica_raw": pre_validation,
                        "ica_solution": ica_validation,
                        "rank": rank,
                    },
                    "superseded_stop_history": {
                        "terminal_code": SUPERSEDED_STOP_REASON,
                        "manifest_id": manifest["manifest_id"],
                        "manifest_sha256": sha256_file(manifest_path),
                        "history_directory": history_name,
                    },
                },
            }
        )
        repaired.pop("stop_or_failure", None)
        repaired["source"]["immutability_validation"] = {
            "unchanged": True,
            "source_sha256_before": source_hash,
            "source_sha256_after": sha256_file(source_path),
            "source_size_bytes_before": source_stat.st_size,
            "source_size_bytes_after": source_path.stat().st_size,
            "source_mtime_ns_before": source_stat.st_mtime_ns,
            "source_mtime_ns_after": source_path.stat().st_mtime_ns,
        }
        if not repaired["source"]["immutability_validation"]["unchanged"] or any(
            repaired["source"]["immutability_validation"][a]
            != repaired["source"]["immutability_validation"][b]
            for a, b in (
                ("source_sha256_before", "source_sha256_after"),
                ("source_size_bytes_before", "source_size_bytes_after"),
                ("source_mtime_ns_before", "source_mtime_ns_after"),
            )
        ):
            raise RuntimeError("Raw EDF changed during artifact-only repair.")

        artifacts = sorted(
            [
                _descriptor(detector_out, work_dir, "detector_evidence"),
                _descriptor(component_out, work_dir, "ica_component_evidence"),
                _descriptor(raw_out, work_dir, "post_ica_continuous"),
                _descriptor(ica_out, work_dir, "ica_solution"),
                _descriptor(ledger_out, work_dir, "stage_ledger"),
            ],
            key=lambda row: row["relative_path"],
        )
        repaired["artifacts"] = artifacts
        repaired["storage"] = {
            "input_bytes": repaired["source"]["source_size_bytes"],
            "output_bytes_excluding_manifest": sum(row["size_bytes"] for row in artifacts),
            "signal_output_bytes": sum(
                row["size_bytes"] for row in artifacts if row["kind"] in {"post_ica_continuous", "ica_solution"}
            ),
            "fif_format": config["derivatives"]["fif_format"],
        }
        repaired["manifest_id"] = deterministic_manifest_id(repaired)
        atomic_write_json(work_dir / "manifest.json", repaired, output_root, raw_root)

        history_dir.mkdir(parents=True, exist_ok=True)
        os.replace(final_dir, history_destination)
        try:
            os.replace(work_dir, final_dir)
        except Exception:
            os.replace(history_destination, final_dir)
            raise
        return {
            "source_filename": source_path.name,
            "status": "complete",
            "reason": "repaired_from_saved_pre_ica_and_ica",
            "result_dir": final_dir.as_posix(),
            "manifest_id": repaired["manifest_id"],
            "final_exclusions": list(route["final_exclusions"]),
            "total_elapsed_seconds": time.perf_counter() - started,
            "earlier_preprocessing_rerun": False,
            "ica_refit": False,
        }
    finally:
        del raw, ica
        gc.collect()
