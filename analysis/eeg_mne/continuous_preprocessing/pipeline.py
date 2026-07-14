"""One-recording production continuous-preprocessing orchestration.

The pipeline executes the accepted 20 stages in order, keeps detector,
reference, interpolation, ICA proposal, and application evidence separate, and
publishes one atomic self-contained result directory. Objective scientific or
technical stops are terminal per-recording results; unexpected derivative
storage/provenance failures remain failed or incomplete without affecting the
next recording.

Inputs are one immutable EDF, the tracked configuration, shared run
provenance, and factual continuous-surface context. Outputs are FIF/JSON files
only below the selected authorized ignored root. Raw EDFs are re-hashed after
processing. This module never writes EDF, constructs epochs, runs AutoReject,
computes CSD, changes event eligibility, or decides participant inclusion.
"""

from __future__ import annotations

import contextlib
import gc
import io
import time
from pathlib import Path
from typing import Any, Callable, Mapping, TypeVar

import mne

from .contracts import (
    EEG_TARGET_CHANNELS,
    MASTOID_CHANNELS,
    PIPELINE_VERSION,
    SCALP_SOURCE_CHANNELS,
    sha256_bytes,
    sha256_file,
)
from .ica import (
    apply_ica_exclusions,
    decide_component_route,
    estimate_eeg_rank,
    fit_ica,
    prepare_ica_fitting_copy,
    propose_eog_components,
    score_eog_components,
)
from .source import (
    apply_channel_types,
    apply_montage,
    read_and_validate_source,
    verify_source_unchanged,
)
from .stages import (
    ObjectiveStop,
    accepted_global_bad_union,
    apply_fir_filter,
    apply_line_noise_branch,
    apply_reference,
    calculate_reference_sources,
    interpolate_global_bads,
    prepare_detector_input,
    run_pyprep_detection,
    validate_finite_continuous,
)
from .storage import (
    ResultStore,
    StageLedger,
    artifact_descriptor,
    atomic_write_json,
    deterministic_manifest_id,
    parse_recording_id,
    peak_rss_bytes,
    recording_key,
    recording_provenance,
    utc_now,
    validate_saved_ica,
    validate_saved_raw,
    write_ica_atomic,
    write_raw_fif_atomic,
)


T = TypeVar("T")


@contextlib.contextmanager
def quiet_third_party_output() -> Any:
    """Capture verbose PyPREP/MNE progress while retaining a compact hash."""

    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr), mne.use_log_level(
        "ERROR"
    ):
        yield stdout, stderr


def _run_stage(
    ledger: StageLedger,
    stage: str,
    operation: Callable[[], T],
    evidence: Callable[[T], Mapping[str, Any]] | None = None,
) -> T:
    """Run one ledger stage and flush its evidence on success."""

    ledger.start(stage)
    result = operation()
    ledger.complete(stage, evidence(result) if evidence else result)  # type: ignore[arg-type]
    return result


def _artifact(path: Path, result_dir: Path, kind: str) -> dict[str, Any]:
    """Add a stable artifact kind to the standard hash/size descriptor."""

    return {"kind": kind, **artifact_descriptor(path, result_dir)}


def _manifest_base(
    *,
    source_path: Path,
    source_sha256: str,
    source_evidence: Mapping[str, Any] | None,
    provenance: Mapping[str, Any],
    selection_context: Mapping[str, Any],
    status: str,
    elapsed_seconds: float,
) -> dict[str, Any]:
    """Create common manifest fields for complete, stopped, and failed results."""

    source = dict(source_evidence or {})
    source.setdefault("source_filename", source_path.name)
    source.setdefault("source_path", source_path.as_posix())
    source.setdefault("source_sha256", source_sha256)
    return {
        "schema_version": 1,
        "pipeline_version": PIPELINE_VERSION,
        "created_at": utc_now(),
        "status": status,
        "source": source,
        "provenance": dict(provenance),
        "continuous_preprocessing_surface_context": dict(selection_context),
        "recording_surface_boundary": {
            "raw_edf_readability": "readable" if source_evidence else "not_confirmed",
            "continuous_preprocessing_eligibility": selection_context.get(
                "continuous_preprocessing_eligibility",
                "selected_validation_cohort_recording",
            ),
            "continuous_preprocessing_status": status,
            "event_source_availability": selection_context.get(
                "event_source_availability", "not_available_to_continuous_pipeline"
            ),
            "accepted_first_pass_epoch_eligibility": "not_evaluated_by_continuous_pipeline",
            "accepted_8905_row_candidate_surface_membership": "not_consumed_by_continuous_pipeline",
            "successful_continuous_preprocessing_implies_epoch_eligibility": False,
        },
        "runtime": {
            "total_elapsed_seconds": elapsed_seconds,
            "process_peak_rss_bytes": peak_rss_bytes(),
        },
    }


def _write_evidence_json(
    work_dir: Path,
    output_root: Path,
    raw_root: Path,
    *,
    detection: Mapping[str, Any] | None,
    union: Mapping[str, Any] | None,
    rank: Mapping[str, Any] | None,
    fit: Mapping[str, Any] | None,
    scores: Mapping[str, Any] | None,
    proposal: Mapping[str, Any] | None,
    route: Mapping[str, Any] | None,
    application: Mapping[str, Any] | None,
) -> list[dict[str, Any]]:
    """Write deterministic detector/component evidence JSON when available."""

    artifacts: list[dict[str, Any]] = []
    if detection is not None:
        detector_path = work_dir / "detector_criteria.json"
        atomic_write_json(
            detector_path,
            {"detection": detection, "accepted_union": union},
            output_root,
            raw_root,
        )
        artifacts.append(_artifact(detector_path, work_dir, "detector_evidence"))
    if any(value is not None for value in (rank, fit, scores, proposal, route, application)):
        component_path = work_dir / "ica_components.json"
        atomic_write_json(
            component_path,
            {
                "rank": rank,
                "fit": fit,
                "scores": scores,
                "proposal": proposal,
                "route": route,
                "application": application,
            },
            output_root,
            raw_root,
        )
        artifacts.append(_artifact(component_path, work_dir, "ica_component_evidence"))
    return artifacts


def _write_signal_and_ica(
    *,
    work_dir: Path,
    output_root: Path,
    raw_root: Path,
    analysis_raw: mne.io.BaseRaw,
    ica: mne.preprocessing.ICA | None,
    rank: Mapping[str, Any] | None,
    final_exclusions: list[int],
    save_post_ica: bool,
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Save and reopen the authorized signal branch and optional ICA object."""

    artifacts: list[dict[str, Any]] = []
    validations: dict[str, Any] = {}
    raw_name = "post_ica_raw.fif" if save_post_ica else "pre_ica_raw.fif"
    raw_path = work_dir / raw_name
    raw_write = write_raw_fif_atomic(
        analysis_raw,
        raw_path,
        fif_format=config["derivatives"]["fif_format"],
        output_root=output_root,
        raw_root=raw_root,
    )
    raw_validation = validate_saved_raw(
        raw_path,
        expected_sample_count=int(analysis_raw.n_times),
        expected_sfreq=float(analysis_raw.info["sfreq"]),
        expected_highpass=float(config["analysis_filter"]["high_pass_hz"]),
        expected_lowpass=float(config["analysis_filter"]["low_pass_hz"]),
        expected_sha256=raw_write["sha256"],
        chunk_seconds=float(config["validation"]["finite_scan_chunk_seconds"]),
    )
    if not raw_validation["valid"]:
        raise RuntimeError(f"Saved continuous derivative failed reopening: {raw_validation}")
    artifacts.append(_artifact(raw_path, work_dir, "post_ica_continuous" if save_post_ica else "pre_ica_continuous"))
    validations[raw_name] = raw_validation

    if ica is not None:
        ica.exclude = list(final_exclusions)
        ica_path = work_dir / "solution-ica.fif"
        ica_write = write_ica_atomic(ica, ica_path, output_root=output_root, raw_root=raw_root)
        ica_validation = validate_saved_ica(
            ica_path,
            expected_components=int(rank["estimated_eeg_rank"]),  # type: ignore[index]
            expected_exclusions=final_exclusions,
            expected_sha256=ica_write["sha256"],
        )
        if not ica_validation["valid"]:
            raise RuntimeError(f"Saved ICA object failed reopening: {ica_validation}")
        artifacts.append(_artifact(ica_path, work_dir, "ica_solution"))
        validations["solution-ica.fif"] = ica_validation
    return artifacts, validations


def _publish_terminal_manifest(
    *,
    store: ResultStore,
    ledger: StageLedger,
    manifest: dict[str, Any],
    artifacts: list[dict[str, Any]],
    raw_root: Path,
) -> dict[str, Any]:
    """Finalize hashes, write manifest atomically, publish, and return summary."""

    ledger_descriptor = _artifact(ledger.path, store.work_dir, "stage_ledger")
    artifact_rows = sorted([*artifacts, ledger_descriptor], key=lambda row: row["relative_path"])
    manifest["artifacts"] = artifact_rows
    manifest["storage"] = {
        "input_bytes": manifest["source"].get("source_size_bytes"),
        "output_bytes_excluding_manifest": sum(row["size_bytes"] for row in artifact_rows),
        "signal_output_bytes": sum(
            row["size_bytes"]
            for row in artifact_rows
            if row["kind"] in {"pre_ica_continuous", "post_ica_continuous", "ica_solution"}
        ),
        "fif_format": manifest.get("derivative_contract", {}).get("fif_format"),
    }
    manifest["manifest_id"] = deterministic_manifest_id(manifest)
    atomic_write_json(store.work_dir / "manifest.json", manifest, store.output_root, raw_root)
    final_dir = store.publish()
    final_manifest = final_dir / "manifest.json"
    return {
        "source_filename": manifest["source"]["source_filename"],
        "status": manifest["status"],
        "reason": (manifest.get("stop_or_failure") or {}).get("code", "completed"),
        "result_dir": final_dir.as_posix(),
        "manifest_sha256": sha256_file(final_manifest),
        "manifest_id": manifest["manifest_id"],
        "total_elapsed_seconds": manifest["runtime"]["total_elapsed_seconds"],
        "process_peak_rss_bytes": manifest["runtime"]["process_peak_rss_bytes"],
        "input_bytes": manifest["storage"]["input_bytes"],
        "output_bytes": manifest["storage"]["output_bytes_excluding_manifest"]
        + final_manifest.stat().st_size,
    }


def process_recording(
    source_path: Path,
    *,
    repo_root: Path,
    raw_root: Path,
    output_root: Path,
    config: Mapping[str, Any],
    run_provenance: Mapping[str, Any],
    selection_context: Mapping[str, Any],
    force: bool = False,
) -> dict[str, Any]:
    """Process, validate, atomically publish, or validly skip one EDF.

    Args:
        source_path: Immutable source EDF.
        repo_root: Repository root.
        raw_root: Configured source directory.
        output_root: Authorized validation or production derivative root.
        config: Validated tracked production configuration.
        run_provenance: Shared config/code/environment provenance.
        selection_context: Deterministic factual cohort-selection evidence.
        force: Explicitly preserve and recompute an existing terminal result.

    Returns:
        Per-recording run summary with terminal status, reason, path, runtime,
        memory, and byte counts.

    Side effects:
        Reads one EDF; writes only inside ``output_root``; publishes one result
        directory atomically; may preserve old local results in history.
    """

    del repo_root  # Included in the signature to keep the public boundary explicit.
    started = time.perf_counter()
    initial_source_hash = sha256_file(source_path)
    provenance = recording_provenance(run_provenance, initial_source_hash)
    store = ResultStore(
        output_root=output_root,
        raw_root=raw_root,
        key=recording_key(source_path.name),
        expected_provenance=provenance,
        force=force,
    )
    cache = store.prepare()
    if cache["action"] == "skip":
        manifest = cache["manifest"]
        return {
            "source_filename": source_path.name,
            "status": "skipped",
            "reason": "valid_terminal_result_unchanged",
            "cached_terminal_status": manifest["status"],
            "result_dir": cache["result_dir"].as_posix(),
            "manifest_id": manifest["manifest_id"],
            "total_elapsed_seconds": time.perf_counter() - started,
            "process_peak_rss_bytes": peak_rss_bytes(),
            "input_bytes": source_path.stat().st_size,
            "output_bytes": sum(row["size_bytes"] for row in manifest["artifacts"]),
        }

    ledger = StageLedger(
        store.work_dir / "stage_ledger.json",
        output_root=output_root,
        raw_root=raw_root,
    )
    raw: mne.io.BaseRaw | None = None
    analysis_raw: mne.io.BaseRaw | None = None
    detector_raw: mne.io.BaseRaw | None = None
    ica_raw: mne.io.BaseRaw | None = None
    ica: mne.preprocessing.ICA | None = None
    source_evidence: dict[str, Any] | None = None
    type_evidence: dict[str, Any] | None = None
    montage_evidence: dict[str, Any] | None = None
    line_evidence: dict[str, Any] | None = None
    filter_evidence: dict[str, Any] | None = None
    detector_input_evidence: dict[str, Any] | None = None
    detection: dict[str, Any] | None = None
    union: dict[str, Any] | None = None
    reference: dict[str, Any] | None = None
    reference_application: dict[str, Any] | None = None
    interpolation: dict[str, Any] | None = None
    post_validation: dict[str, Any] | None = None
    ica_preparation: dict[str, Any] | None = None
    rank: dict[str, Any] | None = None
    fit: dict[str, Any] | None = None
    scores: dict[str, Any] | None = None
    proposal: dict[str, Any] | None = None
    route: dict[str, Any] | None = None
    application: dict[str, Any] | None = None
    artifacts: list[dict[str, Any]] = []
    reopen_validation: dict[str, Any] = {}

    try:
        raw, source_evidence = _run_stage(
            ledger,
            "source_validation_and_read",
            lambda: read_and_validate_source(source_path, raw_root),
            lambda result: result[1],
        )
        if source_evidence["source_sha256"] != initial_source_hash:
            raise ObjectiveStop(
                "source_validation_and_read",
                "source_checksum_or_identity_conflict",
                "Source hash changed between cache assessment and EDF read.",
            )
        type_evidence = _run_stage(
            ledger, "channel_typing", lambda: apply_channel_types(raw)  # type: ignore[arg-type]
        )
        montage_evidence = _run_stage(
            ledger, "montage_application", lambda: apply_montage(raw)  # type: ignore[arg-type]
        )

        def prepare_analysis_line_branch() -> dict[str, Any]:
            nonlocal analysis_raw
            analysis_raw = raw.copy()  # type: ignore[union-attr]
            return apply_line_noise_branch(analysis_raw, config)

        line_evidence = _run_stage(
            ledger, "optional_line_noise_removal", prepare_analysis_line_branch
        )
        filter_evidence = _run_stage(
            ledger,
            "analysis_filtering",
            lambda: apply_fir_filter(
                analysis_raw,  # type: ignore[arg-type]
                high_pass_hz=float(config["analysis_filter"]["high_pass_hz"]),
                low_pass_hz=float(config["analysis_filter"]["low_pass_hz"]),
                settings=config["analysis_filter"],
            ),
        )
        detector_raw, detector_input_evidence = _run_stage(
            ledger,
            "detector_input_preparation",
            lambda: prepare_detector_input(raw, config),  # type: ignore[arg-type]
            lambda result: result[1],
        )

        def detect() -> dict[str, Any]:
            with quiet_third_party_output() as captured:
                result = run_pyprep_detection(detector_raw, config)  # type: ignore[arg-type]
            log_text = captured[0].getvalue() + captured[1].getvalue()
            result["captured_third_party_log_bytes"] = len(log_text.encode("utf-8"))
            result["captured_third_party_log_sha256"] = sha256_bytes(log_text.encode("utf-8"))
            return result

        detection = _run_stage(ledger, "pyprep_criterion_detection", detect)
        union = _run_stage(
            ledger, "accepted_global_bad_union", lambda: accepted_global_bad_union(detection, config)
        )
        del detector_raw
        detector_raw = None
        gc.collect()

        reference = _run_stage(
            ledger,
            "reference_estimation",
            lambda: calculate_reference_sources(union["accepted_global_bads"], config),
        )
        reference_application = _run_stage(
            ledger,
            "reference_application",
            lambda: apply_reference(analysis_raw, reference),  # type: ignore[arg-type]
        )
        interpolation = _run_stage(
            ledger,
            "interpolation",
            lambda: interpolate_global_bads(
                analysis_raw, union["accepted_global_bads"], config  # type: ignore[arg-type]
            ),
        )
        post_validation = _run_stage(
            ledger,
            "post_interpolation_validation",
            lambda: validate_finite_continuous(
                analysis_raw,  # type: ignore[arg-type]
                float(config["validation"]["finite_scan_chunk_seconds"]),
            ),
        )
        ica_raw, ica_preparation = _run_stage(
            ledger,
            "ica_fitting_copy_preparation",
            lambda: prepare_ica_fitting_copy(
                raw,  # type: ignore[arg-type]
                global_bads=union["accepted_global_bads"],
                reference=reference,
                config=config,
            ),
            lambda result: result[1],
        )
        del raw
        raw = None
        gc.collect()
        rank = _run_stage(
            ledger,
            "rank_estimation",
            lambda: estimate_eeg_rank(ica_raw, config),  # type: ignore[arg-type]
        )
        ica, fit = _run_stage(
            ledger,
            "ica_fit",
            lambda: fit_ica(ica_raw, rank, config),  # type: ignore[arg-type]
            lambda result: result[1],
        )
        scores = _run_stage(
            ledger, "component_scoring", lambda: score_eog_components(ica, ica_raw, config)  # type: ignore[arg-type]
        )
        proposal = _run_stage(
            ledger,
            "automatic_component_proposal",
            lambda: propose_eog_components(scores, config),
        )

        def decide_and_apply() -> dict[str, Any]:
            nonlocal route, application
            route = decide_component_route(
                parse_recording_id(source_path.name), scores, proposal, config
            )
            if route["automatic_application_authorized"]:
                application = apply_ica_exclusions(ica, analysis_raw, route)  # type: ignore[arg-type]
            else:
                ica.exclude = []  # type: ignore[union-attr]
                application = {
                    "authorized": False,
                    "final_exclusions": [],
                    "exclusion_count": 0,
                    "ica_application_performed": False,
                    "reason": route["reason"],
                }
            return {"route": route, "application": application}

        _run_stage(ledger, "exception_or_continuation_decision", decide_and_apply)
        del ica_raw
        ica_raw = None
        gc.collect()

        def write_derivatives() -> dict[str, Any]:
            nonlocal artifacts, reopen_validation
            artifacts.extend(
                _write_evidence_json(
                    store.work_dir,
                    output_root,
                    raw_root,
                    detection=detection,
                    union=union,
                    rank=rank,
                    fit=fit,
                    scores=scores,
                    proposal=proposal,
                    route=route,
                    application=application,
                )
            )
            signal_artifacts, validations = _write_signal_and_ica(
                work_dir=store.work_dir,
                output_root=output_root,
                raw_root=raw_root,
                analysis_raw=analysis_raw,  # type: ignore[arg-type]
                ica=ica,
                rank=rank,
                final_exclusions=list(route["final_exclusions"]),  # type: ignore[index]
                save_post_ica=bool(route["automatic_application_authorized"]),  # type: ignore[index]
                config=config,
            )
            artifacts.extend(signal_artifacts)
            reopen_validation.update(validations)
            if parse_recording_id(source_path.name) == 86 and (
                store.work_dir / "post_ica_raw.fif"
            ).exists():
                raise RuntimeError("ID 86 produced an unauthorized post-ICA derivative.")
            return {
                "artifact_count_before_ledger_and_manifest": len(artifacts),
                "reopen_validation": validations,
                "post_ica_derivative_written": bool(
                    route["automatic_application_authorized"]  # type: ignore[index]
                ),
                "pre_ica_review_derivative_written": not bool(
                    route["automatic_application_authorized"]  # type: ignore[index]
                ),
            }

        _run_stage(ledger, "derivative_write", write_derivatives)
        source_immutability = verify_source_unchanged(source_path, source_evidence)
        source_evidence["immutability_validation"] = source_immutability
        terminal_status = (
            "complete" if route["automatic_application_authorized"] else "stopped"  # type: ignore[index]
        )
        _run_stage(
            ledger,
            "manifest_finalization",
            lambda: {
                "terminal_status": terminal_status,
                "source_immutability_passed": source_immutability["unchanged"],
                "derivatives_reopened_and_validated": all(
                    validation["valid"] for validation in reopen_validation.values()
                ),
            },
        )
        ledger.finalize(terminal_status)
        elapsed = time.perf_counter() - started
        manifest = _manifest_base(
            source_path=source_path,
            source_sha256=initial_source_hash,
            source_evidence=source_evidence,
            provenance=provenance,
            selection_context=selection_context,
            status=terminal_status,
            elapsed_seconds=elapsed,
        )
        manifest.update(
            {
                "channel_contract": {
                    "typing": type_evidence,
                    "montage": montage_evidence,
                    "detector_sources": list(SCALP_SOURCE_CHANNELS),
                    "reference_sources_before_bad_exclusion": list(SCALP_SOURCE_CHANNELS),
                    "eeg_targets": list(EEG_TARGET_CHANNELS),
                    "m1_m2_source_excluded_target_retained": list(MASTOID_CHANNELS),
                    "m1_m2_primary_interpolation_candidates": False,
                },
                "line_noise": line_evidence,
                "analysis_filter": filter_evidence,
                "detector": {
                    "input": detector_input_evidence,
                    "criteria": detection,
                    "accepted_union": union,
                },
                "reference": {
                    "estimation": reference,
                    "application": reference_application,
                },
                "interpolation": interpolation,
                "post_interpolation_validation": post_validation,
                "ica": {
                    "preparation": ica_preparation,
                    "rank": rank,
                    "fit": fit,
                    "scores": scores,
                    "proposal": proposal,
                    "route": route,
                    "application": application,
                },
                "derivative_contract": {
                    "fif_format": config["derivatives"]["fif_format"],
                    "ordinary_signal_derivative": config["derivatives"][
                        "ordinary_signal_derivative"
                    ],
                    "ica_object_retained": True,
                    "pre_ica_retained_only_for_review_stop": terminal_status == "stopped",
                    "alternate_branch_reconstructible_from_source_code_config_and_ica": True,
                    "component_time_series_written": False,
                },
                "derivative_reopen_validation": reopen_validation,
            }
        )
        if terminal_status == "stopped":
            manifest["stop_or_failure"] = {
                "kind": "objective_stop",
                "code": route["reason"],  # type: ignore[index]
                "stage": "exception_or_continuation_decision",
                "detail": "Valid pre-ICA derivative and ICA evidence retained; no post-ICA derivative authorized.",
            }
        result = _publish_terminal_manifest(
            store=store,
            ledger=ledger,
            manifest=manifest,
            artifacts=artifacts,
            raw_root=raw_root,
        )
        return result

    except Exception as error:
        active_stage = ledger.active_stage
        if active_stage is None:
            # A post-ledger publication failure leaves the incomplete directory
            # intact so it can never be mistaken for a valid terminal result.
            raise
        objective = isinstance(error, ObjectiveStop)
        stage_is_scientific_or_algorithmic = active_stage not in {
            "derivative_write",
            "manifest_finalization",
        }
        status = "stopped" if objective or stage_is_scientific_or_algorithmic else "failed"
        if objective:
            code = error.code
            detail = error.detail
            error_evidence = error.evidence
        else:
            code = f"{active_stage}_exception"
            detail = f"{type(error).__name__}: {error}"
            error_evidence = {"exception_type": type(error).__name__}
        if status == "stopped":
            ledger.stop(active_stage, code, detail, error_evidence)
        else:
            ledger.fail(active_stage, code, detail, error_evidence)

        # Once the valid analysis-filtered/reference/interpolated branch exists,
        # retain it for an ICA-stage review stop. No partial signal is saved for
        # earlier detector/reference/interpolation stops.
        if (
            status == "stopped"
            and analysis_raw is not None
            and active_stage
            in {
                "rank_estimation",
                "ica_fit",
                "component_scoring",
                "automatic_component_proposal",
                "exception_or_continuation_decision",
            }
        ):
            artifacts.extend(
                _write_evidence_json(
                    store.work_dir,
                    output_root,
                    raw_root,
                    detection=detection,
                    union=union,
                    rank=rank,
                    fit=fit,
                    scores=scores,
                    proposal=proposal,
                    route=route,
                    application=application,
                )
            )
            ica_for_review = ica if fit is not None else None
            if ica_for_review is not None:
                ica_for_review.exclude = []
            signal_artifacts, reopen_validation = _write_signal_and_ica(
                work_dir=store.work_dir,
                output_root=output_root,
                raw_root=raw_root,
                analysis_raw=analysis_raw,
                ica=ica_for_review,
                rank=rank,
                final_exclusions=[],
                save_post_ica=False,
                config=config,
            )
            artifacts.extend(signal_artifacts)

        if source_evidence is not None:
            source_evidence["immutability_validation"] = verify_source_unchanged(
                source_path, source_evidence
            )
        ledger.finalize(status)
        elapsed = time.perf_counter() - started
        manifest = _manifest_base(
            source_path=source_path,
            source_sha256=initial_source_hash,
            source_evidence=source_evidence,
            provenance=provenance,
            selection_context=selection_context,
            status=status,
            elapsed_seconds=elapsed,
        )
        manifest.update(
            {
                "stop_or_failure": {
                    "kind": "objective_stop" if status == "stopped" else "failure",
                    "code": code,
                    "stage": active_stage,
                    "detail": detail,
                    "evidence": error_evidence,
                },
                "channel_contract": {
                    "typing": type_evidence,
                    "montage": montage_evidence,
                    "detector_sources": list(SCALP_SOURCE_CHANNELS),
                    "eeg_targets": list(EEG_TARGET_CHANNELS),
                    "m1_m2_source_excluded_target_retained": list(MASTOID_CHANNELS),
                },
                "detector": {"input": detector_input_evidence, "criteria": detection, "accepted_union": union},
                "reference": {"estimation": reference, "application": reference_application},
                "interpolation": interpolation,
                "ica": {
                    "preparation": ica_preparation,
                    "rank": rank,
                    "fit": fit,
                    "scores": scores,
                    "proposal": proposal,
                    "route": route,
                    "application": application,
                },
                "derivative_contract": {
                    "fif_format": config["derivatives"]["fif_format"],
                    "ordinary_signal_derivative": config["derivatives"][
                        "ordinary_signal_derivative"
                    ],
                    "pre_ica_review_derivative_written": any(
                        row["kind"] == "pre_ica_continuous" for row in artifacts
                    ),
                    "post_ica_derivative_written": False,
                    "component_time_series_written": False,
                },
                "derivative_reopen_validation": reopen_validation,
            }
        )
        return _publish_terminal_manifest(
            store=store,
            ledger=ledger,
            manifest=manifest,
            artifacts=artifacts,
            raw_root=raw_root,
        )
    finally:
        del raw, analysis_raw, detector_raw, ica_raw, ica
        gc.collect()
