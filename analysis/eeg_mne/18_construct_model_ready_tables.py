#!/usr/bin/env python3
"""Construct and validate accepted DEMI EEG model-ready tables.

This is the production stage-18 driver. It consumes only the complete accepted
stage-17 behavioural-lineage and persisted ROI-feature Parquets. It selects
five predeclared ROI/window rows per accepted trial, derives the owner-accepted
accuracy-rating and objective-error predictor representations, publishes
deterministic hypothesis/scope views and machine-readable registries, and
validates predictor-only fixed-effect matrix rank.

Inputs:
    ``_Data/eeg/features_v1`` manifest, validation, behavioural lineage, ROI
    features, definitions, and source-immutability evidence; the frozen
    ``bdat2.rds`` hash declared by stage 17; and the tracked
    ``model_table_config_v1.yaml``.

Outputs:
    Ignored ``_Data/eeg/model_tables_v1`` authoritative general Parquet table,
    deterministic views, configuration snapshot, role/estimand/factor/reason
    registries, count and rank summaries, source immutability evidence,
    machine-readable validation, manifest, state, log, and concise Markdown.

This stage explicitly does not open TFR arrays to reconstruct features, change
eligibility, standardize or bin predictors, impute objective error, instantiate
priors, run prior-predictive simulation, compile or fit a model, calculate an
estimand/interval/test, inspect EEG associations, run sensitivities or
influence analyses, fit a GAMM, apply CSD, or interpret scientific results.

Safety:
    Every EEG value is copied directly from one persisted stage-17 ROI row.
    Accepted feature hashes and TFR size/mtime immutability are checked before
    construction and again before publication. All files are built in a
    sibling staging directory; the fully reopened and validated directory is
    atomically renamed into place. A complete unchanged rerun hashes and
    reopens the accepted output without rewriting it.

Run from the repository root:

    tools/run_model_tables_v1.sh

Use ``--verify-current`` for a validation-only unchanged rerun.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from event_epoch_eligibility import stable_frame_hash
from model_table_construction import (
    EXPECTED_GENERAL_ROWS,
    MODEL_TABLE_NAMESPACE,
    MODEL_TABLE_SCHEMA_VERSION,
    archive_existing_output,
    atomic_write_json,
    atomic_write_parquet,
    atomic_write_text,
    build_general_table,
    controlled_reason_registry,
    count_summary,
    deterministic_views,
    estimand_registry,
    factor_coding_registry,
    file_descriptor,
    frame_descriptor,
    load_and_validate_model_table_config,
    make_staging_directory,
    model_matrix_rank_validation,
    model_role_registry,
    publish_staging_directory,
    require_ignored_output_root,
    sha256_file,
    validate_general_and_views,
)


CONFIG_PATH = Path("analysis/eeg_mne/model_table_config_v1.yaml")
FEATURE_ROOT = Path("_Data/eeg/features_v1")
OUTPUT_ROOT = Path("_Data/eeg/model_tables_v1")
BDAT_RDS = Path("_Private/old_rds/bdat2.rds")
RUN_MANIFEST = "model_table_run_manifest.json"
RUN_STATE = "model_table_run_state.json"
VALIDATION_JSON = "validation/model_table_validation.json"
RANK_PARQUET = "validation/model_matrix_rank.parquet"
SOURCE_IMMUTABILITY_JSON = "source_immutability.json"
SUMMARY_JSON = "model_table_summary.json"
SUMMARY_MD = "model_table_summary.md"
RUN_LOG = "logs/production.log"


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def repo_root_from_script() -> Path:
    """Return and validate the repository root containing stage 18."""

    root = Path(__file__).resolve().parents[2]
    if not (root / ".git").exists() or not (root / CONFIG_PATH).is_file():
        raise RuntimeError("could_not_locate_demi_repository_root")
    return root


def parse_args() -> argparse.Namespace:
    """Parse the validation-only unchanged-reuse option."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify-current",
        action="store_true",
        help="Require and fully validate an unchanged complete stage-18 output.",
    )
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    """Read one required JSON object."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"json_authority_not_object:{path}")
    return value


def json_fingerprint(value: Any) -> str:
    """Hash one JSON-serializable authority seed deterministically."""

    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_descriptor(path: Path, descriptor: Mapping[str, Any]) -> None:
    """Require a file to match its recorded size, hash, and optional row count."""

    if not path.is_file() or path.stat().st_size != int(descriptor["size_bytes"]):
        raise RuntimeError(f"file_size_mismatch:{path}")
    if sha256_file(path) != descriptor["sha256"]:
        raise RuntimeError(f"file_hash_mismatch:{path}")
    if descriptor.get("format") == "parquet" and "row_count" in descriptor:
        if pq.read_metadata(path).num_rows != int(descriptor["row_count"]):
            raise RuntimeError(f"parquet_row_count_mismatch:{path}")


def hashed_snapshot(paths: Iterable[Path], *, root: Path) -> list[dict[str, Any]]:
    """Return stable relative-path, size, mtime, and SHA evidence for files."""

    rows = []
    for path in sorted({Path(item) for item in paths}, key=lambda item: item.as_posix()):
        stat = path.stat()
        rows.append(
            {
                "relative_path": path.relative_to(root).as_posix(),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": sha256_file(path),
            }
        )
    return rows


def size_mtime_snapshot_from_saved(
    saved: Sequence[Mapping[str, Any]], *, repo_root: Path
) -> list[dict[str, Any]]:
    """Re-stat an accepted large-array inventory without reading array content."""

    rows = []
    for item in saved:
        path = repo_root / str(item["relative_path"])
        stat = path.stat()
        rows.append(
            {
                "mtime_ns": stat.st_mtime_ns,
                "relative_path": str(item["relative_path"]),
                "size_bytes": stat.st_size,
            }
        )
    return rows


def load_feature_authority(
    repo_root: Path, config: Mapping[str, Any]
) -> dict[str, Any]:
    """Load and fail-closed validate the accepted stage-17 authority.

    Every stage-17 manifest output is hash-validated. Only the persisted
    behavioural-lineage and ROI-feature Parquets are loaded as table inputs.
    The TFR array inventory is checked for size/mtime identity but no TFR array
    is opened.
    """

    feature_root = repo_root / FEATURE_ROOT
    manifest_path = feature_root / "feature_run_manifest.json"
    accepted = config["accepted_source"]
    if sha256_file(manifest_path) != accepted["feature_run_manifest_sha256"]:
        raise RuntimeError("accepted_feature_run_manifest_hash_mismatch")
    manifest = read_json(manifest_path)
    if (
        manifest.get("status") != "complete"
        or manifest.get("feature_namespace") != "features_v1"
        or manifest.get("authority_fingerprint") != accepted["feature_authority_fingerprint"]
    ):
        raise RuntimeError("accepted_feature_manifest_identity_mismatch")

    output_paths: list[Path] = []
    for label, descriptor in manifest["outputs"].items():
        path = feature_root / descriptor["relative_path"]
        validate_descriptor(path, descriptor)
        output_paths.append(path)
    lineage_path = feature_root / manifest["outputs"]["behavioural_lineage"]["relative_path"]
    roi_path = feature_root / manifest["outputs"]["roi_features"]["relative_path"]
    validation_path = feature_root / manifest["outputs"]["validation"]["relative_path"]
    if (
        sha256_file(lineage_path) != accepted["behavioural_lineage_sha256"]
        or sha256_file(roi_path) != accepted["roi_features_sha256"]
        or sha256_file(validation_path) != accepted["feature_validation_sha256"]
    ):
        raise RuntimeError("accepted_feature_primary_output_hash_mismatch")
    validation = read_json(validation_path)
    if validation.get("status") != "pass" or validation.get("source_tfr_unchanged") is not True:
        raise RuntimeError("accepted_feature_validation_not_passed")

    tfr_manifest_path = repo_root / "_Data/eeg/tfr_v1/tfr_run_manifest.json"
    if sha256_file(tfr_manifest_path) != accepted["tfr_run_manifest_sha256"]:
        raise RuntimeError("accepted_tfr_run_manifest_hash_mismatch")
    feature_tfr_immutability = read_json(feature_root / "source_tfr_immutability.json")
    if feature_tfr_immutability.get("status") != "unchanged":
        raise RuntimeError("accepted_feature_tfr_immutability_not_unchanged")
    tfr_current = size_mtime_snapshot_from_saved(
        feature_tfr_immutability["after"], repo_root=repo_root
    )
    if not (
        tfr_current == feature_tfr_immutability["before"]
        == feature_tfr_immutability["after"]
    ):
        raise RuntimeError("accepted_tfr_array_inventory_changed")

    bdat_path = repo_root / BDAT_RDS
    if (
        manifest["frozen_behaviour"]["path"] != BDAT_RDS.as_posix()
        or manifest["frozen_behaviour"]["sha256"] != accepted["frozen_behaviour_sha256"]
        or sha256_file(bdat_path) != accepted["frozen_behaviour_sha256"]
    ):
        raise RuntimeError("accepted_frozen_behaviour_hash_mismatch")

    lineage = pd.read_parquet(lineage_path)
    roi = pd.read_parquet(roi_path)
    if (
        len(lineage) != 8_798
        or lineage["canonical_event_key"].nunique() != 8_798
        or len(roi) != 114_374
        or not np.isfinite(roi[["value_db", "value_log_power"]]).all().all()
    ):
        raise RuntimeError("accepted_feature_table_surface_mismatch")
    return {
        "root": feature_root,
        "manifest": manifest,
        "manifest_path": manifest_path,
        "lineage": lineage,
        "lineage_path": lineage_path,
        "roi": roi,
        "roi_path": roi_path,
        "validation_path": validation_path,
        "feature_source_paths": [manifest_path, *output_paths],
        "tfr_saved_snapshot": feature_tfr_immutability["after"],
        "tfr_current_snapshot": tfr_current,
        "bdat_path": bdat_path,
    }


def code_authority(repo_root: Path, config_sha256: str) -> tuple[dict[str, Any], str]:
    """Return tracked stage-18 code descriptors and their combined fingerprint."""

    paths = [
        repo_root / CONFIG_PATH,
        repo_root / "analysis/eeg_mne/model_table_construction.py",
        Path(__file__).resolve(),
    ]
    descriptors = {
        path.relative_to(repo_root).as_posix(): {
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in paths
    }
    seed = {"config_sha256": config_sha256, "code_files": descriptors}
    return descriptors, json_fingerprint(seed)


def source_snapshots(authority: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    """Capture feature, TFR, and frozen-behaviour source evidence."""

    return {
        "feature": hashed_snapshot(authority["feature_source_paths"], root=repo_root),
        "tfr": size_mtime_snapshot_from_saved(
            authority["tfr_saved_snapshot"], repo_root=repo_root
        ),
        "behaviour": hashed_snapshot([authority["bdat_path"]], root=repo_root),
    }


def write_summary_markdown(summary: Mapping[str, Any]) -> str:
    """Return concise, neutral stage-18 production summary Markdown."""

    return f"""# DEMI EEG model tables v1 summary

Status: **{summary['status']}**  
Completed: {summary['completed_at']}

- Accepted trials: {summary['accepted_trials']:,}
- Participants/recordings: {summary['participants']} / {summary['recordings']}
- General table rows: {summary['general_rows']:,}
- Primary blocks-1--5 trials: {summary['primary_trials']:,}
- Final-overt bridge trials: {summary['bridge_trials']:,}
- Primary estimands registered: {summary['primary_estimand_count']}
- Supplementary estimands registered: {summary['supplementary_estimand_count']}
- Finite primary-overt objective-error trials: {summary['primary_finite_error_trials']:,}
- Source ROI value identity: exact
- Source feature/TFR/behaviour mutation: none

The model-ready table contains exactly five predeclared ROI/window roles per
accepted trial. Accuracy rating uses the frozen rating-point within-/between-
participant decomposition. Objective error is available through bounded
overt-only theta and alpha views. No prior, model, estimate, contrast,
interval, test, association, sensitivity, CSD derivative, GAMM, or scientific
interpretation was produced.
"""


def write_run_log(lines: Sequence[str]) -> str:
    """Return a deterministic newline-terminated production log."""

    return "\n".join(lines) + "\n"


def validate_persisted_output(
    root: Path, manifest: Mapping[str, Any]
) -> dict[str, Any]:
    """Hash/reopen every declared output and check the persisted core surface."""

    for descriptor in manifest["outputs"].values():
        validate_descriptor(root / descriptor["relative_path"], descriptor)
    table = pd.read_parquet(root / "model_ready_general.parquet")
    validation = read_json(root / VALIDATION_JSON)
    estimands = pd.read_parquet(root / "registries/estimand_registry.parquet")
    ranks = pd.read_parquet(root / RANK_PARQUET)
    if (
        len(table) != EXPECTED_GENERAL_ROWS
        or table["model_row_key"].nunique() != EXPECTED_GENERAL_ROWS
        or validation.get("status") != "pass"
        or len(estimands.loc[estimands["primary_family_member"]]) != 9
        or ranks["rank"].tolist() != [6, 6, 8]
    ):
        raise RuntimeError("persisted_model_table_reopen_validation_failed")
    for name, expected in manifest["summary"]["view_rows"].items():
        if pq.read_metadata(root / f"views/{name}.parquet").num_rows != int(expected):
            raise RuntimeError(f"persisted_model_view_row_count_mismatch:{name}")
    return {
        "general_rows": len(table),
        "view_count": len(manifest["summary"]["view_rows"]),
        "primary_estimands": 9,
        "rank": [6, 6, 8],
    }


def complete_output_reusable(
    output_root: Path,
    authority_fingerprint: str,
    current_sources: Mapping[str, Any],
) -> dict[str, Any] | None:
    """Validate an unchanged complete output without rewriting any artifact."""

    manifest_path = output_root / RUN_MANIFEST
    if not manifest_path.is_file():
        return None
    manifest = read_json(manifest_path)
    if (
        manifest.get("status") != "complete"
        or manifest.get("authority_fingerprint") != authority_fingerprint
    ):
        return None
    immutability = read_json(output_root / SOURCE_IMMUTABILITY_JSON)
    if not (
        immutability["feature_before"] == immutability["feature_after"]
        == current_sources["feature"]
        and immutability["tfr_before"] == immutability["tfr_after"]
        == current_sources["tfr"]
        and immutability["behaviour_before"] == immutability["behaviour_after"]
        == current_sources["behaviour"]
    ):
        raise RuntimeError("current_source_snapshot_differs_from_completed_model_table_run")
    validate_persisted_output(output_root, manifest)
    result = dict(manifest)
    result["invocation_cache_action"] = "complete_output_reused_after_hash_reopen_validation"
    return result


def persist_frame(
    staging: Path,
    relative_path: str,
    frame: pd.DataFrame,
    outputs: dict[str, dict[str, Any]],
    label: str,
) -> None:
    """Atomically persist and describe one staged Parquet frame."""

    path = staging / relative_path
    atomic_write_parquet(path, frame)
    outputs[label] = frame_descriptor(path, frame, root=staging)


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute preflight, construction, validation, atomic publication, or reuse."""

    repo_root = repo_root_from_script()
    os.chdir(repo_root)
    output_root = repo_root / OUTPUT_ROOT
    require_ignored_output_root(repo_root, output_root)
    config, config_sha256 = load_and_validate_model_table_config(repo_root / CONFIG_PATH)
    source = load_feature_authority(repo_root, config)
    source_before = source_snapshots(source, repo_root)
    code_files, code_fingerprint = code_authority(repo_root, config_sha256)
    authority_seed = {
        "schema_version": MODEL_TABLE_SCHEMA_VERSION,
        "model_table_namespace": MODEL_TABLE_NAMESPACE,
        "model_table_config_sha256": config_sha256,
        "model_table_code_fingerprint": code_fingerprint,
        "feature_authority_fingerprint": source["manifest"]["authority_fingerprint"],
        "feature_run_manifest_sha256": sha256_file(source["manifest_path"]),
        "roi_features_sha256": sha256_file(source["roi_path"]),
        "behavioural_lineage_sha256": sha256_file(source["lineage_path"]),
        "feature_validation_sha256": sha256_file(source["validation_path"]),
        "tfr_run_manifest_sha256": config["accepted_source"]["tfr_run_manifest_sha256"],
        "frozen_behaviour_sha256": sha256_file(source["bdat_path"]),
    }
    authority_fingerprint = json_fingerprint(authority_seed)

    reusable = complete_output_reusable(output_root, authority_fingerprint, source_before)
    if reusable is not None:
        print(
            f"PASS model_tables_v1 unchanged reuse rows={reusable['summary']['general_rows']} "
            f"views={len(reusable['summary']['view_rows'])}",
            flush=True,
        )
        return reusable
    if args.verify_current:
        raise RuntimeError("no_unchanged_complete_model_table_output_to_verify")

    archived = archive_existing_output(output_root)
    staging = make_staging_directory(output_root)
    started = time.perf_counter()
    started_at = utc_now()
    log_lines = [
        f"{started_at} stage18_start authority={authority_fingerprint}",
        "preflight=pass feature_hashes=pass tfr_immutability=pass behaviour_hash=pass",
    ]
    atomic_write_json(
        staging / RUN_STATE,
        {
            "schema_version": MODEL_TABLE_SCHEMA_VERSION,
            "model_table_namespace": MODEL_TABLE_NAMESPACE,
            "status": "running",
            "started_at": started_at,
            "authority_fingerprint": authority_fingerprint,
        },
    )
    try:
        provenance = {
            "model_table_config_sha256": config_sha256,
            "model_table_code_fingerprint": code_fingerprint,
            "source_feature_authority_fingerprint": source["manifest"]["authority_fingerprint"],
            "source_feature_run_manifest_sha256": sha256_file(source["manifest_path"]),
            "source_feature_roi_features_sha256": sha256_file(source["roi_path"]),
            "source_feature_behavioural_lineage_sha256": sha256_file(source["lineage_path"]),
            "source_feature_validation_sha256": sha256_file(source["validation_path"]),
            "source_tfr_run_manifest_sha256_stage18": config["accepted_source"]["tfr_run_manifest_sha256"],
            "source_frozen_behaviour_sha256": sha256_file(source["bdat_path"]),
        }
        table, predictor_facts, selected_source_rows = build_general_table(
            source["lineage"], source["roi"], config, provenance=provenance
        )
        views = deterministic_views(table)
        roles = model_role_registry(config)
        factors = factor_coding_registry(config)
        estimands = estimand_registry()
        reasons = controlled_reason_registry()
        ranks = model_matrix_rank_validation(views)
        validation = validate_general_and_views(
            table, views, selected_source_rows, predictor_facts, estimands, ranks
        )
        counts = count_summary(table, views)
        log_lines.extend(
            [
                f"general_table_rows={len(table)} accepted_trials={table['canonical_event_key'].nunique()}",
                f"views_validated={len(views)} primary_estimands=9 supplementary_estimands=3",
                "source_roi_value_identity=exact model_matrix_rank=6,6,8",
            ]
        )

        outputs: dict[str, dict[str, Any]] = {}
        persist_frame(staging, "model_ready_general.parquet", table, outputs, "general_table")
        for name, frame in views.items():
            persist_frame(staging, f"views/{name}.parquet", frame, outputs, f"view_{name}")
        persist_frame(staging, "registries/model_role_registry.parquet", roles, outputs, "model_role_registry")
        persist_frame(staging, "registries/estimand_registry.parquet", estimands, outputs, "estimand_registry")
        persist_frame(staging, "registries/factor_coding_registry.parquet", factors, outputs, "factor_coding_registry")
        persist_frame(staging, "registries/controlled_reason_registry.parquet", reasons, outputs, "controlled_reason_registry")
        persist_frame(staging, "summaries/model_table_counts.parquet", counts, outputs, "count_summary")
        persist_frame(staging, RANK_PARQUET, ranks, outputs, "model_matrix_rank")

        config_snapshot = staging / "configuration/model_table_config_v1.yaml"
        atomic_write_text(config_snapshot, (repo_root / CONFIG_PATH).read_text(encoding="utf-8"))
        if sha256_file(config_snapshot) != config_sha256:
            raise RuntimeError("model_table_config_snapshot_hash_mismatch")
        outputs["config_snapshot"] = file_descriptor(config_snapshot, root=staging)

        source_after = source_snapshots(source, repo_root)
        if source_after != source_before:
            raise RuntimeError("stage18_source_artifact_changed_during_construction")
        immutability = {
            "status": "unchanged",
            "feature_before": source_before["feature"],
            "feature_after": source_after["feature"],
            "tfr_before": source_before["tfr"],
            "tfr_after": source_after["tfr"],
            "behaviour_before": source_before["behaviour"],
            "behaviour_after": source_after["behaviour"],
            "feature_file_count": len(source_before["feature"]),
            "tfr_array_count": len(source_before["tfr"]),
            "tfr_total_bytes": sum(int(item["size_bytes"]) for item in source_before["tfr"]),
        }
        immutability_path = staging / SOURCE_IMMUTABILITY_JSON
        atomic_write_json(immutability_path, immutability)
        outputs["source_immutability"] = file_descriptor(immutability_path, root=staging)

        validation.update(
            {
                "source_feature_unchanged": True,
                "source_tfr_unchanged": True,
                "source_behaviour_unchanged": True,
                "atomic_directory_publication": True,
                "validated_unchanged_reuse_supported": True,
                "forbidden_outputs_absent": [
                    "priors", "prior_predictive", "models", "fits", "estimates",
                    "contrasts", "intervals", "tests", "associations", "influence",
                    "gamm", "csd", "interpretation",
                ],
            }
        )
        validation_path = staging / VALIDATION_JSON
        atomic_write_json(validation_path, validation)
        outputs["validation"] = file_descriptor(validation_path, root=staging)

        completed_at = utc_now()
        runtime = time.perf_counter() - started
        summary = {
            "status": "complete",
            "completed_at": completed_at,
            "runtime_seconds": runtime,
            "accepted_trials": 8_798,
            "participants": 81,
            "recordings": 81,
            "general_rows": len(table),
            "primary_trials": 7_307,
            "bridge_trials": 738,
            "primary_estimand_count": 9,
            "supplementary_estimand_count": 3,
            "primary_finite_error_trials": 3_703,
            "view_rows": {name: len(frame) for name, frame in views.items()},
            "source_value_identity": "exact_persisted_roi_rows",
            "explicit_nonoperations": [
                "no TFR feature reconstruction", "no eligibility change",
                "no predictor standardization or binning", "no objective-error imputation",
                "no priors or prior-predictive simulation", "no model compilation or fitting",
                "no estimand, interval, test, or association calculation",
                "no sensitivity, influence, GAMM, CSD, or interpretation",
            ],
        }
        summary_json_path = staging / SUMMARY_JSON
        summary_md_path = staging / SUMMARY_MD
        atomic_write_json(summary_json_path, summary)
        atomic_write_text(summary_md_path, write_summary_markdown(summary))
        outputs["summary_json"] = file_descriptor(summary_json_path, root=staging)
        outputs["summary_markdown"] = file_descriptor(summary_md_path, root=staging)
        log_lines.append(f"{completed_at} stage18_complete runtime_seconds={runtime:.3f}")
        log_path = staging / RUN_LOG
        atomic_write_text(log_path, write_run_log(log_lines))
        outputs["production_log"] = file_descriptor(log_path, root=staging)

        run_manifest = {
            "schema_version": MODEL_TABLE_SCHEMA_VERSION,
            "model_table_namespace": MODEL_TABLE_NAMESPACE,
            "status": "complete",
            "started_at": started_at,
            "completed_at": completed_at,
            "authority_fingerprint": authority_fingerprint,
            "authority": authority_seed,
            "code_files": code_files,
            "sources": {
                "feature_root": FEATURE_ROOT.as_posix(),
                "feature_run_manifest_sha256": sha256_file(source["manifest_path"]),
                "roi_features_sha256": sha256_file(source["roi_path"]),
                "behavioural_lineage_sha256": sha256_file(source["lineage_path"]),
                "feature_validation_sha256": sha256_file(source["validation_path"]),
                "tfr_run_manifest_sha256": config["accepted_source"]["tfr_run_manifest_sha256"],
                "frozen_behaviour_sha256": sha256_file(source["bdat_path"]),
                "immutability_path": SOURCE_IMMUTABILITY_JSON,
            },
            "outputs": outputs,
            "summary": summary,
            "archive_of_previous_nonreusable_output": (
                archived.relative_to(repo_root).as_posix() if archived is not None else None
            ),
            "software": {
                "python": platform.python_version(),
                "pandas": pd.__version__,
                "numpy": np.__version__,
                "platform": platform.platform(),
            },
        }
        atomic_write_json(staging / RUN_MANIFEST, run_manifest)
        atomic_write_json(
            staging / RUN_STATE,
            {
                "schema_version": MODEL_TABLE_SCHEMA_VERSION,
                "model_table_namespace": MODEL_TABLE_NAMESPACE,
                "status": "complete",
                "started_at": started_at,
                "completed_at": completed_at,
                "authority_fingerprint": authority_fingerprint,
                "general_rows": len(table),
                "view_count": len(views),
                "runtime_seconds": runtime,
            },
        )
        validate_persisted_output(staging, run_manifest)
        final_sources = source_snapshots(source, repo_root)
        if final_sources != source_before:
            raise RuntimeError("stage18_source_artifact_changed_before_publication")
        publish_staging_directory(staging, output_root)
        validate_persisted_output(output_root, run_manifest)
        print(
            f"PASS model_tables_v1 rows={len(table)} views={len(views)} "
            f"runtime={runtime:.2f}s",
            flush=True,
        )
        return run_manifest
    except Exception as error:
        failure = {
            "schema_version": MODEL_TABLE_SCHEMA_VERSION,
            "model_table_namespace": MODEL_TABLE_NAMESPACE,
            "status": "failed",
            "failed_at": utc_now(),
            "error_type": type(error).__name__,
            "error": str(error),
            "staging_preserved": staging.as_posix(),
        }
        try:
            atomic_write_json(staging / RUN_STATE, failure)
        finally:
            raise


def main() -> None:
    """Run stage 18 and allow fail-closed exceptions to surface."""

    run(parse_args())


if __name__ == "__main__":
    main()
