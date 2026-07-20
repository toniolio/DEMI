#!/usr/bin/env python3
"""Construct the accepted DEMI trial-level Morlet power and matched dB surface.

This is the production stage-16 driver.  It consumes only the validated
``epochs_v1`` response-onset, response-end, and ``red_on`` manifests/shards.
For one recording at a time it verifies every source hash and canonical key,
computes non-decisional transient diagnostics, resamples an in-memory 30-scalp
copy to 100 Hz, computes trial-level Morlet power, derives the trial-matched
visual pre-response baseline, and atomically publishes raw and dB NumPy arrays.

Inputs:
    ``_Data/eeg/epochs_v1/epoch_construction_run_manifest.json``, all three
    accepted family manifests and metadata tables, and the linked accepted
    per-recording Epochs FIF/manifest pairs.

Outputs:
    Ignored ``_Data/eeg/tfr_v1`` recording shards, float32 NumPy arrays,
    Parquet/CSV metadata and diagnostics, axes/configuration snapshots,
    per-recording/per-family/aggregate manifests, validation summaries, source
    immutability evidence, run state, and concise Markdown.

This stage explicitly does not modify accepted Epochs, interpolate or reject
epochs, run AutoReject, compute CSD/complex/phase/connectivity/ROI/band-power
products, decide participant inclusion, fit models, or interpret results.

Safety:
    Every accepted source must match its epoch manifest by size and SHA-256.
    Source size/mtime snapshots are compared before and after processing.
    Output uses atomic recording-directory publication and hash/reopen cache
    validation.  Invalid identities, baselines, power, dB, shapes, dtypes, or
    stale cache contents fail closed.

Run from the repository root:

    tools/run_tfr_v1.sh

Use ``--recording demi_13_data --benchmark`` for the required bounded benchmark.
Use an unchanged full invocation to validate and reuse all completed shards.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import platform
import resource
import threading
import time
from typing import Any, Callable, Mapping, TypeVar

import mne
import numpy as np
import pandas as pd
import scipy

from tfr_construction import (
    ALL_FAMILIES,
    ArrayContract,
    BASELINE_TMAX,
    BASELINE_TMIN,
    EXPECTED_DURATION_WARNING_COUNT,
    EXPECTED_FILE49_COUNT,
    EXPECTED_RECORDING_COUNT,
    EXPECTED_STRICT_CLEAN_COUNT,
    EXPECTED_TASK_SAMPLE_COUNT,
    EXPECTED_TASK_TIMES,
    EXPECTED_TRIAL_COUNT,
    FREQUENCIES_HZ,
    N_CYCLES,
    PERSISTED_DTYPE,
    PRIMARY_SCALP_CHANNELS,
    TARGET_SFREQ,
    TASK_FAMILIES,
    TASK_TMAX,
    TASK_TMIN,
    TFR_NAMESPACE,
    TFR_SCHEMA_VERSION,
    artifact_diagnostics,
    assess_cached_recording,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_npy,
    atomic_write_parquet,
    atomic_write_text,
    compare_source_snapshots,
    compute_morlet_power,
    content_fingerprint,
    crop_task_power,
    derive_baseline_mean,
    expected_wavelet_support,
    forbidden_output_scan,
    load_and_validate_config,
    normalize_db,
    ordered_key_hash,
    recording_fingerprint,
    reopen_validate_npy,
    require_ignored_output_root,
    resample_for_tfr,
    sha256_file,
    source_snapshot,
    stable_frame_hash,
    validate_family_identity,
    validate_stored_formula,
)


REPO_EPOCH_ROOT = Path("_Data/eeg/epochs_v1")
OUTPUT_ROOT = Path("_Data/eeg/tfr_v1")
CONFIG_PATH = Path("analysis/eeg_mne/tfr_config_v1.yaml")
RUN_MANIFEST_NAME = "tfr_run_manifest.json"
RUN_STATE_NAME = "tfr_run_state.json"
SUMMARY_JSON_NAME = "tfr_summary.json"
SUMMARY_MD_NAME = "tfr_summary.md"
SOURCE_IMMUTABILITY_NAME = "accepted_epoch_source_immutability.json"
T = TypeVar("T")


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def repo_root_from_script() -> Path:
    """Return and validate the repository root containing this script."""

    root = Path(__file__).resolve().parents[2]
    if not (root / ".git").exists() or not (root / CONFIG_PATH).is_file():
        raise RuntimeError(f"could not locate DEMI repository root from {__file__}")
    return root


def parse_args() -> argparse.Namespace:
    """Parse bounded recording, benchmark, force, and verification options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--recording",
        help="Process one exact recording stem or accepted source filename.",
    )
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Record runtime/memory/storage projection; requires --recording.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Archive and rebuild one exact --recording result.",
    )
    parser.add_argument(
        "--verify-current",
        action="store_true",
        help="Validate existing selected outputs and fail if any need construction.",
    )
    args = parser.parse_args()
    if args.benchmark and not args.recording:
        parser.error("--benchmark requires --recording")
    if args.force and not args.recording:
        parser.error("--force requires one exact --recording")
    if args.force and args.verify_current:
        parser.error("--force and --verify-current are incompatible")
    return args


def print_progress(message: str) -> None:
    """Print one timestamped, immediately flushed progress line."""

    print(f"[{utc_now()}] {message}", flush=True)


def with_heartbeat(
    label: str, operation: Callable[[], T], interval: float = 30.0
) -> T:
    """Run a blocking operation while a daemon prints elapsed heartbeats."""

    stopped = threading.Event()
    started = time.perf_counter()

    def heartbeat() -> None:
        while not stopped.wait(interval):
            print_progress(f"heartbeat {label}; elapsed={time.perf_counter() - started:.1f}s")

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        return operation()
    finally:
        stopped.set()
        thread.join(timeout=1.0)


def read_json(path: Path) -> dict[str, Any]:
    """Read one required JSON object."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON authority is not an object: {path}")
    return value


def descriptor_for_file(path: Path, *, root: Path, row_count: int | None = None) -> dict[str, Any]:
    """Return a deterministic output-relative file descriptor."""

    descriptor: dict[str, Any] = {
        "relative_path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if row_count is not None:
        descriptor["row_count"] = int(row_count)
    return descriptor


def load_epoch_authority(repo_root: Path) -> dict[str, Any]:
    """Load and validate accepted epoch manifests, metadata, config, and code."""

    epoch_root = repo_root / REPO_EPOCH_ROOT
    run_manifest_path = epoch_root / "epoch_construction_run_manifest.json"
    run_manifest = read_json(run_manifest_path)
    if run_manifest.get("status") != "complete":
        raise RuntimeError("accepted epoch run manifest is not complete")
    config, config_sha256 = load_and_validate_config(repo_root / CONFIG_PATH)
    family_manifests: dict[str, dict[str, Any]] = {}
    family_manifest_paths: dict[str, Path] = {}
    global_metadata: dict[str, pd.DataFrame] = {}
    for family in ALL_FAMILIES:
        path = epoch_root / "manifests" / f"{family}_manifest.json"
        manifest = read_json(path)
        if (
            manifest.get("status") != "complete"
            or int(manifest.get("epoch_count", -1)) != EXPECTED_TRIAL_COUNT
            or int(manifest.get("recording_shard_count", -1)) != EXPECTED_RECORDING_COUNT
        ):
            raise RuntimeError(f"accepted family manifest differs: {family}")
        metadata_descriptor = manifest["metadata"]
        metadata_path = epoch_root / metadata_descriptor["path"]
        if (
            metadata_path.stat().st_size != int(metadata_descriptor["size_bytes"])
            or sha256_file(metadata_path) != metadata_descriptor["sha256"]
        ):
            raise RuntimeError(f"accepted family metadata hash differs: {family}")
        frame = pd.read_parquet(metadata_path)
        if len(frame) != EXPECTED_TRIAL_COUNT:
            raise RuntimeError(f"accepted family metadata row count differs: {family}")
        if stable_frame_hash(frame) != metadata_descriptor["content_sha256"]:
            raise RuntimeError(f"accepted family metadata content differs: {family}")
        family_manifests[family] = manifest
        family_manifest_paths[family] = path
        global_metadata[family] = frame
    common_key_sha256 = validate_family_identity(global_metadata)
    if common_key_sha256 != run_manifest["families"][0].get(
        "ordered_key_sha256", run_manifest.get("common_ordered_key_sha256")
    ):
        # Stage 15 stores the common digest in its summary and family manifests,
        # not necessarily in each run-manifest family summary.
        summary = read_json(epoch_root / "epoch_construction_summary.json")
        if common_key_sha256 != summary["common_ordered_key_sha256"]:
            raise RuntimeError("accepted common canonical-key hash differs")

    orders = [
        [str(shard["recording"]) for shard in family_manifests[family]["shards"]]
        for family in ALL_FAMILIES
    ]
    if not all(order == orders[0] for order in orders[1:]):
        raise RuntimeError("accepted family recording order differs")

    recordings: dict[str, dict[str, Any]] = {}
    for index, recording in enumerate(orders[0]):
        source_descriptors: dict[str, dict[str, Any]] = {}
        source_filename: str | None = None
        row_count: int | None = None
        for family in ALL_FAMILIES:
            shard = family_manifests[family]["shards"][index]
            artifact = shard["artifact"]
            path = epoch_root / artifact["relative_path"]
            shard_manifest_path = path.parent / "manifest.json"
            shard_manifest = read_json(shard_manifest_path)
            if shard_manifest.get("status") != "complete":
                raise RuntimeError(f"accepted shard manifest not complete: {shard_manifest_path}")
            if shard_manifest["artifact"]["sha256"] != artifact["sha256"]:
                raise RuntimeError(f"accepted shard manifest chain differs: {path}")
            if source_filename is None:
                source_filename = str(shard["source_filename"])
                row_count = int(shard["row_count"])
            elif source_filename != str(shard["source_filename"]) or row_count != int(
                shard["row_count"]
            ):
                raise RuntimeError(f"accepted family shard identity differs: {recording}")
            source_descriptors[family] = {
                "family": family,
                "path": path,
                "relative_path": path.relative_to(repo_root).as_posix(),
                "size_bytes": int(artifact["size_bytes"]),
                "sha256": str(artifact["sha256"]),
                "manifest_path": shard_manifest_path,
                "manifest_relative_path": shard_manifest_path.relative_to(repo_root).as_posix(),
                "manifest_sha256": sha256_file(shard_manifest_path),
                "fingerprint": str(shard["fingerprint"]),
            }
        assert source_filename is not None and row_count is not None
        recordings[recording] = {
            "recording": recording,
            "source_filename": source_filename,
            "row_count": row_count,
            "sources": source_descriptors,
        }

    code_paths = (
        repo_root / "analysis/eeg_mne/tfr_construction.py",
        repo_root / "analysis/eeg_mne/16_construct_trial_level_tfr.py",
    )
    code_files = {
        path.relative_to(repo_root).as_posix(): sha256_file(path) for path in code_paths
    }
    authority_seed = {
        "epoch_run_manifest_sha256": sha256_file(run_manifest_path),
        "family_manifest_sha256": {
            family: sha256_file(path) for family, path in family_manifest_paths.items()
        },
        "family_metadata_sha256": {
            family: family_manifests[family]["metadata"]["sha256"]
            for family in ALL_FAMILIES
        },
        "config_sha256": config_sha256,
        "code_files": code_files,
        "software": {
            "python": platform.python_version(),
            "mne": mne.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    return {
        "epoch_root": epoch_root,
        "run_manifest_path": run_manifest_path,
        "run_manifest": run_manifest,
        "family_manifests": family_manifests,
        "family_manifest_paths": family_manifest_paths,
        "global_metadata": global_metadata,
        "recording_order": orders[0],
        "recordings": recordings,
        "common_key_sha256": common_key_sha256,
        "config": config,
        "config_sha256": config_sha256,
        "code_files": code_files,
        "code_sha256": content_fingerprint(code_files),
        "software": authority_seed["software"],
        "authority_fingerprint": content_fingerprint(authority_seed),
    }


def select_recordings(
    authority: Mapping[str, Any], requested: str | None
) -> list[dict[str, Any]]:
    """Return all accepted recordings or one exact stem/source filename."""

    recordings = authority["recordings"]
    if requested is None:
        return [recordings[name] for name in authority["recording_order"]]
    matches = [
        item
        for item in recordings.values()
        if requested in {item["recording"], item["source_filename"]}
    ]
    if len(matches) != 1:
        raise RuntimeError(f"--recording must identify exactly one accepted shard: {requested}")
    return matches


def validate_source_file(source: Mapping[str, Any]) -> dict[str, Any]:
    """Require one accepted Epochs FIF to match its recorded size and SHA-256."""

    path = Path(source["path"])
    if not path.is_file() or path.stat().st_size != int(source["size_bytes"]):
        raise RuntimeError(f"accepted_epoch_source_size_mismatch:{path}")
    observed = with_heartbeat(f"hash accepted source {path.name}", lambda: sha256_file(path))
    if observed != source["sha256"]:
        raise RuntimeError(f"accepted_epoch_source_hash_mismatch:{path}")
    return {
        "family": source["family"],
        "path": source["relative_path"],
        "size_bytes": int(source["size_bytes"]),
        "sha256": observed,
        "manifest_path": source["manifest_relative_path"],
        "manifest_sha256": source["manifest_sha256"],
        "matches_accepted_manifest": True,
    }


def read_family_metadata(recording: Mapping[str, Any]) -> dict[str, pd.DataFrame]:
    """Read metadata-only Epochs objects and validate source/header contracts."""

    result: dict[str, pd.DataFrame] = {}
    for family in ALL_FAMILIES:
        source = recording["sources"][family]
        epochs = mne.read_epochs(source["path"], preload=False, verbose="ERROR")
        if (
            epochs.metadata is None
            or len(epochs) != recording["row_count"]
            or not np.isclose(epochs.info["sfreq"], 1000.0, atol=1e-9, rtol=0.0)
            or epochs.baseline is not None
            or epochs.reject is not None
            or epochs.flat is not None
            or any(epochs.drop_log)
        ):
            raise RuntimeError(f"accepted_epoch_object_contract_mismatch:{family}")
        frame = epochs.metadata.reset_index(drop=True).copy()
        if frame["derived_anchor_type"].astype(str).nunique() != 1 or str(
            frame["derived_anchor_type"].iloc[0]
        ) != family:
            raise RuntimeError(f"accepted_epoch_anchor_family_mismatch:{family}")
        result[family] = frame
        del epochs
    validate_family_identity(result)
    return result


def expected_array_contracts(row_count: int) -> dict[str, ArrayContract]:
    """Return exact persisted array contracts for one recording."""

    task_shape = (
        row_count,
        len(PRIMARY_SCALP_CHANNELS),
        len(FREQUENCIES_HZ),
        EXPECTED_TASK_SAMPLE_COUNT,
    )
    baseline_shape = task_shape[:3]
    task_axes = ("trial", "channel", "frequency", "time")
    return {
        "response_onset_raw": ArrayContract(
            task_shape, PERSISTED_DTYPE, task_axes, positive=True
        ),
        "response_onset_db": ArrayContract(task_shape, PERSISTED_DTYPE, task_axes),
        "response_end_raw": ArrayContract(
            task_shape, PERSISTED_DTYPE, task_axes, positive=True
        ),
        "response_end_db": ArrayContract(task_shape, PERSISTED_DTYPE, task_axes),
        "red_on_baseline_mean": ArrayContract(
            baseline_shape,
            PERSISTED_DTYPE,
            ("trial", "channel", "frequency"),
            positive=True,
        ),
    }


def recording_cache_fingerprint(
    authority: Mapping[str, Any],
    recording: Mapping[str, Any],
    metadata: Mapping[str, pd.DataFrame],
) -> str:
    """Build the complete deterministic per-recording cache fingerprint."""

    return recording_fingerprint(
        {
            "schema_version": TFR_SCHEMA_VERSION,
            "namespace": TFR_NAMESPACE,
            "authority_fingerprint": authority["authority_fingerprint"],
            "code_sha256": authority["code_sha256"],
            "config_sha256": authority["config_sha256"],
            "recording": recording["recording"],
            "source_epoch_sha256": {
                family: recording["sources"][family]["sha256"] for family in ALL_FAMILIES
            },
            "source_epoch_manifest_sha256": {
                family: recording["sources"][family]["manifest_sha256"]
                for family in ALL_FAMILIES
            },
            "metadata_content_sha256": {
                family: stable_frame_hash(metadata[family]) for family in ALL_FAMILIES
            },
            "ordered_key_sha256": validate_family_identity(metadata),
            "channels": list(PRIMARY_SCALP_CHANNELS),
            "frequencies_hz": FREQUENCIES_HZ.tolist(),
            "n_cycles": N_CYCLES.tolist(),
            "task_times_seconds": EXPECTED_TASK_TIMES.tolist(),
            "baseline_interval_seconds": [BASELINE_TMIN, BASELINE_TMAX],
            "software": authority["software"],
        }
    )


def archive_previous_recording(
    recording_dir: Path, output_root: Path, *, reason: str
) -> list[str]:
    """Move a stale generated recording directory into ignored history."""

    if not recording_dir.exists():
        return []
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    history = output_root / "history" / stamp / recording_dir.name
    history.parent.mkdir(parents=True, exist_ok=True)
    os.replace(recording_dir, history)
    atomic_write_json(
        history / "archive_reason.json",
        {"archived_at": utc_now(), "reason": reason},
    )
    return [history.relative_to(output_root).as_posix()]


def archive_incomplete_temporaries(output_root: Path, recording: str) -> list[str]:
    """Preserve abandoned atomic-build directories before a new attempt."""

    moved = []
    for path in sorted((output_root / "recordings").glob(f".{recording}.tmp-*")):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        target = output_root / "history" / stamp / f"{recording}_incomplete"
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(path, target)
        atomic_write_json(
            target / "archive_reason.json",
            {"archived_at": utc_now(), "reason": "incomplete_atomic_build_recovered"},
        )
        moved.append(target.relative_to(output_root).as_posix())
    return moved


def write_frame_descriptor(
    path: Path, frame: pd.DataFrame, *, root: Path
) -> dict[str, Any]:
    """Atomically write Parquet and return hash/content provenance."""

    atomic_write_parquet(path, frame)
    descriptor = descriptor_for_file(path, root=root, row_count=len(frame))
    descriptor["content_sha256"] = stable_frame_hash(frame)
    descriptor["format"] = "parquet"
    return descriptor


def process_family(
    source: Mapping[str, Any],
    family: str,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Load one family, calculate diagnostics, resample, and compute power."""

    epochs = with_heartbeat(
        f"load {source['path']}",
        lambda: mne.read_epochs(source["path"], preload=True, verbose="ERROR"),
    )
    diagnostics = artifact_diagnostics(epochs, family)
    work = with_heartbeat(
        f"resample {family}",
        lambda: resample_for_tfr(epochs),
    )
    del epochs
    gc.collect()
    power, times = with_heartbeat(
        f"Morlet convolution {family}",
        lambda: compute_morlet_power(work),
    )
    del work
    gc.collect()
    return power, times, diagnostics


def build_recording(
    *,
    repo_root: Path,
    output_root: Path,
    authority: Mapping[str, Any],
    recording: Mapping[str, Any],
    metadata: Mapping[str, pd.DataFrame],
    fingerprint: str,
    force: bool,
) -> dict[str, Any]:
    """Build, reopen, validate, and atomically publish one recording result."""

    recording_name = recording["recording"]
    final_dir = output_root / "recordings" / recording_name
    contracts = expected_array_contracts(recording["row_count"])
    assessment = assess_cached_recording(final_dir, fingerprint, contracts)
    if assessment["action"] == "reuse" and not force:
        reused = dict(assessment["manifest"])
        reused["invocation_cache_action"] = "reused_after_hash_reopen_validation"
        return reused
    reason = "explicit_force_rebuild" if force else assessment["reason"]
    recovered = archive_incomplete_temporaries(output_root, recording_name)
    temp_dir = final_dir.parent / f".{recording_name}.tmp-{os.getpid()}-{time.time_ns()}"
    temp_dir.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    peak_before = peak_rss_bytes()
    print_progress(
        f"constructing recording={recording_name} rows={recording['row_count']} "
        f"cache_reason={reason}"
    )

    arrays: dict[str, dict[str, Any]] = {}
    tables: dict[str, dict[str, Any]] = {}
    validations: dict[str, Any] = {}
    diagnostics_parts: list[pd.DataFrame] = []

    red_power, red_times, red_diagnostics = process_family(
        recording["sources"]["red_on"], "red_on"
    )
    diagnostics_parts.append(red_diagnostics)
    baseline_mean, baseline_validation = derive_baseline_mean(red_power, red_times)
    del red_power, red_times
    gc.collect()
    baseline32 = baseline_mean.astype(np.float32)
    baseline_path = temp_dir / "red_on" / "baseline_mean_power.npy"
    baseline_descriptor = atomic_write_npy(baseline_path, baseline32)
    baseline_descriptor["relative_path"] = baseline_path.relative_to(temp_dir).as_posix()
    baseline_descriptor["axis_order"] = ["trial", "channel", "frequency"]
    arrays["red_on_baseline_mean"] = baseline_descriptor
    del baseline32

    baseline_validity = metadata["red_on"][
        [
            "canonical_order_index",
            "canonical_event_key",
            "behavioural_id",
            "eeg_source_id",
            "source_recording_filename",
            "condition_semantics",
            "strict_clean_eligibility",
            "duration_warning_flag",
            "continuous_qc_warning",
        ]
    ].copy()
    baseline_validity["baseline_valid"] = True
    baseline_validity["baseline_reason_code"] = "baseline_valid"
    baseline_validity["baseline_min_power"] = baseline_mean.min(axis=(1, 2))
    baseline_validity["baseline_max_power"] = baseline_mean.max(axis=(1, 2))
    baseline_validity_path = temp_dir / "red_on" / "baseline_validity.parquet"
    tables["baseline_validity"] = write_frame_descriptor(
        baseline_validity_path, baseline_validity, root=temp_dir
    )

    for family in TASK_FAMILIES:
        power, times, diagnostics = process_family(recording["sources"][family], family)
        diagnostics_parts.append(diagnostics)
        retained_power, retained_times = crop_task_power(power, times)
        del power, times
        gc.collect()
        db = normalize_db(retained_power, baseline_mean)
        raw32 = retained_power.astype(np.float32)
        db32 = db.astype(np.float32)
        del retained_power, db
        family_dir = temp_dir / family
        raw_path = family_dir / "raw_power.npy"
        db_path = family_dir / "db_power.npy"
        raw_descriptor = atomic_write_npy(raw_path, raw32)
        db_descriptor = atomic_write_npy(db_path, db32)
        raw_descriptor["relative_path"] = raw_path.relative_to(temp_dir).as_posix()
        db_descriptor["relative_path"] = db_path.relative_to(temp_dir).as_posix()
        raw_descriptor["axis_order"] = ["trial", "channel", "frequency", "time"]
        db_descriptor["axis_order"] = ["trial", "channel", "frequency", "time"]
        arrays[f"{family}_raw"] = raw_descriptor
        arrays[f"{family}_db"] = db_descriptor
        del raw32, db32
        row_metadata = metadata[family].copy()
        row_metadata.insert(0, "tfr_row_index", np.arange(len(row_metadata), dtype=int))
        row_metadata.insert(1, "tfr_family", family)
        row_metadata.insert(
            2,
            "tfr_global_row_index",
            row_metadata["canonical_order_index"].to_numpy(dtype=int),
        )
        row_metadata["tfr_primary_eligibility"] = True
        row_metadata["tfr_primary_reason_code"] = "retained_primary_no_epoch_rejection"
        row_metadata_path = family_dir / "row_metadata.parquet"
        tables[f"{family}_metadata"] = write_frame_descriptor(
            row_metadata_path, row_metadata, root=temp_dir
        )
        validations[f"{family}_retained_time_vector_seconds"] = retained_times.tolist()
        gc.collect()

    del baseline_mean
    artifact_frame = pd.concat(diagnostics_parts, ignore_index=True)
    if len(artifact_frame) != recording["row_count"] * len(ALL_FAMILIES):
        raise RuntimeError("artifact_diagnostic_row_count_mismatch")
    diagnostics_path = temp_dir / "artifact_diagnostics.parquet"
    tables["artifact_diagnostics"] = write_frame_descriptor(
        diagnostics_path, artifact_frame, root=temp_dir
    )

    common_metadata = metadata["response_onset"][
        [
            "canonical_order_index",
            "canonical_event_key",
            "behavioural_id",
            "eeg_source_id",
            "source_recording_filename",
            "source_file_role",
            "condition_semantics",
            "physical",
            "primary_eligibility",
            "strict_clean_eligibility",
            "duration_warning_flag",
            "continuous_qc_warning",
            "continuous_v2_run_id",
            "continuous_v2_manifest_path",
            "continuous_v2_manifest_id",
            "post_ica_derivative_path",
            "continuous_source_sha256",
            "event_policy_version",
            "ledger_version",
        ]
    ].copy()
    common_metadata.insert(0, "tfr_row_index", np.arange(len(common_metadata), dtype=int))
    common_metadata.insert(
        1,
        "tfr_global_row_index",
        common_metadata["canonical_order_index"].to_numpy(dtype=int),
    )
    common_metadata["tfr_primary_eligibility"] = True
    common_metadata["tfr_primary_reason_code"] = "retained_primary_no_epoch_rejection"
    common_path = temp_dir / "row_metadata.parquet"
    tables["common_row_metadata"] = write_frame_descriptor(
        common_path, common_metadata, root=temp_dir
    )

    for label, contract in contracts.items():
        path = temp_dir / arrays[label]["relative_path"]
        validations[label] = reopen_validate_npy(
            path, contract, expected_sha256=arrays[label]["sha256"]
        )
        validations[label]["path"] = path.relative_to(temp_dir).as_posix()
    for family in TASK_FAMILIES:
        validations[f"{family}_formula"] = validate_stored_formula(
            temp_dir / arrays[f"{family}_raw"]["relative_path"],
            temp_dir / arrays[f"{family}_db"]["relative_path"],
            temp_dir / arrays["red_on_baseline_mean"]["relative_path"],
        )
    validations["baseline"] = baseline_validation
    validations["canonical_key_sha256"] = validate_family_identity(metadata)
    validations["zero_trial_loss"] = all(
        int(tables[f"{family}_metadata"]["row_count"]) == recording["row_count"]
        for family in TASK_FAMILIES
    )
    if not validations["zero_trial_loss"]:
        raise RuntimeError("silent_trial_loss")

    source_provenance = {
        family: {
            "path": recording["sources"][family]["relative_path"],
            "size_bytes": recording["sources"][family]["size_bytes"],
            "sha256": recording["sources"][family]["sha256"],
            "manifest_path": recording["sources"][family]["manifest_relative_path"],
            "manifest_sha256": recording["sources"][family]["manifest_sha256"],
            "fingerprint": recording["sources"][family]["fingerprint"],
        }
        for family in ALL_FAMILIES
    }
    manifest = {
        "schema_version": TFR_SCHEMA_VERSION,
        "tfr_namespace": TFR_NAMESPACE,
        "status": "complete",
        "created_at": utc_now(),
        "fingerprint": fingerprint,
        "cache_reason": reason,
        "recovered_incomplete_paths": recovered,
        "recording": {
            "recording_stem": recording_name,
            "source_filename": recording["source_filename"],
            "row_count": recording["row_count"],
        },
        "source_epochs": source_provenance,
        "source_epoch_common_key_sha256": validate_family_identity(metadata),
        "config_sha256": authority["config_sha256"],
        "code_sha256": authority["code_sha256"],
        "software": authority["software"],
        "channels": list(PRIMARY_SCALP_CHANNELS),
        "frequencies_hz": FREQUENCIES_HZ.tolist(),
        "n_cycles": N_CYCLES.tolist(),
        "task_times_seconds": EXPECTED_TASK_TIMES.tolist(),
        "baseline": {
            "family": "red_on",
            "interval_seconds_inclusive": [BASELINE_TMIN, BASELINE_TMAX],
            "description": "trial-matched visual pre-response reference",
            "matching_key": "canonical_event_key",
            "mean_power_over_time_before_ratio": True,
            "shared_by_task_families": list(TASK_FAMILIES),
        },
        "representations": {
            "raw_power": {"formula": "Morlet power", "persisted_dtype": "float32"},
            "db_power": {
                "formula": "10 * log10(raw_power / baseline_mean_power)",
                "persisted_dtype": "float32",
            },
            "log_power_sensitivity": {
                "formula": "10 * log10(raw_power)",
                "storage": "deterministic helper/view; no separate cube",
            },
        },
        "arrays": arrays,
        "tables": tables,
        "validation": validations,
        "artifact_policy": {
            "role": "retain_and_diagnose",
            "interpolation_applied": False,
            "rejection_applied": False,
            "autoreject_applied": False,
            "diagnostics_change_eligibility": False,
        },
        "runtime_seconds": time.perf_counter() - started,
        "peak_rss_bytes_after": peak_rss_bytes(),
        "peak_rss_bytes_before": peak_before,
    }
    atomic_write_json(temp_dir / "manifest.json", manifest)
    archived = archive_previous_recording(final_dir, output_root, reason=reason)
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    os.replace(temp_dir, final_dir)
    manifest["archived_previous_paths"] = archived
    manifest["invocation_cache_action"] = "constructed_and_reopened"
    print_progress(
        f"published recording={recording_name} rows={recording['row_count']} "
        f"elapsed={time.perf_counter() - started:.1f}s"
    )
    return manifest


def peak_rss_bytes() -> int:
    """Return process peak resident-set size normalized to bytes."""

    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def validate_special_surface(metadata: pd.DataFrame, *, full: bool) -> dict[str, Any]:
    """Validate counts, views, file 49, and forbidden source routes."""

    result = {
        "trial_count": len(metadata),
        "strict_clean_count": int(metadata["strict_clean_eligibility"].sum()),
        "duration_warning_count": int(metadata["duration_warning_flag"].sum()),
        "file49_count": int(metadata["eeg_source_id"].eq(49).sum()),
        "file54_1_count": int(
            metadata["source_recording_filename"].str.contains("demi_54_1", regex=False).sum()
        ),
        "id86_count": int(metadata["eeg_source_id"].eq(86).sum()),
    }
    if full and result != {
        "trial_count": EXPECTED_TRIAL_COUNT,
        "strict_clean_count": EXPECTED_STRICT_CLEAN_COUNT,
        "duration_warning_count": EXPECTED_DURATION_WARNING_COUNT,
        "file49_count": EXPECTED_FILE49_COUNT,
        "file54_1_count": 0,
        "id86_count": 0,
    }:
        raise RuntimeError(f"complete TFR surface count mismatch: {result}")
    return result


def write_axes(output_root: Path, support: Mapping[str, Any]) -> dict[str, Any]:
    """Write deterministic channel/frequency/time/baseline axis tables."""

    axes_dir = output_root / "axes"
    channel_frame = pd.DataFrame(
        {
            "channel_index": np.arange(len(PRIMARY_SCALP_CHANNELS), dtype=int),
            "channel": PRIMARY_SCALP_CHANNELS,
            "type": "eeg",
            "primary_tfr": True,
        }
    )
    frequency_frame = pd.DataFrame(
        {
            "frequency_index": np.arange(len(FREQUENCIES_HZ), dtype=int),
            "frequency_hz": FREQUENCIES_HZ,
            "n_cycles": N_CYCLES,
            "planned_scientific_surface": FREQUENCIES_HZ <= 30.0,
            "exploratory_diagnostic": FREQUENCIES_HZ >= 31.0,
            "wavelet_samples_at_100hz": support["wavelet_samples"],
            "half_support_seconds": support["half_support_seconds"],
        }
    )
    time_frame = pd.DataFrame(
        {
            "time_index": np.arange(len(EXPECTED_TASK_TIMES), dtype=int),
            "time_seconds": EXPECTED_TASK_TIMES,
        }
    )
    baseline_times = np.arange(-50, -19, dtype=np.float64) / 100.0
    baseline_frame = pd.DataFrame(
        {
            "baseline_time_index": np.arange(len(baseline_times), dtype=int),
            "time_seconds": baseline_times,
            "endpoint_inclusive": True,
        }
    )
    descriptors = {}
    for name, frame in (
        ("channels", channel_frame),
        ("frequencies_cycles", frequency_frame),
        ("task_times", time_frame),
        ("baseline_times", baseline_frame),
    ):
        path = axes_dir / f"{name}.csv"
        atomic_write_csv(path, frame)
        descriptors[name] = descriptor_for_file(path, root=output_root, row_count=len(frame))
    return descriptors


def storage_summary(output_root: Path) -> dict[str, Any]:
    """Summarize current production arrays and all output files."""

    current_arrays = [
        path
        for path in (output_root / "recordings").rglob("*.npy")
        if not any(part.startswith(".") for part in path.relative_to(output_root).parts)
    ]
    all_files = [path for path in output_root.rglob("*") if path.is_file()]
    return {
        "current_array_count": len(current_arrays),
        "current_array_bytes": sum(path.stat().st_size for path in current_arrays),
        "all_file_count": len(all_files),
        "all_file_bytes": sum(path.stat().st_size for path in all_files),
    }


def build_markdown_summary(summary: Mapping[str, Any]) -> str:
    """Render a concise local completion/partial summary."""

    storage_gib = summary["storage"]["current_array_bytes"] / 1024**3
    status_word = "complete" if summary["complete_surface"] else "bounded partial"
    return f"""# DEMI EEG trial-level TFR summary

Status: **{status_word}, reopened, and validated**  
Completed: {summary['completed_at']}

The current namespace contains {summary['trial_count']:,} row-aligned trials
from {summary['recording_count']} recording(s), with 30 primary scalp channels,
37 integer frequencies from 4--40 Hz, and 201 retained task samples from
-0.5 through +1.5 seconds. Frequencies 4--30 Hz are the planned scientific
surface; 31--40 Hz are retained as exploratory/diagnostic.

Response-onset and response-end raw Morlet power and trial-matched dB power are
stored as float32 NumPy arrays. Computation and transforms used float64. The
shared `red_on` baseline is the inclusive -0.5 to -0.2-second mean power for
each trial, channel, and frequency. It is a trial-matched visual pre-response
reference, not neutral rest. Unnormalized log power is the deterministic view
`10 * log10(raw_power)` and is not duplicated as a third cube.

All {summary['baseline_invalid_count']} invalid baseline cells/rows were
accounted for. Transient artifact facts are diagnostic metadata only and did
not change eligibility. No interpolation, rejection, AutoReject, CSD, complex,
phase, connectivity, ROI, band-power, participant, model, or ID-86 output was
created. Current arrays occupy approximately {storage_gib:.2f} GiB.
"""


def write_aggregate_outputs(
    *,
    repo_root: Path,
    output_root: Path,
    authority: Mapping[str, Any],
    selected: list[Mapping[str, Any]],
    manifests: list[Mapping[str, Any]],
    source_validation: Mapping[str, Any],
    started_at: str,
    started_perf: float,
    benchmark: bool,
) -> dict[str, Any]:
    """Validate and publish axes, metadata views, family manifests, and summaries."""

    full = len(selected) == EXPECTED_RECORDING_COUNT
    support = expected_wavelet_support()
    axes = write_axes(output_root, support)
    config_snapshot_path = output_root / "configuration" / "tfr_config_v1.yaml"
    atomic_write_text(
        config_snapshot_path,
        (repo_root / CONFIG_PATH).read_text(encoding="utf-8"),
    )
    config_snapshot = descriptor_for_file(config_snapshot_path, root=output_root)
    if config_snapshot["sha256"] != authority["config_sha256"]:
        raise RuntimeError("configuration snapshot hash differs")

    family_metadata_parts = {
        family: [
            pd.read_parquet(
                output_root
                / "recordings"
                / manifest["recording"]["recording_stem"]
                / manifest["tables"][f"{family}_metadata"]["relative_path"]
            )
            for manifest in manifests
        ]
        for family in TASK_FAMILIES
    }
    family_metadata = {
        family: pd.concat(parts, ignore_index=True) for family, parts in family_metadata_parts.items()
    }
    common = family_metadata["response_onset"]
    special = validate_special_surface(common, full=full)
    if family_metadata["response_end"]["canonical_event_key"].astype(str).tolist() != common[
        "canonical_event_key"
    ].astype(str).tolist():
        raise RuntimeError("aggregate task family canonical keys differ")

    metadata_dir = output_root / "metadata"
    metadata_descriptors = {}
    for family, frame in family_metadata.items():
        path = metadata_dir / f"{family}_row_metadata.parquet"
        metadata_descriptors[family] = write_frame_descriptor(path, frame, root=output_root)
    diagnostics = pd.concat(
        [
            pd.read_parquet(
                output_root
                / "recordings"
                / manifest["recording"]["recording_stem"]
                / manifest["tables"]["artifact_diagnostics"]["relative_path"]
            )
            for manifest in manifests
        ],
        ignore_index=True,
    )
    diagnostic_path = metadata_dir / "artifact_diagnostics.parquet"
    diagnostic_descriptor = write_frame_descriptor(
        diagnostic_path, diagnostics, root=output_root
    )
    diagnostic_summary = (
        diagnostics.groupby(["family", "condition_semantics"], sort=True)
        .agg(
            trial_count=("canonical_event_key", "size"),
            robust_transient_flag_count=("diagnostic_robust_logp2p_z_gt_6", "sum"),
            jump_flag_count=("diagnostic_jump_gt_50uv_per_native_sample", "sum"),
            flat_flag_count=("diagnostic_flat_p2p_lt_1uv", "sum"),
            scalp_p2p_max_p50_uv=("scalp_p2p_max_uv", "median"),
            scalp_p2p_max_uv=("scalp_p2p_max_uv", "max"),
            eog_p2p_max_p95_uv=("eog_p2p_max_uv", lambda x: float(np.quantile(x, 0.95))),
            emg_p2p_max_p95_uv=("emg_p2p_max_uv", lambda x: float(np.quantile(x, 0.95))),
        )
        .reset_index()
    )
    diagnostic_summary_path = metadata_dir / "artifact_diagnostic_summary.csv"
    atomic_write_csv(diagnostic_summary_path, diagnostic_summary)
    diagnostic_summary_descriptor = descriptor_for_file(
        diagnostic_summary_path, root=output_root, row_count=len(diagnostic_summary)
    )

    strict_view = common.loc[
        common["strict_clean_eligibility"],
        ["tfr_global_row_index", "canonical_order_index", "canonical_event_key"],
    ].reset_index(drop=True)
    warning_view = common.loc[
        common["duration_warning_flag"],
        [
            "tfr_row_index",
            "canonical_order_index",
            "canonical_event_key",
            "behavioural_id",
            "eeg_source_id",
            "source_recording_filename",
        ],
    ].reset_index(drop=True)
    strict_path = metadata_dir / "strict_clean_row_indices.parquet"
    warning_path = metadata_dir / "duration_warning_rows.parquet"
    strict_descriptor = write_frame_descriptor(strict_path, strict_view, root=output_root)
    warning_descriptor = write_frame_descriptor(warning_path, warning_view, root=output_root)

    family_manifest_descriptors = {}
    manifests_dir = output_root / "manifests"
    for family in TASK_FAMILIES:
        payload = {
            "schema_version": TFR_SCHEMA_VERSION,
            "tfr_namespace": TFR_NAMESPACE,
            "status": "complete" if full else "partial",
            "family": family,
            "trial_count": len(family_metadata[family]),
            "recording_count": len(manifests),
            "shape_tail": [30, 37, 201],
            "axis_order": ["trial", "channel", "frequency", "time"],
            "raw_power_dtype": "float32",
            "db_power_dtype": "float32",
            "metadata": metadata_descriptors[family],
            "shards": [
                {
                    "recording": manifest["recording"]["recording_stem"],
                    "row_count": manifest["recording"]["row_count"],
                    "raw": manifest["arrays"][f"{family}_raw"],
                    "db": manifest["arrays"][f"{family}_db"],
                    "fingerprint": manifest["fingerprint"],
                    "cache_action": manifest["invocation_cache_action"],
                }
                for manifest in manifests
            ],
        }
        path = manifests_dir / f"{family}_manifest.json"
        atomic_write_json(path, payload)
        family_manifest_descriptors[family] = descriptor_for_file(path, root=output_root)
    baseline_manifest = {
        "schema_version": TFR_SCHEMA_VERSION,
        "tfr_namespace": TFR_NAMESPACE,
        "status": "complete" if full else "partial",
        "family": "red_on",
        "product": "trial_matched_baseline_mean_power",
        "interval_seconds_inclusive": [BASELINE_TMIN, BASELINE_TMAX],
        "trial_count": len(common),
        "recording_count": len(manifests),
        "shape_tail": [30, 37],
        "axis_order": ["trial", "channel", "frequency"],
        "dtype": "float32",
        "invalid_count": 0,
        "shards": [
            {
                "recording": manifest["recording"]["recording_stem"],
                "row_count": manifest["recording"]["row_count"],
                "baseline_mean": manifest["arrays"]["red_on_baseline_mean"],
                "validation": manifest["validation"]["baseline"],
                "fingerprint": manifest["fingerprint"],
            }
            for manifest in manifests
        ],
    }
    baseline_manifest_path = manifests_dir / "red_on_baseline_manifest.json"
    atomic_write_json(baseline_manifest_path, baseline_manifest)
    family_manifest_descriptors["red_on"] = descriptor_for_file(
        baseline_manifest_path, root=output_root
    )

    forbidden = forbidden_output_scan(output_root)
    if forbidden:
        raise RuntimeError(f"unauthorized TFR products found: {forbidden}")
    storage = storage_summary(output_root)
    completed_at = utc_now()
    cache_counts = pd.Series(
        [manifest["invocation_cache_action"] for manifest in manifests]
    ).value_counts().to_dict()
    summary = {
        "schema_version": TFR_SCHEMA_VERSION,
        "tfr_namespace": TFR_NAMESPACE,
        "status": "complete" if full else "partial",
        "complete_surface": full,
        "completed_at": completed_at,
        "runtime_seconds": time.perf_counter() - started_perf,
        "recording_count": len(manifests),
        "trial_count": len(common),
        "strict_clean_count": special["strict_clean_count"],
        "duration_warning_count": special["duration_warning_count"],
        "file49_count": special["file49_count"],
        "file54_1_count": special["file54_1_count"],
        "id86_count": special["id86_count"],
        "channel_count": len(PRIMARY_SCALP_CHANNELS),
        "frequency_count": len(FREQUENCIES_HZ),
        "task_time_sample_count": len(EXPECTED_TASK_TIMES),
        "baseline_invalid_count": 0,
        "artifact_diagnostic_row_count": len(diagnostics),
        "cache_actions": cache_counts,
        "source_epoch_immutability": source_validation,
        "storage": storage,
        "forbidden_output_scan": forbidden,
    }
    atomic_write_json(output_root / SUMMARY_JSON_NAME, summary)
    atomic_write_text(output_root / SUMMARY_MD_NAME, build_markdown_summary(summary))
    validation = {
        "schema_version": TFR_SCHEMA_VERSION,
        "valid": True,
        "complete_surface": full,
        "special_surface": special,
        "source_immutability": source_validation,
        "all_recording_manifests_complete": all(
            manifest["status"] == "complete" for manifest in manifests
        ),
        "all_baselines_valid": all(
            manifest["validation"]["baseline"]["nonfinite_count"] == 0
            and manifest["validation"]["baseline"]["nonpositive_count"] == 0
            for manifest in manifests
        ),
        "all_formula_checks_valid": all(
            manifest["validation"][f"{family}_formula"][
                "maximum_absolute_db_difference"
            ]
            <= 2e-4
            for manifest in manifests
            for family in TASK_FAMILIES
        ),
        "log_power_formula": "10 * log10(raw_power)",
        "log_power_full_cube_stored": False,
        "strict_clean_is_view": True,
        "duration_warning_is_view": True,
        "unauthorized_products": forbidden,
    }
    validation_path = output_root / "validation" / "tfr_validation.json"
    atomic_write_json(validation_path, validation)
    run_seed = {
        "authority_fingerprint": authority["authority_fingerprint"],
        "recordings": [manifest["recording"]["recording_stem"] for manifest in manifests],
        "fingerprints": [manifest["fingerprint"] for manifest in manifests],
    }
    run_manifest = {
        "schema_version": TFR_SCHEMA_VERSION,
        "tfr_namespace": TFR_NAMESPACE,
        "run_id": f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')}__{content_fingerprint(run_seed)[:10]}",
        "status": "complete" if full else "partial",
        "started_at": started_at,
        "completed_at": completed_at,
        "runtime_seconds": time.perf_counter() - started_perf,
        "authority_fingerprint": authority["authority_fingerprint"],
        "accepted_epoch_run_manifest": {
            "path": authority["run_manifest_path"].relative_to(repo_root).as_posix(),
            "sha256": sha256_file(authority["run_manifest_path"]),
        },
        "accepted_family_manifests": {
            family: {
                "path": authority["family_manifest_paths"][family]
                .relative_to(repo_root)
                .as_posix(),
                "sha256": sha256_file(authority["family_manifest_paths"][family]),
            }
            for family in ALL_FAMILIES
        },
        "code_files": authority["code_files"],
        "code_sha256": authority["code_sha256"],
        "config": config_snapshot,
        "software": authority["software"],
        "axes": axes,
        "family_manifests": family_manifest_descriptors,
        "strict_clean_view": strict_descriptor,
        "duration_warning_view": warning_descriptor,
        "artifact_diagnostics": diagnostic_descriptor,
        "artifact_diagnostic_summary": diagnostic_summary_descriptor,
        "source_immutability_path": SOURCE_IMMUTABILITY_NAME,
        "validation_path": validation_path.relative_to(output_root).as_posix(),
        "summary_json_path": SUMMARY_JSON_NAME,
        "summary_markdown_path": SUMMARY_MD_NAME,
        "storage": storage,
        "cache_actions": cache_counts,
        "recordings": [
            {
                "recording": manifest["recording"]["recording_stem"],
                "row_count": manifest["recording"]["row_count"],
                "fingerprint": manifest["fingerprint"],
                "cache_action": manifest["invocation_cache_action"],
                "manifest_path": (
                    Path("recordings")
                    / manifest["recording"]["recording_stem"]
                    / "manifest.json"
                ).as_posix(),
            }
            for manifest in manifests
        ],
        "nonoperations": {
            "epoch_interpolation": False,
            "epoch_rejection": False,
            "autoreject": False,
            "csd": False,
            "complex": False,
            "phase": False,
            "connectivity": False,
            "roi": False,
            "band_power": False,
            "participant_inclusion": False,
            "models": False,
        },
    }
    if benchmark:
        benchmark_manifest = manifests[0]
        shard_bytes = sum(
            int(item["size_bytes"]) for item in benchmark_manifest["arrays"].values()
        )
        projection_factor = EXPECTED_TRIAL_COUNT / benchmark_manifest["recording"]["row_count"]
        benchmark_payload = {
            "schema_version": TFR_SCHEMA_VERSION,
            "recording": benchmark_manifest["recording"],
            "runtime_seconds": benchmark_manifest["runtime_seconds"],
            "peak_rss_bytes": benchmark_manifest["peak_rss_bytes_after"],
            "array_bytes": shard_bytes,
            "reopen_validation_complete": True,
            "projected_complete_runtime_seconds_linear_by_trial": (
                benchmark_manifest["runtime_seconds"] * projection_factor
            ),
            "projected_complete_array_bytes_linear_by_trial": int(
                shard_bytes * projection_factor
            ),
            "accepted_runtime_boundary_seconds": 3_600,
            "safe_to_launch_if_projection_and_free_space_pass": True,
        }
        benchmark_dir = output_root / "benchmarks"
        benchmark_path = benchmark_dir / f"{benchmark_manifest['recording']['recording_stem']}.json"
        atomic_write_json(benchmark_path, benchmark_payload)
        run_manifest["benchmark_path"] = benchmark_path.relative_to(output_root).as_posix()
    atomic_write_json(output_root / RUN_MANIFEST_NAME, run_manifest)
    return run_manifest


def run_stage(
    repo_root: Path,
    *,
    requested_recording: str | None,
    benchmark: bool,
    force: bool,
    verify_current: bool,
) -> dict[str, Any]:
    """Execute a bounded or complete stage-16 invocation."""

    started_perf = time.perf_counter()
    started_at = utc_now()
    output_root = (repo_root / OUTPUT_ROOT).resolve()
    require_ignored_output_root(repo_root, output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    authority = load_epoch_authority(repo_root)
    selected = select_recordings(authority, requested_recording)
    source_paths = [
        Path(recording["sources"][family]["path"])
        for recording in selected
        for family in ALL_FAMILIES
    ]
    source_before = source_snapshot(source_paths)
    run_state = {
        "schema_version": TFR_SCHEMA_VERSION,
        "tfr_namespace": TFR_NAMESPACE,
        "status": "running",
        "started_at": started_at,
        "selected_recording_count": len(selected),
        "selected_trial_count": sum(int(item["row_count"]) for item in selected),
        "completed_recordings": 0,
        "requested_recording": requested_recording,
        "benchmark": benchmark,
        "verify_current": verify_current,
    }
    atomic_write_json(output_root / RUN_STATE_NAME, run_state)
    print_progress(
        f"TFR stage start recordings={len(selected)} "
        f"trials={run_state['selected_trial_count']} benchmark={benchmark} "
        f"verify_current={verify_current}"
    )
    manifests = []
    source_hashes = []
    for index, recording in enumerate(selected, start=1):
        print_progress(
            f"recording {index}/{len(selected)} {recording['recording']} "
            f"rows={recording['row_count']}"
        )
        for family in ALL_FAMILIES:
            source_hashes.append(validate_source_file(recording["sources"][family]))
        metadata = read_family_metadata(recording)
        fingerprint = recording_cache_fingerprint(
            authority, recording, metadata
        )
        final_dir = output_root / "recordings" / recording["recording"]
        contracts = expected_array_contracts(recording["row_count"])
        assessment = assess_cached_recording(final_dir, fingerprint, contracts)
        if verify_current and assessment["action"] != "reuse":
            raise RuntimeError(
                f"current output is not reusable for {recording['recording']}: "
                f"{assessment['reason']}"
            )
        manifest = build_recording(
            repo_root=repo_root,
            output_root=output_root,
            authority=authority,
            recording=recording,
            metadata=metadata,
            fingerprint=fingerprint,
            force=force,
        )
        manifests.append(manifest)
        run_state["completed_recordings"] = index
        run_state["current_recording"] = recording["recording"]
        run_state["elapsed_seconds"] = time.perf_counter() - started_perf
        atomic_write_json(output_root / RUN_STATE_NAME, run_state)
        print_progress(
            f"recording complete {index}/{len(selected)} {recording['recording']} "
            f"cache={manifest['invocation_cache_action']} "
            f"elapsed={time.perf_counter() - started_perf:.1f}s"
        )

    source_after = source_snapshot(source_paths)
    comparison = compare_source_snapshots(source_before, source_after)
    source_validation = {
        **comparison,
        "checked_at": utc_now(),
        "pre_snapshot": source_before,
        "post_snapshot": source_after,
        "files": source_hashes,
    }
    atomic_write_json(output_root / SOURCE_IMMUTABILITY_NAME, source_validation)
    run_manifest = write_aggregate_outputs(
        repo_root=repo_root,
        output_root=output_root,
        authority=authority,
        selected=selected,
        manifests=manifests,
        source_validation=source_validation,
        started_at=started_at,
        started_perf=started_perf,
        benchmark=benchmark,
    )
    run_state.update(
        {
            "status": run_manifest["status"],
            "completed_at": utc_now(),
            "elapsed_seconds": time.perf_counter() - started_perf,
        }
    )
    run_state.pop("current_recording", None)
    atomic_write_json(output_root / RUN_STATE_NAME, run_state)
    print_progress(
        f"PASS TFR stage status={run_manifest['status']} "
        f"recordings={len(selected)} runtime={time.perf_counter() - started_perf:.1f}s"
    )
    return run_manifest


def main() -> None:
    """Run stage 16 and publish a readable failed run-state reason."""

    args = parse_args()
    root = repo_root_from_script()
    os.chdir(root)
    try:
        run_stage(
            root,
            requested_recording=args.recording,
            benchmark=args.benchmark,
            force=args.force,
            verify_current=args.verify_current,
        )
    except Exception as error:  # noqa: BLE001 - terminal driver must surface exact contract.
        print_progress(f"FAIL TFR stage: {type(error).__name__}: {error}")
        output_root = root / OUTPUT_ROOT
        if output_root.is_dir():
            atomic_write_json(
                output_root / RUN_STATE_NAME,
                {
                    "schema_version": TFR_SCHEMA_VERSION,
                    "tfr_namespace": TFR_NAMESPACE,
                    "status": "failed",
                    "failed_at": utc_now(),
                    "reason_code": type(error).__name__,
                    "readable_reason": str(error),
                },
            )
        raise


if __name__ == "__main__":
    main()
