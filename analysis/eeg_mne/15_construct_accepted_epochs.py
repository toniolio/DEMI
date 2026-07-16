"""Construct the accepted DEMI response and ``red_on`` EEG epoch families.

This is the active epoch-construction stage for the Python/MNE reanalysis.  It
exists to turn the accepted event/epoch eligibility ledger into three separate,
independently traceable voltage-epoch families without making later artifact,
spectral, channel, participant, ROI, or statistical decisions.

Inputs:
    ``_Data/eeg/event_epoch_eligibility_v1/`` canonical parquet ledger, its
    run manifest and summary, and only the linked authoritative
    ``continuous_v2`` post-ICA FIF/per-recording manifests.  Anchor timing and
    raw event provenance come exclusively from ledger fields; FIF annotations
    are never used to rediscover events.

Outputs:
    Bounded per-recording MNE Epochs FIF shards for ``response_onset``,
    ``response_end``, and ``red_on`` below ``_Data/eeg/epochs_v1/``; canonical
    family metadata; strict-clean and duration-warning row views; per-shard,
    per-family, and aggregate manifests; machine-readable validation and source
    immutability evidence; count tables; and a concise Markdown summary.

This script explicitly does not alter source FIFs, resample, apply a voltage
baseline, reject amplitudes or annotations, run AutoReject, apply CSD, compute
TFR/band power/spectral normalization, exclude participants, or construct an
ID-86 pre-ICA sensitivity.  Existing terminal shards are reused only when
their complete provenance fingerprint, size, hash, reopen validation, ordered
keys, and metadata match.  Incomplete or stale generated shards are preserved
under local history before rebuilding.

Safety boundaries:
    The output root must be the ignored versioned ``epochs_v1`` namespace.
    Every selected row and family must construct successfully.  Any missing
    anchor, half-sample ambiguity, source mismatch, MNE drop, out-of-bounds
    interval, count difference, or artifact validation failure stops the run.
    Source derivative size, modification time, and manifest checksum are
    verified unchanged before final acceptance.

Run from the repository root with immediate output flushing:

    PYTHONUNBUFFERED=1 MPLCONFIGDIR=/tmp/demi-mpl-cache \
      PATH="$(pwd)/.venv/bin:$PATH" \
      python3 analysis/eeg_mne/15_construct_accepted_epochs.py

The tracked ``tools/run_epochs_v1.sh`` launcher supplies these settings,
checks storage/ignore boundaries, keeps macOS awake when possible, and records
the terminal log.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time
from typing import Any, TypeVar

import mne
import numpy as np
import pandas as pd

from epoch_construction import (
    EPOCH_NAMESPACE,
    EPOCH_SCHEMA_VERSION,
    EXPECTED_DURATION_WARNING_COUNT,
    EXPECTED_EPOCH_COUNT,
    EXPECTED_FILE49_COUNT,
    EXPECTED_SFREQ,
    EXPECTED_STRICT_CLEAN_COUNT,
    FAMILIES,
    EpochFamily,
    assess_cached_epoch_shard,
    assert_identical_family_keys,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_parquet,
    atomic_write_text,
    build_family_metadata,
    compare_source_snapshots,
    content_fingerprint,
    create_epochs,
    forbidden_product_scan,
    ordered_key_hash,
    reopen_and_validate_epochs,
    select_epoch_rows,
    sha256_file,
    source_snapshot,
    stable_frame_hash,
    write_epochs_atomic,
)


LEDGER_DIR = Path("_Data/eeg/event_epoch_eligibility_v1")
LEDGER_PATH = LEDGER_DIR / "event_epoch_eligibility_ledger.parquet"
LEDGER_MANIFEST_PATH = LEDGER_DIR / "event_epoch_eligibility_run_manifest.json"
LEDGER_SUMMARY_PATH = LEDGER_DIR / "event_epoch_eligibility_summary.json"
CONTINUOUS_ROOT = Path("_Data/eeg/mne_preprocessing/continuous_v2")
OUTPUT_ROOT = Path("_Data/eeg/epochs_v1")
RUN_MANIFEST_NAME = "epoch_construction_run_manifest.json"
RUN_STATE_NAME = "epoch_construction_run_state.json"
SUMMARY_JSON_NAME = "epoch_construction_summary.json"
SUMMARY_MD_NAME = "epoch_construction_summary.md"
COUNT_SUMMARY_NAME = "epoch_count_summary.csv"
SOURCE_VALIDATION_NAME = "source_continuous_immutability.json"

T = TypeVar("T")


def utc_now() -> str:
    """Return a timezone-aware UTC timestamp with microsecond precision.

    Returns:
        ISO-8601 timestamp ending in ``Z``.

    Side effects:
        Reads the system clock.
    """

    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace(
        "+00:00", "Z"
    )


def repo_root_from_script() -> Path:
    """Return and validate the repository root containing this stage.

    Returns:
        Absolute repository root.

    Raises:
        RuntimeError: If expected project files are absent.

    Side effects:
        Resolves filesystem paths only.
    """

    root = Path(__file__).resolve().parents[2]
    if not (root / "analysis/eeg_mne/14_build_event_epoch_eligibility_ledger.py").is_file():
        raise RuntimeError(f"could not identify DEMI repository root from {__file__}")
    return root


def parse_args() -> argparse.Namespace:
    """Parse the intentionally small stage-15 command-line surface.

    Returns:
        Namespace containing the explicit ``--force`` flag.

    Side effects:
        Reads ``sys.argv`` and may print argparse help.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Preserve every existing generated shard in local history and rebuild "
            "the complete accepted surface. Ordinary reruns validate and reuse."
        ),
    )
    return parser.parse_args()


def print_progress(message: str) -> None:
    """Print one timestamped progress line with immediate flushing.

    Args:
        message: Human-readable progress detail.

    Side effects:
        Writes to standard output.
    """

    print(f"[{utc_now()}] {message}", flush=True)


def with_heartbeat(label: str, operation: Callable[[], T], interval: float = 30.0) -> T:
    """Run a blocking operation while emitting periodic elapsed-time updates.

    Args:
        label: Current recording/family operation shown to the user.
        operation: Zero-argument blocking callable.
        interval: Heartbeat interval in seconds.

    Returns:
        The callable's return value.

    Side effects:
        Starts one daemon thread, reads the clock, and prints heartbeats while
        ``operation`` is active.
    """

    stopped = threading.Event()
    started = time.perf_counter()

    def heartbeat() -> None:
        """Print elapsed time until the blocking operation completes."""

        while not stopped.wait(interval):
            print_progress(f"heartbeat {label}; elapsed={time.perf_counter() - started:.1f}s")

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        return operation()
    finally:
        stopped.set()
        thread.join(timeout=1.0)


def require_ignored_output_root(repo_root: Path, output_root: Path) -> None:
    """Require the exact versioned local root and confirm git ignores it.

    Args:
        repo_root: Absolute repository root.
        output_root: Proposed absolute epoch namespace.

    Raises:
        RuntimeError: If the root differs or is not ignored.

    Side effects:
        Runs read-only ``git check-ignore``.
    """

    expected = (repo_root / OUTPUT_ROOT).resolve()
    if output_root.resolve() != expected:
        raise RuntimeError(f"epoch outputs must use exactly {expected}")
    relative = output_root.resolve().relative_to(repo_root).as_posix()
    result = subprocess.run(
        ["git", "check-ignore", "--quiet", relative],
        cwd=repo_root,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"epoch output root is not ignored by git: {relative}")


def require_authoritative_input_path(
    path: Path, *, repo_root: Path, continuous_root: Path
) -> Path:
    """Resolve one linked derivative and require it below ``continuous_v2``.

    Args:
        path: Repository-relative path carried by the ledger.
        repo_root: Absolute repository root.
        continuous_root: Absolute authoritative continuous namespace.

    Returns:
        Resolved existing FIF path.

    Raises:
        RuntimeError: If the path escapes, is absent, or is not a post-ICA FIF.

    Side effects:
        Resolves and stats the path.
    """

    resolved = (repo_root / path).resolve() if not path.is_absolute() else path.resolve()
    if continuous_root.resolve() not in resolved.parents:
        raise RuntimeError(f"linked derivative escapes continuous_v2: {path}")
    if not resolved.is_file() or resolved.name != "post_ica_raw.fif":
        raise RuntimeError(f"linked authoritative post-ICA FIF is absent: {resolved}")
    return resolved


def load_input_authority(repo_root: Path) -> dict[str, Any]:
    """Load, hash, and validate the canonical ledger authority.

    Args:
        repo_root: Absolute repository root.

    Returns:
        Ledger, selected rows, accepted hashes/manifests, and code provenance.

    Raises:
        RuntimeError: For a missing input or any ledger hash/count mismatch.

    Side effects:
        Reads parquet/JSON and hashes tracked stage-15 code files.
    """

    paths = {
        "ledger": repo_root / LEDGER_PATH,
        "ledger_manifest": repo_root / LEDGER_MANIFEST_PATH,
        "ledger_summary": repo_root / LEDGER_SUMMARY_PATH,
    }
    missing = [label for label, path in paths.items() if not path.is_file()]
    if missing:
        raise RuntimeError(f"required epoch input authorities are absent: {missing}")
    ledger = pd.read_parquet(paths["ledger"])
    ledger_manifest = json.loads(paths["ledger_manifest"].read_text(encoding="utf-8"))
    ledger_summary = json.loads(paths["ledger_summary"].read_text(encoding="utf-8"))
    observed_content_hash = stable_frame_hash(ledger)
    accepted_content_hashes = {
        str(ledger_manifest.get("canonical_ledger_content_sha256", "")),
        str(ledger_summary.get("canonical_ledger_content_sha256", "")),
    }
    if accepted_content_hashes != {observed_content_hash}:
        raise RuntimeError(
            f"canonical ledger content hash differs: observed={observed_content_hash}, "
            f"accepted={accepted_content_hashes}"
        )
    selected = select_epoch_rows(ledger)
    code_paths = [
        repo_root / "analysis/eeg_mne/epoch_construction.py",
        repo_root / "analysis/eeg_mne/15_construct_accepted_epochs.py",
    ]
    code_files = {
        path.relative_to(repo_root).as_posix(): sha256_file(path) for path in code_paths
    }
    return {
        "ledger": ledger,
        "selected": selected,
        "ledger_manifest": ledger_manifest,
        "ledger_summary": ledger_summary,
        "ledger_content_sha256": observed_content_hash,
        "input_files": {
            label: {
                "path": path.relative_to(repo_root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for label, path in paths.items()
        },
        "code_files": code_files,
        "code_sha256": content_fingerprint(code_files),
    }


def one_value(rows: pd.DataFrame, column: str) -> Any:
    """Return the unique nonmissing value of a recording-level ledger field.

    Args:
        rows: Selected ledger rows for one derivative.
        column: Field required to be constant.

    Returns:
        Unique scalar value.

    Raises:
        RuntimeError: If the field is not constant and nonmissing.

    Side effects:
        None.
    """

    values = rows[column].dropna().unique()
    if len(values) != 1:
        raise RuntimeError(f"recording rows do not have exactly one {column}: {values}")
    return values[0]


def validate_continuous_link(
    rows: pd.DataFrame,
    *,
    repo_root: Path,
    continuous_root: Path,
) -> dict[str, Any]:
    """Validate one selected recording's ledger-to-manifest-to-FIF chain.

    Args:
        rows: Selected canonical rows sharing one derivative path.
        repo_root: Absolute repository root.
        continuous_root: Absolute authoritative continuous-v2 root.

    Returns:
        Paths and accepted source/artifact identifiers required by a shard.

    Raises:
        RuntimeError: For any path, status, manifest ID, artifact, size, source,
            or ledger/manifest mismatch.

    Side effects:
        Reads the linked continuous JSON manifest and stats its post-ICA FIF.
    """

    derivative_relative = Path(str(one_value(rows, "continuous_derivative_path")))
    derivative_path = require_authoritative_input_path(
        derivative_relative,
        repo_root=repo_root,
        continuous_root=continuous_root,
    )
    manifest_relative = Path(str(one_value(rows, "continuous_v2_manifest_path")))
    manifest_path = (repo_root / manifest_relative).resolve()
    if manifest_path.parent != derivative_path.parent or not manifest_path.is_file():
        raise RuntimeError(f"continuous manifest does not own linked derivative: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise RuntimeError(f"linked continuous result is not complete: {manifest_path}")
    if str(manifest.get("manifest_id")) != str(one_value(rows, "continuous_v2_manifest_id")):
        raise RuntimeError(f"continuous manifest ID differs from ledger: {manifest_path}")
    artifacts = {
        str(item.get("kind")): item
        for item in manifest.get("artifacts", [])
        if isinstance(item, dict)
    }
    artifact = artifacts.get("post_ica_continuous")
    if artifact is None:
        raise RuntimeError(f"linked manifest has no post-ICA artifact: {manifest_path}")
    artifact_path = (manifest_path.parent / str(artifact.get("relative_path"))).resolve()
    if artifact_path != derivative_path:
        raise RuntimeError(f"ledger and manifest post-ICA paths differ: {derivative_path}")
    if int(artifact.get("size_bytes", -1)) != derivative_path.stat().st_size:
        raise RuntimeError(f"post-ICA source size differs from manifest: {derivative_path}")
    source = manifest.get("source", {})
    if str(source.get("source_filename")) != str(one_value(rows, "source_filename")):
        raise RuntimeError(f"continuous source filename differs from ledger: {manifest_path}")
    if str(source.get("source_sha256")) != str(one_value(rows, "continuous_source_sha256")):
        raise RuntimeError(f"continuous raw-source hash differs from ledger: {manifest_path}")
    return {
        "recording_stem": derivative_path.parent.name,
        "source_filename": str(one_value(rows, "source_filename")),
        "derivative_path": derivative_path,
        "derivative_relative_path": derivative_path.relative_to(repo_root).as_posix(),
        "derivative_sha256": str(artifact.get("sha256")),
        "derivative_size_bytes": int(artifact.get("size_bytes")),
        "continuous_manifest_path": manifest_path,
        "continuous_manifest_relative_path": manifest_path.relative_to(repo_root).as_posix(),
        "continuous_manifest_id": str(manifest.get("manifest_id")),
        "continuous_manifest_sha256": sha256_file(manifest_path),
        "source_raw_sha256": str(source.get("source_sha256")),
        "expected_sample_count": int(source.get("sample_count")),
        "expected_sfreq": float(source.get("sampling_frequency_hz")),
        "completion_class": str(manifest.get("completion_class")),
        "qc_warnings": manifest.get("qc_warnings", []),
    }


def archive_previous_generated_result(
    *,
    artifact_path: Path,
    manifest_path: Path,
    output_root: Path,
    family_name: str,
    recording_stem: str,
    reason: str,
) -> list[str]:
    """Preserve stale/incomplete generated files in local history.

    Args:
        artifact_path: Existing shard artifact path, if present.
        manifest_path: Existing shard manifest path, if present.
        output_root: Versioned epoch namespace.
        family_name: Epoch family owning the generated files.
        recording_stem: Stable continuous recording stem.
        reason: Cache assessment or explicit-force reason.

    Returns:
        Output-relative paths of files moved into history.

    Side effects:
        Creates a history directory and atomically moves existing generated
        epoch files.  Source continuous files are never in scope.
    """

    existing = [path for path in (artifact_path, manifest_path) if path.exists()]
    if not existing:
        return []
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    history = output_root / "history" / stamp / family_name / recording_stem
    history.mkdir(parents=True, exist_ok=False)
    moved: list[str] = []
    for path in existing:
        target = history / path.name
        os.replace(path, target)
        moved.append(target.relative_to(output_root).as_posix())
    atomic_write_json(
        history / "archive_reason.json",
        {"archived_at": utc_now(), "reason": reason, "moved": moved},
    )
    return moved


def validate_reopened_raw(raw: mne.io.BaseRaw, link: dict[str, Any]) -> None:
    """Validate the linked continuous Raw header before epoch construction.

    Args:
        raw: Reopened authoritative post-ICA Raw object.
        link: Validated continuous linkage record.

    Raises:
        RuntimeError: For sampling, sample-count, first-sample, filename, bad
            list, or data-loading contract differences.

    Side effects:
        Reads Raw header attributes only.
    """

    problems: list[str] = []
    if not np.isclose(raw.info["sfreq"], EXPECTED_SFREQ, atol=1e-9, rtol=0.0):
        problems.append(f"sfreq={raw.info['sfreq']}")
    if raw.n_times != link["expected_sample_count"]:
        problems.append(f"n_times={raw.n_times} expected={link['expected_sample_count']}")
    if raw.first_samp != 0:
        problems.append(f"first_samp={raw.first_samp}")
    observed_filenames = tuple(
        Path(path).resolve() for path in raw.filenames if path is not None
    )
    if observed_filenames != (Path(link["derivative_path"]).resolve(),):
        problems.append(f"filenames={raw.filenames}")
    if raw.info["bads"]:
        problems.append(f"operational_bads={raw.info['bads']}")
    if raw.preload:
        problems.append("raw_unexpectedly_preloaded")
    if problems:
        raise RuntimeError(f"continuous Raw contract differs for {link['source_filename']}: {problems}")


def shard_fingerprint(
    *,
    authority: dict[str, Any],
    link: dict[str, Any],
    family: EpochFamily,
    metadata: pd.DataFrame,
) -> str:
    """Build the complete deterministic cache fingerprint for one shard.

    Args:
        authority: Canonical ledger and code provenance.
        link: Validated continuous derivative provenance.
        family: Accepted epoch family contract.
        metadata: Fully derived expected shard metadata.

    Returns:
        SHA-256 JSON content fingerprint.

    Side effects:
        None.
    """

    return content_fingerprint(
        {
            "epoch_schema_version": EPOCH_SCHEMA_VERSION,
            "epoch_namespace": EPOCH_NAMESPACE,
            "ledger_content_sha256": authority["ledger_content_sha256"],
            "code_sha256": authority["code_sha256"],
            "continuous_manifest_id": link["continuous_manifest_id"],
            "continuous_derivative_sha256": link["derivative_sha256"],
            "family": {
                "name": family.name,
                "event_id": family.event_id,
                "tmin": family.tmin,
                "tmax": family.tmax,
                "sfreq": EXPECTED_SFREQ,
                "baseline": None,
            },
            "ordered_key_sha256": ordered_key_hash(metadata["canonical_event_key"]),
            "metadata_content_sha256": stable_frame_hash(metadata),
        }
    )


def build_or_reuse_shard(
    *,
    raw: mne.io.BaseRaw,
    rows: pd.DataFrame,
    family: EpochFamily,
    authority: dict[str, Any],
    link: dict[str, Any],
    output_root: Path,
    force: bool,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Construct or validate/reuse one recording/family Epochs shard.

    Args:
        raw: Reopened authoritative continuous Raw.
        rows: Canonical selected rows for this recording.
        family: Accepted family contract.
        authority: Ledger/code provenance.
        link: Continuous derivative provenance.
        output_root: Versioned epoch namespace.
        force: Whether to preserve and rebuild matching terminal results.

    Returns:
        Expected metadata and a complete per-shard manifest.

    Side effects:
        May read/hash/reopen an existing epoch artifact, move stale generated
        files to history, preload accepted segments, write/reopen a new Epochs
        FIF, and atomically publish a terminal manifest.
    """

    metadata = build_family_metadata(
        rows,
        family,
        sfreq=float(raw.info["sfreq"]),
        first_samp=raw.first_samp,
        last_samp=raw.last_samp,
    )
    fingerprint = shard_fingerprint(
        authority=authority,
        link=link,
        family=family,
        metadata=metadata,
    )
    shard_dir = output_root / "families" / family.name / "recordings" / link["recording_stem"]
    artifact_path = shard_dir / f"{family.name}-epo.fif"
    manifest_path = shard_dir / "manifest.json"
    assessment = assess_cached_epoch_shard(manifest_path, artifact_path, fingerprint)

    if assessment["action"] == "reuse" and not force:
        cached = assessment["manifest"]
        validation = with_heartbeat(
            f"reopen cached {link['recording_stem']} {family.name}",
            lambda: reopen_and_validate_epochs(
                artifact_path,
                metadata,
                family,
                expected_sha256=cached["artifact"]["sha256"],
            ),
        )
        if cached.get("validation", {}).get("metadata_content_sha256") != validation.get(
            "metadata_content_sha256"
        ):
            raise RuntimeError(f"cached metadata validation differs: {artifact_path}")
        reused = dict(cached)
        reused["cache_action"] = "reused_after_hash_and_reopen_validation"
        reused["last_revalidated_at"] = utc_now()
        return metadata, reused

    reason = "explicit_force_rebuild" if force else assessment["reason"]
    archived = archive_previous_generated_result(
        artifact_path=artifact_path,
        manifest_path=manifest_path,
        output_root=output_root,
        family_name=family.name,
        recording_stem=link["recording_stem"],
        reason=reason,
    )
    print_progress(
        f"constructing {link['recording_stem']} family={family.name}; "
        f"epochs={len(metadata)}; cache_reason={reason}"
    )
    started = time.perf_counter()
    epochs = with_heartbeat(
        f"preload {link['recording_stem']} {family.name}",
        lambda: create_epochs(raw, metadata, family),
    )
    artifact = with_heartbeat(
        f"write {link['recording_stem']} {family.name}",
        lambda: write_epochs_atomic(epochs, artifact_path),
    )
    del epochs
    gc.collect()
    validation = with_heartbeat(
        f"reopen new {link['recording_stem']} {family.name}",
        lambda: reopen_and_validate_epochs(
            artifact_path,
            metadata,
            family,
            expected_sha256=artifact["sha256"],
        ),
    )
    artifact["relative_path"] = artifact_path.relative_to(output_root).as_posix()
    manifest = {
        "schema_version": EPOCH_SCHEMA_VERSION,
        "epoch_namespace": EPOCH_NAMESPACE,
        "status": "complete",
        "created_at": utc_now(),
        "cache_action": "constructed_and_reopened",
        "cache_reason": reason,
        "archived_previous_generated_files": archived,
        "fingerprint": fingerprint,
        "family": {
            "name": family.name,
            "event_id": family.event_id,
            "tmin_seconds": family.tmin,
            "tmax_seconds": family.tmax,
            "scientific_role": family.scientific_role,
            "endpoint_inclusive": True,
            "baseline": None,
        },
        "recording": {
            "recording_stem": link["recording_stem"],
            "source_filename": link["source_filename"],
            "continuous_derivative_path": link["derivative_relative_path"],
            "continuous_derivative_sha256": link["derivative_sha256"],
            "continuous_manifest_path": link["continuous_manifest_relative_path"],
            "continuous_manifest_id": link["continuous_manifest_id"],
        },
        "row_count": len(metadata),
        "strict_clean_count": int(metadata["strict_clean_eligibility"].sum()),
        "duration_warning_count": int(metadata["duration_warning_flag"].sum()),
        "canonical_order_min": int(metadata["canonical_order_index"].min()),
        "canonical_order_max": int(metadata["canonical_order_index"].max()),
        "ordered_key_sha256": ordered_key_hash(metadata["canonical_event_key"]),
        "metadata_content_sha256": stable_frame_hash(metadata),
        "artifact": artifact,
        "validation": validation,
        "runtime_seconds": time.perf_counter() - started,
    }
    atomic_write_json(manifest_path, manifest)
    return metadata, manifest


def validate_global_family_metadata(
    selected: pd.DataFrame, family_metadata: dict[str, pd.DataFrame]
) -> str:
    """Validate count/order/special-case contracts across complete families.

    Args:
        selected: Exact canonical selected ledger rows.
        family_metadata: Concatenated metadata for all three families.

    Returns:
        Common ordered canonical-key digest.

    Raises:
        RuntimeError: For any count, order, warning, special-route, or family
            difference.

    Side effects:
        None.
    """

    common_hash = assert_identical_family_keys(family_metadata)
    selected_keys = selected["canonical_event_key"].astype(str).tolist()
    for family in FAMILIES:
        frame = family_metadata[family.name]
        if len(frame) != EXPECTED_EPOCH_COUNT:
            raise RuntimeError(f"{family.name} has {len(frame)} rows")
        if frame["canonical_event_key"].astype(str).tolist() != selected_keys:
            raise RuntimeError(f"{family.name} differs from canonical ledger ordering")
        if not np.array_equal(
            frame["canonical_order_index"].to_numpy(dtype=int),
            np.arange(EXPECTED_EPOCH_COUNT),
        ):
            raise RuntimeError(f"{family.name} canonical indices are not contiguous")
        if int(frame["strict_clean_eligibility"].sum()) != EXPECTED_STRICT_CLEAN_COUNT:
            raise RuntimeError(f"{family.name} strict-clean count differs")
        if int(frame["duration_warning_flag"].sum()) != EXPECTED_DURATION_WARNING_COUNT:
            raise RuntimeError(f"{family.name} duration-warning count differs")
        file49 = frame[frame["eeg_source_id"].eq(49)]
        if len(file49) != EXPECTED_FILE49_COUNT or not file49[
            "continuous_qc_warning"
        ].str.contains("global_bad_proportion_above_25_percent", regex=False).all():
            raise RuntimeError(f"{family.name} file 49 warning/count differs")
        if frame["eeg_source_id"].eq(86).any():
            raise RuntimeError(f"{family.name} contains ID 86")
        if frame["source_recording_filename"].str.contains("demi_54_1", regex=False).any():
            raise RuntimeError(f"{family.name} contains file 54_1")
    return common_hash


def verify_source_derivatives(
    links: list[dict[str, Any]], before: list[dict[str, Any]]
) -> dict[str, Any]:
    """Hash every used source FIF and prove its size/mtime surface unchanged.

    Args:
        links: Validated recording linkage records in construction order.
        before: Pre-construction source size/mtime snapshot.

    Returns:
        Machine-readable immutability and checksum evidence.

    Raises:
        RuntimeError: If a source changed or differs from its accepted manifest.

    Side effects:
        Stats and sequentially hashes every used post-ICA source FIF while
        printing recording-level progress.
    """

    after = source_snapshot([link["derivative_path"] for link in links])
    comparison = compare_source_snapshots(before, after)
    hashes = []
    for index, link in enumerate(links, start=1):
        print_progress(
            f"source immutability hash {index}/{len(links)} recording={link['recording_stem']}"
        )
        observed = with_heartbeat(
            f"hash source {link['recording_stem']}",
            lambda path=link["derivative_path"]: sha256_file(path),
        )
        if observed != link["derivative_sha256"]:
            raise RuntimeError(
                f"source derivative hash differs from accepted manifest: {link['derivative_path']}"
            )
        hashes.append(
            {
                "recording_stem": link["recording_stem"],
                "path": link["derivative_relative_path"],
                "size_bytes": link["derivative_size_bytes"],
                "sha256": observed,
                "matches_accepted_continuous_manifest": True,
            }
        )
    return {
        **comparison,
        "checked_at": utc_now(),
        "pre_snapshot": before,
        "post_snapshot": after,
        "files": hashes,
    }


def output_storage_summary(output_root: Path) -> dict[str, Any]:
    """Summarize generated epoch namespace storage by file type and family.

    Args:
        output_root: Versioned epoch namespace.

    Returns:
        Total bytes/count plus per-family current FIF bytes/count.

    Side effects:
        Walks and stats output files.
    """

    all_files = [path for path in output_root.rglob("*") if path.is_file()]
    family_rows = {}
    for family in FAMILIES:
        files = list(
            (output_root / "families" / family.name / "recordings").rglob(
                "*-epo.fif"
            )
        )
        family_rows[family.name] = {
            "artifact_count": len(files),
            "bytes": sum(path.stat().st_size for path in files),
        }
    return {
        "file_count": len(all_files),
        "total_bytes": sum(path.stat().st_size for path in all_files),
        "families": family_rows,
    }


def build_markdown_summary(summary: dict[str, Any]) -> str:
    """Render the concise local epoch-construction summary.

    Args:
        summary: Final machine-readable summary payload.

    Returns:
        Markdown text ending in a newline.

    Side effects:
        None.
    """

    family_lines = []
    for family in summary["families"]:
        family_lines.append(
            f"- `{family['family']}`: {family['epoch_count']:,} epochs; "
            f"{family['strict_clean_count']:,} strict-clean; "
            f"{family['duration_warning_count']} duration warnings; "
            f"{family['sample_count_per_epoch']:,} inclusive samples per epoch."
        )
    storage_gib = summary["storage"]["total_bytes"] / (1024**3)
    return f"""# Accepted DEMI EEG epoch construction summary

Status: **complete and reopened**  
Completed: {summary['completed_at']}

## Constructed surface

{chr(10).join(family_lines)}

All families retain native 1000 Hz sampling, use exact `tmin = -1.5 s`,
include both requested endpoints, and have `baseline=None`. No row was dropped.
The strict-clean view is metadata/indexing over the same saved signal surface.

## Exceptional routes

- File 49 contributes {summary['file49_epoch_count']} rows per family with its
  8/30 interpolation warning retained.
- File 54_1 contributes no rows under the accepted split-file event policy.
- ID 86 contributes no rows because no accepted post-ICA derivative exists.

## Validation boundary

Every Epochs shard was reopened with preload and checked for ordered canonical
keys, event samples, metadata, finite data, sampling, inclusive endpoints,
empty drop logs, and absent baseline/reject/flat/detrend settings. The
{summary['source_immutability']['file_count']} linked continuous FIFs retained
their size, modification time, and accepted manifest checksum. The current
namespace occupies approximately {storage_gib:.2f} GiB.

No AutoReject, amplitude rejection, CSD, resampling, TFR, band-power, spectral
normalization, ROI definition, participant exclusion, statistical modeling, or
ID-86 sensitivity product was created.
"""


def construct_epochs(repo_root: Path, *, force: bool = False) -> dict[str, Any]:
    """Execute the complete accepted three-family epoch construction.

    Args:
        repo_root: Absolute repository root.
        force: Preserve and rebuild all existing generated shards when true.

    Returns:
        Final aggregate run manifest.

    Raises:
        RuntimeError: If any accepted construction or validation contract fails.

    Side effects:
        Reads authoritative ledger/continuous inputs and writes only the ignored
        versioned ``epochs_v1`` namespace with visible progress and heartbeats.
    """

    started_wall = time.perf_counter()
    started_at = utc_now()
    output_root = (repo_root / OUTPUT_ROOT).resolve()
    continuous_root = (repo_root / CONTINUOUS_ROOT).resolve()
    require_ignored_output_root(repo_root, output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    authority = load_input_authority(repo_root)
    selected = authority["selected"]
    groups = list(selected.groupby("continuous_derivative_path", sort=False))
    links = [
        validate_continuous_link(
            rows,
            repo_root=repo_root,
            continuous_root=continuous_root,
        )
        for _, rows in groups
    ]
    if len(links) != 81:
        raise RuntimeError(f"selected epoch surface uses {len(links)} recordings, expected 81")
    source_before = source_snapshot([link["derivative_path"] for link in links])
    run_seed = {
        "epoch_namespace": EPOCH_NAMESPACE,
        "ledger_content_sha256": authority["ledger_content_sha256"],
        "code_sha256": authority["code_sha256"],
        "continuous_v2_run_id": str(one_value(selected, "continuous_v2_run_id")),
    }
    run_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')}__{content_fingerprint(run_seed)[:10]}"
    run_state = {
        "schema_version": EPOCH_SCHEMA_VERSION,
        "run_id": run_id,
        "status": "running",
        "started_at": started_at,
        "force": force,
        "file_count": len(groups),
        "completed_shards": 0,
        "cumulative_epochs": {family.name: 0 for family in FAMILIES},
    }
    atomic_write_json(output_root / RUN_STATE_NAME, run_state)
    print_progress(
        f"accepted epoch run {run_id}; files={len(groups)}; rows={len(selected)}; force={force}"
    )

    family_parts: dict[str, list[pd.DataFrame]] = {
        family.name: [] for family in FAMILIES
    }
    family_shards: dict[str, list[dict[str, Any]]] = {
        family.name: [] for family in FAMILIES
    }
    cumulative = {family.name: 0 for family in FAMILIES}

    for file_index, ((_, rows), link) in enumerate(zip(groups, links, strict=True), start=1):
        file_started = time.perf_counter()
        print_progress(
            f"opening recording={link['recording_stem']} source={link['source_filename']} "
            f"file={file_index}/{len(groups)} rows={len(rows)}"
        )
        raw = mne.io.read_raw_fif(link["derivative_path"], preload=False, verbose="ERROR")
        try:
            validate_reopened_raw(raw, link)
            for family in FAMILIES:
                print_progress(
                    f"current recording={link['recording_stem']} family={family.name} "
                    f"file={file_index}/{len(groups)} cumulative={cumulative[family.name]}"
                )
                metadata, manifest = build_or_reuse_shard(
                    raw=raw,
                    rows=rows,
                    family=family,
                    authority=authority,
                    link=link,
                    output_root=output_root,
                    force=force,
                )
                family_parts[family.name].append(metadata)
                family_shards[family.name].append(manifest)
                cumulative[family.name] += len(metadata)
                run_state["completed_shards"] += 1
                run_state["cumulative_epochs"] = cumulative.copy()
                run_state["current"] = {
                    "file_index": file_index,
                    "file_total": len(groups),
                    "recording": link["recording_stem"],
                    "family": family.name,
                }
                run_state["elapsed_seconds"] = time.perf_counter() - started_wall
                atomic_write_json(output_root / RUN_STATE_NAME, run_state)
                print_progress(
                    f"complete recording={link['recording_stem']} family={family.name} "
                    f"file={file_index}/{len(groups)} cumulative={cumulative[family.name]} "
                    f"elapsed={time.perf_counter() - started_wall:.1f}s"
                )
        finally:
            raw.close()
            del raw
            gc.collect()
        print_progress(
            f"recording complete {file_index}/{len(groups)} {link['recording_stem']}; "
            f"file_elapsed={time.perf_counter() - file_started:.1f}s; "
            f"run_elapsed={time.perf_counter() - started_wall:.1f}s"
        )

    family_metadata = {
        name: pd.concat(parts, ignore_index=True) for name, parts in family_parts.items()
    }
    common_key_hash = validate_global_family_metadata(selected, family_metadata)
    metadata_dir = output_root / "metadata"
    selected_path = metadata_dir / "selected_canonical_ledger_rows.parquet"
    atomic_write_parquet(selected_path, selected)
    selected_descriptor = {
        "path": selected_path.relative_to(output_root).as_posix(),
        "size_bytes": selected_path.stat().st_size,
        "sha256": sha256_file(selected_path),
        "row_count": len(selected),
        "content_sha256": stable_frame_hash(selected),
    }

    family_metadata_descriptors: dict[str, dict[str, Any]] = {}
    for family in FAMILIES:
        path = metadata_dir / f"{family.name}_metadata.parquet"
        atomic_write_parquet(path, family_metadata[family.name])
        family_metadata_descriptors[family.name] = {
            "path": path.relative_to(output_root).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "row_count": len(family_metadata[family.name]),
            "content_sha256": stable_frame_hash(family_metadata[family.name]),
        }

    reference_metadata = family_metadata[FAMILIES[0].name]
    strict_view = reference_metadata.loc[
        reference_metadata["strict_clean_eligibility"],
        ["canonical_order_index", "canonical_event_key"],
    ].reset_index(drop=True)
    warning_view = reference_metadata.loc[
        reference_metadata["duration_warning_flag"],
        [
            "canonical_order_index",
            "canonical_event_key",
            "behavioural_id",
            "eeg_source_id",
            "source_recording_filename",
            "event_policy_reason_code",
        ],
    ].reset_index(drop=True)
    strict_parquet = metadata_dir / "strict_clean_row_indices.parquet"
    strict_csv = metadata_dir / "strict_clean_row_indices.csv"
    warning_parquet = metadata_dir / "duration_warning_rows.parquet"
    warning_csv = metadata_dir / "duration_warning_rows.csv"
    atomic_write_parquet(strict_parquet, strict_view)
    atomic_write_csv(strict_csv, strict_view)
    atomic_write_parquet(warning_parquet, warning_view)
    atomic_write_csv(warning_csv, warning_view)

    source_validation = verify_source_derivatives(links, source_before)
    atomic_write_json(output_root / SOURCE_VALIDATION_NAME, source_validation)

    family_summaries = []
    manifests_dir = output_root / "manifests"
    validation_dir = output_root / "validation"
    for family in FAMILIES:
        metadata = family_metadata[family.name]
        shards = family_shards[family.name]
        validations = [shard["validation"] for shard in shards]
        validation_payload = {
            "schema_version": EPOCH_SCHEMA_VERSION,
            "family": family.name,
            "valid": all(item.get("valid") for item in validations),
            "shard_count": len(validations),
            "epoch_count": sum(int(item["epoch_count"]) for item in validations),
            "zero_silent_drops": all(int(item["drop_count"]) == 0 for item in validations),
            "ordered_key_sha256": common_key_hash,
            "shards": validations,
        }
        validation_path = validation_dir / f"{family.name}_validation.json"
        atomic_write_json(validation_path, validation_payload)
        family_manifest = {
            "schema_version": EPOCH_SCHEMA_VERSION,
            "epoch_namespace": EPOCH_NAMESPACE,
            "status": "complete",
            "family": {
                "name": family.name,
                "event_id": family.event_id,
                "tmin_seconds": family.tmin,
                "tmax_seconds": family.tmax,
                "endpoint_inclusive": True,
                "sample_count_per_epoch": int(metadata["epoch_expected_sample_count"].iloc[0]),
                "sampling_frequency_hz": EXPECTED_SFREQ,
                "baseline": None,
            },
            "epoch_count": len(metadata),
            "strict_clean_count": int(metadata["strict_clean_eligibility"].sum()),
            "duration_warning_count": int(metadata["duration_warning_flag"].sum()),
            "recording_shard_count": len(shards),
            "ordered_key_sha256": common_key_hash,
            "metadata": family_metadata_descriptors[family.name],
            "validation_path": validation_path.relative_to(output_root).as_posix(),
            "shards": [
                {
                    "recording": shard["recording"]["recording_stem"],
                    "source_filename": shard["recording"]["source_filename"],
                    "row_count": shard["row_count"],
                    "strict_clean_count": shard["strict_clean_count"],
                    "duration_warning_count": shard["duration_warning_count"],
                    "artifact": shard["artifact"],
                    "fingerprint": shard["fingerprint"],
                    "cache_action": shard["cache_action"],
                }
                for shard in shards
            ],
        }
        family_manifest_path = manifests_dir / f"{family.name}_manifest.json"
        atomic_write_json(family_manifest_path, family_manifest)
        family_summaries.append(
            {
                "family": family.name,
                "epoch_count": len(metadata),
                "strict_clean_count": int(metadata["strict_clean_eligibility"].sum()),
                "duration_warning_count": int(metadata["duration_warning_flag"].sum()),
                "sampling_frequency_hz": EXPECTED_SFREQ,
                "tmin_seconds": family.tmin,
                "tmax_seconds": family.tmax,
                "sample_count_per_epoch": int(metadata["epoch_expected_sample_count"].iloc[0]),
                "endpoint_inclusive": True,
                "baseline": None,
                "recording_shard_count": len(shards),
                "manifest_path": family_manifest_path.relative_to(output_root).as_posix(),
                "metadata_path": family_metadata_descriptors[family.name]["path"],
                "validation_path": validation_path.relative_to(output_root).as_posix(),
            }
        )

    count_summary = pd.DataFrame(family_summaries)
    atomic_write_csv(output_root / COUNT_SUMMARY_NAME, count_summary)
    forbidden = forbidden_product_scan(output_root)
    if forbidden:
        raise RuntimeError(f"unauthorized products found in epoch namespace: {forbidden}")
    storage = output_storage_summary(output_root)
    completed_at = utc_now()
    summary = {
        "schema_version": EPOCH_SCHEMA_VERSION,
        "epoch_namespace": EPOCH_NAMESPACE,
        "status": "complete",
        "completed_at": completed_at,
        "runtime_seconds": time.perf_counter() - started_wall,
        "ledger_content_sha256": authority["ledger_content_sha256"],
        "continuous_v2_run_id": str(one_value(selected, "continuous_v2_run_id")),
        "ordinary_selected_row_count": EXPECTED_EPOCH_COUNT,
        "strict_clean_row_count": EXPECTED_STRICT_CLEAN_COUNT,
        "duration_warning_row_count": EXPECTED_DURATION_WARNING_COUNT,
        "file49_epoch_count": EXPECTED_FILE49_COUNT,
        "file54_1_epoch_count": 0,
        "id86_epoch_count": 0,
        "common_ordered_key_sha256": common_key_hash,
        "families": family_summaries,
        "contract_nonoperations": {
            "voltage_baseline_applied": False,
            "amplitude_rejection_applied": False,
            "annotation_rejection_applied": False,
            "autoreject_applied": False,
            "csd_applied": False,
            "resampling_applied": False,
            "tfr_computed": False,
            "band_power_computed": False,
            "spectral_normalization_applied": False,
            "participant_exclusion_applied": False,
            "id86_sensitivity_created": False,
        },
        "source_immutability": {
            "unchanged": source_validation["unchanged"],
            "file_count": source_validation["file_count"],
            "total_bytes": source_validation["total_bytes"],
            "validation_path": SOURCE_VALIDATION_NAME,
        },
        "storage": storage,
        "forbidden_product_scan": forbidden,
    }
    atomic_write_json(output_root / SUMMARY_JSON_NAME, summary)
    atomic_write_text(output_root / SUMMARY_MD_NAME, build_markdown_summary(summary))
    run_manifest = {
        "schema_version": EPOCH_SCHEMA_VERSION,
        "epoch_namespace": EPOCH_NAMESPACE,
        "run_id": run_id,
        "status": "complete",
        "started_at": started_at,
        "completed_at": completed_at,
        "runtime_seconds": time.perf_counter() - started_wall,
        "force": force,
        "input_authorities": authority["input_files"],
        "ledger_content_sha256": authority["ledger_content_sha256"],
        "code_files": authority["code_files"],
        "code_sha256": authority["code_sha256"],
        "continuous_v2_run_id": str(one_value(selected, "continuous_v2_run_id")),
        "selected_ledger_metadata": selected_descriptor,
        "strict_clean_views": [
            strict_parquet.relative_to(output_root).as_posix(),
            strict_csv.relative_to(output_root).as_posix(),
        ],
        "duration_warning_views": [
            warning_parquet.relative_to(output_root).as_posix(),
            warning_csv.relative_to(output_root).as_posix(),
        ],
        "families": family_summaries,
        "source_immutability_path": SOURCE_VALIDATION_NAME,
        "summary_json_path": SUMMARY_JSON_NAME,
        "summary_markdown_path": SUMMARY_MD_NAME,
        "count_summary_path": COUNT_SUMMARY_NAME,
        "storage": storage,
    }
    atomic_write_json(output_root / RUN_MANIFEST_NAME, run_manifest)
    run_state.update(
        {
            "status": "complete",
            "completed_at": completed_at,
            "elapsed_seconds": time.perf_counter() - started_wall,
            "cumulative_epochs": cumulative,
            "completed_shards": len(groups) * len(FAMILIES),
        }
    )
    run_state.pop("current", None)
    atomic_write_json(output_root / RUN_STATE_NAME, run_state)
    print_progress(
        f"PASS accepted epochs: onset={cumulative['response_onset']} "
        f"end={cumulative['response_end']} red_on={cumulative['red_on']} "
        f"runtime={time.perf_counter() - started_wall:.1f}s storage={storage['total_bytes']} bytes"
    )
    return run_manifest


def main() -> None:
    """Run stage 15, print a concise failure, and return a nonzero exit code.

    Side effects:
        Parses CLI arguments, changes to the repository root, constructs local
        outputs, and prints terminal status.
    """

    args = parse_args()
    root = repo_root_from_script()
    os.chdir(root)
    try:
        construct_epochs(root, force=args.force)
    except Exception as error:  # noqa: BLE001 - terminal driver must report the exact failed contract.
        print_progress(f"FAIL accepted epoch construction: {type(error).__name__}: {error}")
        raise


if __name__ == "__main__":
    main()
