"""Atomic derivative storage, ledgers, validation, and cache provenance.

This module implements the saved-derivative contract for DEMI continuous
preprocessing. It writes only new files below one configured ignored,
versioned root, publishes each recording directory atomically, preserves
interrupted/stale/superseded directories in local history, and validates all
terminal artifacts before they may be skipped on a later run.

Signal-bearing outputs are MNE FIF ``Raw`` files and ICA objects. JSON evidence
contains stage status, hashes, runtime, and validation results. The module
prohibits EDF, epoch, AutoReject, and CSD-like output names and never writes
inside the raw source directory.
"""

from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import re
import resource
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import mne
import numpy as np

from .contracts import (
    ALL_RETAINED_CHANNELS,
    CHANNEL_TYPE_BY_NAME,
    EEG_TARGET_CHANNELS,
    PIPELINE_VERSION,
    STAGE_NAMES,
    canonical_json_bytes,
    sha256_bytes,
    sha256_file,
)


TERMINAL_CACHEABLE_STATUSES = {"complete", "stopped"}
BLOCKED_OUTPUT_TOKENS = ("epoch", "autoreject", "csd")
ENVIRONMENT_PACKAGES = (
    "mne",
    "numpy",
    "scipy",
    "pandas",
    "pyprep",
    "scikit-learn",
    "PyYAML",
)


def utc_now() -> str:
    """Return a timezone-explicit UTC timestamp for operational evidence."""

    return datetime.now(timezone.utc).isoformat()


def history_stamp() -> str:
    """Return a filesystem-safe UTC stamp used only for preserved history."""

    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")


def peak_rss_bytes() -> int:
    """Return the process peak resident set size in bytes where practical."""

    maximum = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # macOS reports bytes; Linux reports KiB. The project runs on macOS, while
    # this branch keeps public synthetic tests meaningful on Linux as well.
    return maximum if sys.platform == "darwin" else maximum * 1024


def recording_key(source_filename: str) -> str:
    """Create a stable lowercase result-directory key from an EDF filename."""

    stem = Path(source_filename).stem.lower()
    key = re.sub(r"[^a-z0-9]+", "_", stem).strip("_")
    if not key or key.startswith("."):
        raise ValueError(f"Cannot create a safe recording key from {source_filename!r}.")
    return key


def parse_recording_id(source_filename: str) -> int:
    """Parse the numeric EEG recording identifier without inferring identity."""

    match = re.match(r"^demi_(\d{1,3})(?:_|\s|\.|$)", source_filename, flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Unrecognized DEMI EDF filename: {source_filename}")
    return int(match.group(1))


def validate_output_root(repo_root: Path, output_root: Path, configured_relative: str) -> None:
    """Require one exact authorized ignored versioned output directory.

    Args:
        repo_root: Repository root.
        output_root: Proposed resolved output root.
        configured_relative: Exact tracked configuration value.

    Returns:
        None.

    Side effects:
        Resolves paths; writes nothing.
    """

    expected = (repo_root / configured_relative).resolve()
    if output_root.resolve() != expected:
        raise ValueError(f"Output root must be exactly {expected}; found {output_root.resolve()}.")
    authorized = {
        "_Data/eeg/mne_preprocessing/continuous_validation_v1",
        "_Data/eeg/mne_preprocessing/continuous_v1",
        "_Data/eeg/mne_preprocessing/continuous_v2",
    }
    if configured_relative not in authorized:
        raise ValueError(
            "Only the authorized continuous_validation_v1, continuous_v1, or continuous_v2 root is permitted."
        )


def validate_derivative_path(path: Path, output_root: Path, raw_root: Path) -> None:
    """Reject raw-path, EDF, epoch, AutoReject, CSD, and escaped output paths.

    Args:
        path: Proposed output file.
        output_root: Authorized versioned validation or production root.
        raw_root: Immutable source EDF root.

    Returns:
        None.

    Side effects:
        Resolves paths; writes nothing.
    """

    resolved = path.resolve()
    root = output_root.resolve()
    raw = raw_root.resolve()
    if resolved == raw or raw in resolved.parents:
        raise ValueError(f"Derivative path may not be inside the raw EDF directory: {path}")
    if resolved != root and root not in resolved.parents:
        raise ValueError(f"Derivative path escapes the authorized output root: {path}")
    if resolved.suffix.lower() == ".edf":
        raise ValueError("Writing repaired or derived EDF files is prohibited.")
    lowered = resolved.name.lower()
    if any(token in lowered for token in BLOCKED_OUTPUT_TOKENS):
        raise ValueError(f"Prohibited epoch/AutoReject/CSD-like output name: {path.name}")


def atomic_write_json(path: Path, value: Any, output_root: Path, raw_root: Path) -> None:
    """Write deterministic JSON through a sibling temporary file and rename.

    Args:
        path: Final JSON path.
        value: JSON-compatible value; NaN is rejected.
        output_root: Authorized root.
        raw_root: Immutable source root.

    Returns:
        None.

    Side effects:
        Creates one temporary JSON file and atomically renames it to ``path``.
    """

    validate_derivative_path(path, output_root, raw_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}.json")
    validate_derivative_path(temporary, output_root, raw_root)
    payload = canonical_json_bytes(value) + b"\n"
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_write_text(path: Path, value: str, output_root: Path, raw_root: Path) -> None:
    """Write UTF-8 text through a sibling temporary file and atomic rename.

    Args:
        path: Final text path.
        value: Complete UTF-8 text.
        output_root: Authorized output root.
        raw_root: Immutable source root.

    Returns:
        None.

    Side effects:
        Creates and fsyncs one temporary file, then renames it to ``path``.
    """

    validate_derivative_path(path, output_root, raw_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}{path.suffix}")
    validate_derivative_path(temporary, output_root, raw_root)
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def artifact_descriptor(path: Path, result_dir: Path) -> dict[str, Any]:
    """Return a relative-path, byte-size, and SHA-256 artifact descriptor."""

    return {
        "relative_path": path.relative_to(result_dir).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def write_raw_fif_atomic(
    raw: mne.io.BaseRaw,
    path: Path,
    *,
    fif_format: str,
    output_root: Path,
    raw_root: Path,
) -> dict[str, Any]:
    """Atomically save one continuous MNE Raw FIF derivative.

    Args:
        raw: Continuous recording to save.
        path: Final path, which must end in ``_raw.fif``.
        fif_format: Explicit MNE precision (``single`` or ``double``).
        output_root: Authorized root.
        raw_root: Immutable source root.

    Returns:
        Size, checksum, precision, sample count, and sampling rate evidence.

    Side effects:
        Writes a temporary FIF and atomically renames it.
    """

    validate_derivative_path(path, output_root, raw_root)
    if not path.name.endswith("_raw.fif"):
        raise ValueError("Continuous FIF derivatives must end in _raw.fif.")
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite an existing derivative: {path}")
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}_raw.fif")
    validate_derivative_path(temporary, output_root, raw_root)
    raw.save(temporary, fmt=fif_format, overwrite=False, split_size="2GB", verbose="ERROR")
    os.replace(temporary, path)
    return {
        "path": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "fif_format": fif_format,
        "sample_count": int(raw.n_times),
        "sampling_frequency_hz": float(raw.info["sfreq"]),
        "channel_count": len(raw.ch_names),
    }


def write_ica_atomic(
    ica: mne.preprocessing.ICA,
    path: Path,
    *,
    output_root: Path,
    raw_root: Path,
) -> dict[str, Any]:
    """Atomically save one MNE ICA solution.

    Args:
        ica: Fitted ICA object.
        path: Final path ending in ``-ica.fif``.
        output_root: Authorized root.
        raw_root: Immutable source root.

    Returns:
        Size, checksum, component count, and recorded exclusion evidence.

    Side effects:
        Writes a temporary ICA FIF and atomically renames it.
    """

    validate_derivative_path(path, output_root, raw_root)
    if not path.name.endswith("-ica.fif"):
        raise ValueError("ICA objects must end in -ica.fif.")
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite an existing ICA object: {path}")
    temporary = path.with_name(f".{path.stem}.tmp-{uuid.uuid4().hex}-ica.fif")
    validate_derivative_path(temporary, output_root, raw_root)
    ica.save(temporary, overwrite=False, verbose="ERROR")
    os.replace(temporary, path)
    return {
        "path": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "n_components": int(ica.n_components_),
        "exclude": [int(index) for index in ica.exclude],
    }


@dataclass
class StageLedger:
    """Flush-safe state machine for the 20 accepted production stages."""

    path: Path
    output_root: Path
    raw_root: Path
    stages: Sequence[str] = STAGE_NAMES
    overall_status: str = "initialized"
    records: list[dict[str, Any]] = field(init=False)
    _active_stage: str | None = field(default=None, init=False)
    _active_started: float | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        """Initialize pending stages and flush the first ledger snapshot."""

        if tuple(self.stages) != STAGE_NAMES:
            raise ValueError("Stage ledger must expose the accepted 20-stage sequence.")
        self.records = [
            {"stage_index": index, "stage": stage, "status": "pending"}
            for index, stage in enumerate(self.stages, start=1)
        ]
        self.flush()

    def _record(self, stage: str) -> dict[str, Any]:
        """Return one stage record or raise for an unknown name."""

        for record in self.records:
            if record["stage"] == stage:
                return record
        raise KeyError(f"Unknown production stage: {stage}")

    @property
    def active_stage(self) -> str | None:
        """Return the currently running stage, if any."""

        return self._active_stage

    def start(self, stage: str) -> None:
        """Start the next pending stage and flush before expensive work."""

        if self._active_stage is not None:
            raise RuntimeError(f"Stage {self._active_stage} is already running.")
        record = self._record(stage)
        pending = next((row for row in self.records if row["status"] == "pending"), None)
        if record["status"] != "pending" or pending is not record:
            raise RuntimeError(f"Stage transition is out of order: {stage}")
        self._active_stage = stage
        self._active_started = time.perf_counter()
        self.overall_status = "running"
        record.update({"status": "running", "started_at": utc_now()})
        self.flush()

    def complete(self, stage: str, evidence: Mapping[str, Any] | None = None) -> None:
        """Complete the active stage with elapsed time, memory, and evidence."""

        if self._active_stage != stage or self._active_started is None:
            raise RuntimeError(f"Cannot complete inactive stage {stage}.")
        record = self._record(stage)
        record.update(
            {
                "status": "complete",
                "finished_at": utc_now(),
                "elapsed_seconds": time.perf_counter() - self._active_started,
                "process_peak_rss_bytes": peak_rss_bytes(),
                "evidence": dict(evidence or {}),
            }
        )
        self._active_stage = None
        self._active_started = None
        self.flush()

    def stop(self, stage: str, code: str, detail: str, evidence: Mapping[str, Any] | None = None) -> None:
        """Mark an objective stop and all later stages as not reached."""

        self._terminal(stage, "stopped", code, detail, evidence)

    def fail(self, stage: str, code: str, detail: str, evidence: Mapping[str, Any] | None = None) -> None:
        """Mark an implementation failure and all later stages as not reached."""

        self._terminal(stage, "failed", code, detail, evidence)

    def _terminal(
        self,
        stage: str,
        status: str,
        code: str,
        detail: str,
        evidence: Mapping[str, Any] | None,
    ) -> None:
        """Implement a validated terminal transition."""

        if self._active_stage != stage or self._active_started is None:
            raise RuntimeError(f"Cannot mark inactive stage {stage} as {status}.")
        record = self._record(stage)
        record.update(
            {
                "status": status,
                "finished_at": utc_now(),
                "elapsed_seconds": time.perf_counter() - self._active_started,
                "process_peak_rss_bytes": peak_rss_bytes(),
                "terminal_code": code,
                "terminal_detail": detail,
                "evidence": dict(evidence or {}),
            }
        )
        reached_terminal = False
        for later in self.records:
            if later is record:
                reached_terminal = True
                continue
            if reached_terminal and later["status"] == "pending":
                later["status"] = "not_reached"
                later["reason"] = f"recording_{status}_at_{stage}"
        self.overall_status = status
        self._active_stage = None
        self._active_started = None
        self.flush()

    def finalize(self, status: str) -> None:
        """Set a terminal overall status after all required stages settle."""

        if status not in {"complete", "stopped", "failed"}:
            raise ValueError(f"Invalid ledger terminal status: {status}")
        if self._active_stage is not None:
            raise RuntimeError("Cannot finalize a ledger with a running stage.")
        if status == "complete" and any(row["status"] != "complete" for row in self.records):
            raise RuntimeError("A complete ledger must have all 20 stages complete.")
        self.overall_status = status
        self.flush()

    def payload(self) -> dict[str, Any]:
        """Return the current JSON-compatible ledger payload."""

        return {
            "schema_version": 1,
            "pipeline_version": PIPELINE_VERSION,
            "overall_status": self.overall_status,
            "stage_count": len(self.records),
            "stages": self.records,
        }

    def flush(self) -> None:
        """Atomically persist the current ledger state."""

        atomic_write_json(self.path, self.payload(), self.output_root, self.raw_root)


def software_environment() -> dict[str, Any]:
    """Return exact Python/platform/package versions used for cache validity."""

    packages: dict[str, str | None] = {}
    for package in ENVIRONMENT_PACKAGES:
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    return {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": packages,
    }


def git_evidence(repo_root: Path) -> dict[str, Any]:
    """Record the current commit and tracked-worktree status without mutation."""

    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--short", "--untracked-files=no"],
        cwd=repo_root,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout.splitlines()
    return {"commit": commit, "tracked_worktree_changes": status, "tracked_worktree_clean": not status}


def code_provenance(repo_root: Path, code_paths: Sequence[Path]) -> dict[str, Any]:
    """Hash every production code file used to validate resumable outputs."""

    rows = []
    for path in sorted((path.resolve() for path in code_paths), key=lambda item: item.as_posix()):
        if not path.is_file():
            raise FileNotFoundError(f"Production code provenance path is missing: {path}")
        rows.append(
            {
                "path": path.relative_to(repo_root.resolve()).as_posix(),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    return {"files": rows, "sha256": sha256_bytes(canonical_json_bytes(rows))}


def run_provenance(
    repo_root: Path,
    config: Mapping[str, Any],
    code_paths: Sequence[Path],
    execution_context: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build config/code/environment hashes shared by every selected recording."""

    environment = software_environment()
    code = code_provenance(repo_root, code_paths)
    git = git_evidence(repo_root)
    execution = dict(execution_context or {})
    return {
        "pipeline_version": PIPELINE_VERSION,
        "config_path": config["_provenance"]["path"],
        "config_sha256": config["_provenance"]["sha256"],
        "code": code,
        "environment": environment,
        "environment_sha256": sha256_bytes(canonical_json_bytes(environment)),
        "execution": execution,
        "execution_sha256": sha256_bytes(canonical_json_bytes(execution)),
        "git": git,
    }


def recording_provenance(run: Mapping[str, Any], source_sha256: str) -> dict[str, Any]:
    """Create one stable source/config/code/environment cache fingerprint."""

    stable = {
        "pipeline_version": run["pipeline_version"],
        "source_sha256": source_sha256,
        "config_sha256": run["config_sha256"],
        "code_sha256": run["code"]["sha256"],
        "environment_sha256": run["environment_sha256"],
        "execution_sha256": run.get(
            "execution_sha256", sha256_bytes(canonical_json_bytes({}))
        ),
    }
    return {**stable, "fingerprint": sha256_bytes(canonical_json_bytes(stable))}


def deterministic_manifest_id(manifest: Mapping[str, Any]) -> str:
    """Hash stable scientific/provenance fields, excluding runtime timestamps."""

    stable = {
        "schema_version": manifest["schema_version"],
        "pipeline_version": manifest["pipeline_version"],
        "source": manifest["source"],
        "provenance": manifest["provenance"],
        "status": manifest["status"],
        "completion_class": manifest.get("completion_class"),
        "qc_warnings": manifest.get("qc_warnings", []),
        "stop_or_failure": manifest.get("stop_or_failure"),
        "channel_contract": manifest.get("channel_contract"),
        "detector": manifest.get("detector"),
        "reference": manifest.get("reference"),
        "interpolation": manifest.get("interpolation"),
        "ica": manifest.get("ica"),
    }
    return sha256_bytes(canonical_json_bytes(stable))


def validate_saved_raw(
    path: Path,
    *,
    expected_sample_count: int,
    expected_sfreq: float,
    expected_highpass: float,
    expected_lowpass: float,
    expected_sha256: str,
    chunk_seconds: float,
) -> dict[str, Any]:
    """Reopen and fully validate a saved continuous FIF derivative.

    Args:
        path: Saved FIF.
        expected_sample_count: Source-aligned sample count.
        expected_sfreq: Working sampling frequency.
        expected_highpass: Analysis high-pass metadata.
        expected_lowpass: Analysis low-pass metadata.
        expected_sha256: Checksum recorded after the atomic write.
        chunk_seconds: Bounded nonfinite scan chunk length.

    Returns:
        Reopen, channel, montage, sample, filter, reference, bad-list, finite,
        and checksum validation evidence.

    Side effects:
        Opens and scans the saved FIF read-only.
    """

    actual_hash = sha256_file(path)
    raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
    observed_types = dict(zip(raw.ch_names, raw.get_channel_types()))
    errors: list[str] = []
    if tuple(raw.ch_names) != ALL_RETAINED_CHANNELS:
        errors.append("channel_names")
    if observed_types != dict(CHANNEL_TYPE_BY_NAME):
        errors.append("channel_types")
    if raw.n_times != expected_sample_count:
        errors.append("sample_count")
    if not np.isclose(raw.info["sfreq"], expected_sfreq):
        errors.append("sampling_frequency")
    if not np.isclose(raw.info["highpass"], expected_highpass):
        errors.append("highpass_metadata")
    if not np.isclose(raw.info["lowpass"], expected_lowpass):
        errors.append("lowpass_metadata")
    if raw.info["bads"]:
        errors.append("operational_bad_list")
    if int(raw.info["custom_ref_applied"]) == 0:
        errors.append("reference_metadata")
    montage = raw.get_montage()
    montage_names = set(montage.ch_names) if montage else set()
    if not set(EEG_TARGET_CHANNELS).issubset(montage_names):
        errors.append("montage")
    if not all(channel in raw.ch_names for channel in ("M1", "M2")):
        errors.append("m1_m2_retention")
    step = max(1, int(round(chunk_seconds * raw.info["sfreq"])))
    nonfinite = 0
    chunks = 0
    for start in range(0, raw.n_times, step):
        data = raw.get_data(start=start, stop=min(raw.n_times, start + step))
        nonfinite += int(np.size(data) - np.isfinite(data).sum())
        chunks += 1
    if nonfinite:
        errors.append("nonfinite_samples")
    if actual_hash != expected_sha256:
        errors.append("output_checksum")
    return {
        "reopened_with_mne": True,
        "valid": not errors,
        "errors": errors,
        "channel_names_and_types_valid": "channel_names" not in errors and "channel_types" not in errors,
        "m1_m2_retained": "m1_m2_retention" not in errors,
        "montage_available": "montage" not in errors,
        "sample_count": int(raw.n_times),
        "sampling_frequency_hz": float(raw.info["sfreq"]),
        "highpass_hz": float(raw.info["highpass"]),
        "lowpass_hz": float(raw.info["lowpass"]),
        "custom_reference_applied": int(raw.info["custom_ref_applied"]) != 0,
        "operational_bad_list": list(raw.info["bads"]),
        "nonfinite_sample_count": nonfinite,
        "chunks_scanned": chunks,
        "sha256": actual_hash,
        "checksum_matches": actual_hash == expected_sha256,
    }


def validate_saved_ica(
    path: Path,
    *,
    expected_components: int,
    expected_exclusions: Sequence[int],
    expected_sha256: str,
) -> dict[str, Any]:
    """Reopen an ICA object and compare rank, exclusions, and checksum."""

    actual_hash = sha256_file(path)
    ica = mne.preprocessing.read_ica(path, verbose="ERROR")
    observed_exclusions = [int(index) for index in ica.exclude]
    errors: list[str] = []
    if int(ica.n_components_) != int(expected_components):
        errors.append("component_count")
    if observed_exclusions != list(expected_exclusions):
        errors.append("exclusion_ledger")
    if actual_hash != expected_sha256:
        errors.append("output_checksum")
    return {
        "reopened_with_mne": True,
        "valid": not errors,
        "errors": errors,
        "n_components": int(ica.n_components_),
        "exclude": observed_exclusions,
        "sha256": actual_hash,
        "checksum_matches": actual_hash == expected_sha256,
    }


def assess_terminal_result(result_dir: Path, expected_provenance: Mapping[str, Any]) -> dict[str, Any]:
    """Determine whether a published result is valid, stale, corrupt, or failed.

    Args:
        result_dir: Published self-contained recording directory.
        expected_provenance: Current source/config/code/environment fingerprint.

    Returns:
        Assessment with status, reason, and optional parsed manifest.

    Side effects:
        Reads JSON and hashes every recorded artifact. Writes nothing.
    """

    manifest_path = result_dir / "manifest.json"
    if not manifest_path.is_file():
        return {"cache_status": "incomplete", "reason": "manifest_missing"}
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return {"cache_status": "invalid", "reason": f"manifest_unreadable:{type(error).__name__}"}
    if manifest.get("provenance", {}).get("fingerprint") != expected_provenance["fingerprint"]:
        return {"cache_status": "stale", "reason": "provenance_fingerprint_changed", "manifest": manifest}
    if manifest.get("status") not in TERMINAL_CACHEABLE_STATUSES:
        return {"cache_status": "retry", "reason": f"noncacheable_status:{manifest.get('status')}", "manifest": manifest}
    for artifact in manifest.get("artifacts", []):
        path = result_dir / artifact["relative_path"]
        if not path.is_file():
            return {"cache_status": "invalid", "reason": f"artifact_missing:{artifact['relative_path']}", "manifest": manifest}
        if path.stat().st_size != artifact["size_bytes"] or sha256_file(path) != artifact["sha256"]:
            return {"cache_status": "invalid", "reason": f"artifact_hash_or_size_mismatch:{artifact['relative_path']}", "manifest": manifest}
    return {"cache_status": "valid", "reason": "terminal_result_and_hashes_valid", "manifest": manifest}


@dataclass
class ResultStore:
    """Manage one atomic self-contained per-recording result directory."""

    output_root: Path
    raw_root: Path
    key: str
    expected_provenance: Mapping[str, Any]
    force: bool = False
    recordings_dir: Path = field(init=False)
    history_dir: Path = field(init=False)
    final_dir: Path = field(init=False)
    work_dir: Path = field(init=False)

    def __post_init__(self) -> None:
        """Derive and validate all storage paths without writing them."""

        self.recordings_dir = self.output_root / "recordings"
        self.history_dir = self.output_root / "history"
        self.final_dir = self.recordings_dir / self.key
        self.work_dir = self.recordings_dir / f".{self.key}.incomplete"
        for path in (self.recordings_dir, self.history_dir, self.final_dir, self.work_dir):
            validate_derivative_path(path, self.output_root, self.raw_root)

    def _preserve(self, path: Path, reason: str) -> Path:
        """Atomically move an old local result into versioned history."""

        self.history_dir.mkdir(parents=True, exist_ok=True)
        destination = self.history_dir / f"{self.key}__{reason}__{history_stamp()}"
        validate_derivative_path(destination, self.output_root, self.raw_root)
        os.replace(path, destination)
        return destination

    def prepare(self) -> dict[str, Any]:
        """Skip a valid terminal result or create a clean incomplete workspace.

        Returns:
            Cache decision and paths. ``action`` is ``skip`` or ``process``.

        Side effects:
            Creates output directories. Preserves old incomplete, stale,
            invalid, failed, or explicitly force-recomputed results in history.
    """

        self.recordings_dir.mkdir(parents=True, exist_ok=True)
        history: list[str] = []
        if self.work_dir.exists():
            history.append(self._preserve(self.work_dir, "incomplete").as_posix())
        assessment: dict[str, Any] | None = None
        if self.final_dir.exists():
            assessment = assess_terminal_result(self.final_dir, self.expected_provenance)
            if assessment["cache_status"] == "valid" and not self.force:
                return {
                    "action": "skip",
                    "reason": assessment["reason"],
                    "result_dir": self.final_dir,
                    "manifest": assessment["manifest"],
                    "preserved_history": history,
                }
            reason = "force" if self.force else assessment["cache_status"]
            history.append(self._preserve(self.final_dir, reason).as_posix())
        self.work_dir.mkdir(parents=False, exist_ok=False)
        return {
            "action": "process",
            "reason": "force_recompute" if self.force else "no_valid_terminal_cache",
            "result_dir": self.work_dir,
            "previous_assessment": assessment,
            "preserved_history": history,
        }

    def publish(self) -> Path:
        """Atomically publish the complete work directory.

        Returns:
            Final recording directory.

        Side effects:
            Renames the incomplete directory after validating ``manifest.json``.
        """

        if self.final_dir.exists():
            raise FileExistsError(f"Final result unexpectedly exists: {self.final_dir}")
        if not (self.work_dir / "manifest.json").is_file():
            raise RuntimeError("Cannot publish a result without manifest.json.")
        os.replace(self.work_dir, self.final_dir)
        return self.final_dir
