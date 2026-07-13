"""Source validation, channel typing, and montage stages for DEMI EEG.

The functions in this module open one EDF read-only, validate its original
fixed-width signal labels, load the continuous recording into memory, apply the
accepted channel types, and attach the approved ``standard_1005`` montage.
They return explicit evidence dictionaries for the per-recording manifest.

No function writes, renames, repairs, or otherwise mutates a source EDF. This
module does not filter data, identify bad channels, construct epochs, compute
CSD, run AutoReject, or make recording-inclusion decisions.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import mne
import numpy as np

from channel_qc import read_edf_signal_headers

from .contracts import (
    ALL_RETAINED_CHANNELS,
    CHANNEL_TYPE_BY_NAME,
    EEG_TARGET_CHANNELS,
    sha256_file,
)


def validate_source_path(path: Path, raw_root: Path) -> None:
    """Require an EDF directly below the configured immutable raw directory.

    Args:
        path: Proposed source recording.
        raw_root: Configured raw EDF directory.

    Returns:
        None.

    Raises:
        ValueError: If the path escapes the raw root or is not an EDF.

    Side effects:
        Resolves paths and inspects file existence; writes nothing.
    """

    resolved_root = raw_root.resolve()
    resolved_path = path.resolve()
    if resolved_path.parent != resolved_root:
        raise ValueError(f"Source EDF must be directly inside {resolved_root}: {path}")
    if resolved_path.suffix.lower() != ".edf":
        raise ValueError(f"Source recording must be an EDF: {path}")
    if not resolved_path.is_file():
        raise FileNotFoundError(f"Source EDF is missing: {path}")


def _physical_signal_labels(path: Path) -> list[str]:
    """Return original EDF signal labels excluding the annotation signal."""

    labels = [row.label.strip() for row in read_edf_signal_headers(path)]
    return [label for label in labels if label.upper() != "EDF ANNOTATIONS"]


def validate_required_channels(labels: list[str]) -> dict[str, Any]:
    """Validate the exact 37-channel DEMI recording surface.

    Args:
        labels: Original physical signal labels from the EDF header.

    Returns:
        Manifest evidence for exact names, counts, and order.

    Raises:
        ValueError: If any required channel is absent, duplicated, unexpected,
            misspelled, or reordered.

    Side effects:
        None.
    """

    counts = Counter(labels)
    duplicates = sorted(label for label, count in counts.items() if count > 1)
    missing = sorted(set(ALL_RETAINED_CHANNELS) - set(labels))
    unexpected = sorted(set(labels) - set(ALL_RETAINED_CHANNELS))
    if duplicates or missing or unexpected or len(labels) != len(ALL_RETAINED_CHANNELS):
        raise ValueError(
            "Invalid required channel surface: "
            f"missing={missing}, duplicated={duplicates}, unexpected={unexpected}, "
            f"observed_count={len(labels)}"
        )
    if tuple(labels) != ALL_RETAINED_CHANNELS:
        raise ValueError(
            "Required channel order differs from the accepted EDF contract: "
            f"observed={labels}"
        )
    return {
        "required_channel_count": len(ALL_RETAINED_CHANNELS),
        "observed_channel_count": len(labels),
        "channel_names": list(labels),
        "missing_channels": [],
        "duplicated_channels": [],
        "unexpected_channels": [],
        "channel_order_matches_contract": True,
    }


def read_and_validate_source(path: Path, raw_root: Path) -> tuple[mne.io.BaseRaw, dict[str, Any]]:
    """Hash, validate, and preload one immutable source EDF.

    Args:
        path: Source EDF.
        raw_root: Configured raw EDF directory.

    Returns:
        ``(raw, evidence)`` where ``raw`` is the full preloaded recording and
        ``evidence`` records source identity, checksum, and acquisition shape.

    Side effects:
        Reads the EDF header and all signal samples. It writes nothing.
    """

    validate_source_path(path, raw_root)
    stat_before = path.stat()
    source_hash = sha256_file(path)
    labels = _physical_signal_labels(path)
    channel_evidence = validate_required_channels(labels)
    raw = mne.io.read_raw_edf(path, preload=True, infer_types=False, verbose="ERROR")
    if raw.ch_names != labels:
        raise ValueError(
            "MNE channel labels do not match the original physical EDF labels: "
            f"header={labels}, mne={raw.ch_names}"
        )
    if not np.isfinite(float(raw.info["sfreq"])) or float(raw.info["sfreq"]) <= 0:
        raise ValueError("Source EDF has an invalid sampling frequency.")
    if raw.n_times < 2:
        raise ValueError("Source EDF has too few samples for continuous preprocessing.")
    evidence = {
        "source_filename": path.name,
        "source_path": path.as_posix(),
        "source_sha256": source_hash,
        "source_size_bytes": stat_before.st_size,
        "source_mtime_ns": stat_before.st_mtime_ns,
        "sampling_frequency_hz": float(raw.info["sfreq"]),
        "sample_count": int(raw.n_times),
        "duration_seconds": float(raw.n_times / raw.info["sfreq"]),
        "mne_reader": "mne.io.read_raw_edf",
        "preloaded": True,
        **channel_evidence,
    }
    return raw, evidence


def apply_channel_types(raw: mne.io.BaseRaw) -> dict[str, Any]:
    """Apply and verify the accepted 32 EEG, 2 EOG, 2 EMG, 1 stim typing.

    Args:
        raw: In-memory full recording.

    Returns:
        Manifest evidence with the exact type assigned to every channel.

    Side effects:
        Mutates only the in-memory MNE ``Raw`` metadata.
    """

    if tuple(raw.ch_names) != ALL_RETAINED_CHANNELS:
        raise ValueError("Channel typing received a recording outside the accepted surface.")
    raw.set_channel_types(dict(CHANNEL_TYPE_BY_NAME), on_unit_change="ignore", verbose="ERROR")
    observed = dict(zip(raw.ch_names, raw.get_channel_types()))
    if observed != dict(CHANNEL_TYPE_BY_NAME):
        raise ValueError(f"Channel typing verification failed: {observed}")
    return {
        "channel_types": observed,
        "type_counts": dict(Counter(observed.values())),
        "m1_m2_retained_as_eeg_targets": all(
            channel in raw.ch_names and observed[channel] == "eeg" for channel in ("M1", "M2")
        ),
    }


def uppercase_standard_1005() -> mne.channels.DigMontage:
    """Return ``standard_1005`` with uppercase labels matching the DEMI EDFs."""

    montage = mne.channels.make_standard_montage("standard_1005")
    rename = {name: name.upper() for name in montage.ch_names}
    if len(set(rename.values())) != len(rename):
        raise ValueError("Uppercasing standard_1005 produced duplicate labels.")
    montage.rename_channels(rename)
    return montage


def apply_montage(raw: mne.io.BaseRaw) -> dict[str, Any]:
    """Attach and validate the approved ``standard_1005`` EEG coordinates.

    Args:
        raw: Typed in-memory full recording.

    Returns:
        Montage validation evidence for every retained EEG target.

    Side effects:
        Mutates only the in-memory MNE montage/digitization metadata.
    """

    montage = uppercase_standard_1005()
    raw.set_montage(montage, on_missing="raise", match_case=True, verbose="ERROR")
    positions = raw.get_montage().get_positions()["ch_pos"] if raw.get_montage() else {}
    missing = [channel for channel in EEG_TARGET_CHANNELS if channel not in positions]
    nonfinite = [
        channel
        for channel in EEG_TARGET_CHANNELS
        if channel in positions and not np.isfinite(np.asarray(positions[channel], dtype=float)).all()
    ]
    if missing or nonfinite:
        raise ValueError(
            f"standard_1005 montage validation failed: missing={missing}, nonfinite={nonfinite}"
        )
    return {
        "montage": "standard_1005",
        "coordinate_frame": "head",
        "eeg_targets_with_positions": list(EEG_TARGET_CHANNELS),
        "position_count": len(EEG_TARGET_CHANNELS),
        "historical_besa_coordinate_source_used": False,
    }


def verify_source_unchanged(path: Path, original: dict[str, Any]) -> dict[str, Any]:
    """Re-hash a source EDF and prove its size and timestamp stayed unchanged.

    Args:
        path: Source EDF.
        original: Evidence returned by :func:`read_and_validate_source`.

    Returns:
        A manifest dictionary with before/after comparisons.

    Raises:
        RuntimeError: If checksum, size, or modification time changed.

    Side effects:
        Reads the EDF to compute SHA-256. It writes nothing.
    """

    stat_after = path.stat()
    hash_after = sha256_file(path)
    evidence = {
        "source_sha256_before": original["source_sha256"],
        "source_sha256_after": hash_after,
        "source_size_bytes_before": original["source_size_bytes"],
        "source_size_bytes_after": stat_after.st_size,
        "source_mtime_ns_before": original["source_mtime_ns"],
        "source_mtime_ns_after": stat_after.st_mtime_ns,
    }
    evidence["unchanged"] = (
        evidence["source_sha256_before"] == evidence["source_sha256_after"]
        and evidence["source_size_bytes_before"] == evidence["source_size_bytes_after"]
        and evidence["source_mtime_ns_before"] == evidence["source_mtime_ns_after"]
    )
    if not evidence["unchanged"]:
        raise RuntimeError(f"Raw EDF changed during preprocessing: {path}")
    return evidence
