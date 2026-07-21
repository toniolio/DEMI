#!/usr/bin/env python3
"""Construct and validate accepted DEMI trial-level EEG features.

This is the production stage-17 driver. It consumes the complete accepted
``tfr_v1`` response-onset and response-end arrays, verifies one-to-one linkage
to the frozen ``bdat2.rds`` behavioral authority, and constructs conventional
fixed-band/fixed-window features one recording at a time. All 30 physical
channel rows are preserved. Predeclared frontal-theta, posterior-alpha, and
task-hand-normalized motor-beta ROIs are then derived exclusively from the
completed channel rows.

Inputs:
    ``_Data/eeg/tfr_v1`` manifests, axes, row metadata, diagnostics, and
    read-only raw/dB arrays; ``_Private/old_rds/bdat2.rds``; the frozen bad
    imagery list; the historical offset export; and the tracked
    ``feature_config_v1.yaml``.

Outputs:
    Ignored ``_Data/eeg/features_v1`` Parquet channel/ROI tables, per-recording
    cache shards, feature/ROI/channel definitions, behavioral lineage and
    analysis-scope views, participant/descriptive summaries, source
    immutability evidence, manifests, validation, state, and Markdown.

This stage explicitly does not rerun convolution, create a baseline, use
31--40 Hz for feature rows, split beta, change eligibility, average
participants, center predictors, build a model design table, apply CSD, fit a
model, calculate a contrast, or interpret EEG effects.

Safety:
    TFR arrays are opened read-only with memory mapping and processed one
    recording at a time. Source array hashes, behavioral keys, accepted counts,
    task-hand mappings, endpoint membership, finite feature values, ROI
    equality, and source size/mtime snapshots fail closed. Recording shards and
    aggregate Parquet files publish atomically. An unchanged complete rerun
    hashes/reopens current artifacts and returns without rewriting them.

Run from the repository root:

    tools/run_features_v1.sh

Use ``--verify-current`` for a validation-only unchanged rerun. A bounded
``--recording demi_70_data`` invocation constructs or validates one cache shard
without publishing a complete aggregate surface.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import time
from typing import Any, Mapping, Sequence

import mne
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from event_epoch_eligibility import stable_frame_hash
from feature_construction import (
    EXPECTED_CHANNEL_COUNT,
    EXPECTED_CHANNEL_ROWS,
    EXPECTED_DURATION_WARNING_COUNT,
    EXPECTED_FILE49_COUNT,
    EXPECTED_PARTICIPANT_COUNT,
    EXPECTED_P4_TRIALS,
    EXPECTED_P60_TRIALS,
    EXPECTED_RECORDING_COUNT,
    EXPECTED_ROI_ROWS,
    EXPECTED_SHORT_0P5_COUNT,
    EXPECTED_SHORT_1P0_COUNT,
    EXPECTED_STRICT_CLEAN_COUNT,
    EXPECTED_TRIAL_COUNT,
    FEATURE_NAMESPACE,
    FEATURE_SCHEMA_VERSION,
    PRIMARY_CHANNELS,
    apply_task_hand_mapping,
    atomic_combine_parquet,
    atomic_write_json,
    atomic_write_parquet,
    atomic_write_text,
    build_channel_rows,
    compare_source_snapshots,
    definition_table,
    derive_roi_rows,
    feature_definitions,
    file_descriptor,
    frame_descriptor,
    inclusive_indices,
    load_and_validate_config,
    require_ignored_output_root,
    roi_definition_table,
    roi_equality_validation,
    sha256_file,
    source_snapshot,
    standard_1005_channel_table,
    summarize_transformed_cells,
)


CONFIG_PATH = Path("analysis/eeg_mne/feature_config_v1.yaml")
TFR_ROOT = Path("_Data/eeg/tfr_v1")
OUTPUT_ROOT = Path("_Data/eeg/features_v1")
BDAT_RDS = Path("_Private/old_rds/bdat2.rds")
BAD_IMAGERY_RDS = Path("_Private/old_rds/bad_imagery_trials.rds")
OFFSET_AUTHORITY = Path("_Data/behavior/event_offsets/event_offsets_old_compatible.csv")
RUN_MANIFEST = "feature_run_manifest.json"
RUN_STATE = "feature_run_state.json"
SUMMARY_JSON = "feature_summary.json"
SUMMARY_MD = "feature_summary.md"
VALIDATION_JSON = "validation/feature_validation.json"
SOURCE_IMMUTABILITY_JSON = "source_tfr_immutability.json"

BEHAVIOURAL_COLUMNS = (
    "participant", "db_id", "handedness", "exp_condition", "feedback_type",
    "figure_set", "block_num", "session_num", "condition", "trial_num",
    "figure_type", "figure_file", "trace_file", "stimulus_mt", "avg_velocity",
    "accuracy_rating", "vividness_rating", "mt", "mt_clip", "vresp", "group",
    "rep", "complexity", "error",
)
ARTIFACT_FLAGS = (
    "diagnostic_robust_logp2p_z_gt_6",
    "diagnostic_jump_gt_50uv_per_native_sample",
    "diagnostic_flat_p2p_lt_1uv",
)


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def print_progress(message: str) -> None:
    """Print one timestamped progress line immediately."""

    print(f"[{utc_now()}] {message}", flush=True)


def repo_root_from_script() -> Path:
    """Return and validate the repository root containing this driver."""

    root = Path(__file__).resolve().parents[2]
    if not (root / ".git").exists() or not (root / CONFIG_PATH).is_file():
        raise RuntimeError("could_not_locate_demi_repository_root")
    return root


def parse_args() -> argparse.Namespace:
    """Parse bounded recording, force, and validation-only options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recording", help="Construct or validate one exact TFR recording stem.")
    parser.add_argument("--force", action="store_true", help="Rebuild one exact --recording cache shard.")
    parser.add_argument(
        "--verify-current",
        action="store_true",
        help="Require and fully validate an unchanged complete current output.",
    )
    args = parser.parse_args()
    if args.force and not args.recording:
        parser.error("--force requires --recording")
    if args.verify_current and (args.force or args.recording):
        parser.error("--verify-current cannot be combined with --force or --recording")
    return args


def read_json(path: Path) -> dict[str, Any]:
    """Read one required JSON object."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"json_authority_not_object:{path}")
    return value


def content_fingerprint(value: Any) -> str:
    """Hash one JSON-serializable authority seed deterministically."""

    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_descriptor(path: Path, descriptor: Mapping[str, Any]) -> None:
    """Require a file to match its declared size and SHA-256."""

    if not path.is_file() or path.stat().st_size != int(descriptor["size_bytes"]):
        raise RuntimeError(f"file_size_mismatch:{path}")
    if sha256_file(path) != descriptor["sha256"]:
        raise RuntimeError(f"file_hash_mismatch:{path}")
    if descriptor.get("format") == "parquet" and "row_count" in descriptor:
        if pq.read_metadata(path).num_rows != int(descriptor["row_count"]):
            raise RuntimeError(f"parquet_row_count_mismatch:{path}")


def load_tfr_authority(repo_root: Path) -> dict[str, Any]:
    """Load and validate accepted TFR manifests, axes, metadata, and recordings."""

    tfr_root = repo_root / TFR_ROOT
    run_path = tfr_root / "tfr_run_manifest.json"
    run = read_json(run_path)
    if (
        run.get("status") != "complete"
        or run.get("tfr_namespace") != "tfr_v1"
        or int(run["storage"]["current_array_count"]) != 405
    ):
        raise RuntimeError("accepted_tfr_run_manifest_mismatch")

    axes: dict[str, pd.DataFrame] = {}
    for key in ("channels", "frequencies_cycles", "task_times"):
        descriptor = run["axes"][key]
        path = tfr_root / descriptor["relative_path"]
        validate_descriptor(path, descriptor)
        axes[key] = pd.read_csv(path)
    if axes["channels"]["channel"].tolist() != list(PRIMARY_CHANNELS):
        raise RuntimeError("accepted_tfr_channel_axis_mismatch")
    frequencies = axes["frequencies_cycles"]["frequency_hz"].to_numpy(dtype=float)
    times = axes["task_times"]["time_seconds"].to_numpy(dtype=float)
    if not np.array_equal(frequencies, np.arange(4.0, 41.0)) or len(times) != 201:
        raise RuntimeError("accepted_tfr_frequency_or_time_axis_mismatch")

    metadata: dict[str, pd.DataFrame] = {}
    family_manifest_paths: dict[str, Path] = {}
    for family in ("response_onset", "response_end"):
        family_descriptor = run["family_manifests"][family]
        family_path = tfr_root / family_descriptor["relative_path"]
        validate_descriptor(family_path, family_descriptor)
        family_manifest_paths[family] = family_path
        family_manifest = read_json(family_path)
        descriptor = family_manifest["metadata"]
        path = tfr_root / descriptor["relative_path"]
        validate_descriptor(path, descriptor)
        frame = pd.read_parquet(path)
        if len(frame) != EXPECTED_TRIAL_COUNT or stable_frame_hash(frame) != descriptor["content_sha256"]:
            raise RuntimeError(f"accepted_tfr_metadata_content_mismatch:{family}")
        metadata[family] = frame
    if not metadata["response_onset"]["canonical_event_key"].equals(
        metadata["response_end"]["canonical_event_key"]
    ):
        raise RuntimeError("accepted_tfr_onset_end_identity_mismatch")

    artifact_descriptor = run["artifact_diagnostics"]
    artifact_path = tfr_root / artifact_descriptor["relative_path"]
    validate_descriptor(artifact_path, artifact_descriptor)
    artifacts = pd.read_parquet(artifact_path)
    if len(artifacts) != EXPECTED_TRIAL_COUNT * 3:
        raise RuntimeError("accepted_tfr_artifact_diagnostic_count_mismatch")

    recordings: list[dict[str, Any]] = []
    source_paths: list[Path] = []
    for aggregate_recording in run["recordings"]:
        manifest_path = tfr_root / aggregate_recording["manifest_path"]
        manifest = read_json(manifest_path)
        if (
            manifest.get("status") != "complete"
            or manifest.get("fingerprint") != aggregate_recording["fingerprint"]
            or int(manifest["recording"]["row_count"]) != int(aggregate_recording["row_count"])
        ):
            raise RuntimeError(f"accepted_tfr_recording_manifest_mismatch:{manifest_path}")
        recording_root = manifest_path.parent
        array_paths = {}
        for label, descriptor in manifest["arrays"].items():
            path = recording_root / descriptor["relative_path"]
            if not path.is_file() or path.stat().st_size != int(descriptor["size_bytes"]):
                raise RuntimeError(f"accepted_tfr_array_size_mismatch:{path}")
            source_paths.append(path)
            array_paths[label] = path
        recordings.append(
            {
                "recording_stem": manifest["recording"]["recording_stem"],
                "source_filename": manifest["recording"]["source_filename"],
                "row_count": int(manifest["recording"]["row_count"]),
                "fingerprint": manifest["fingerprint"],
                "manifest_path": manifest_path,
                "manifest_sha256": sha256_file(manifest_path),
                "manifest": manifest,
                "array_paths": array_paths,
            }
        )
    if len(recordings) != EXPECTED_RECORDING_COUNT or sum(item["row_count"] for item in recordings) != EXPECTED_TRIAL_COUNT:
        raise RuntimeError("accepted_tfr_recording_surface_mismatch")
    return {
        "root": tfr_root,
        "run": run,
        "run_path": run_path,
        "run_sha256": sha256_file(run_path),
        "metadata": metadata,
        "artifacts": artifacts,
        "recordings": recordings,
        "source_paths": source_paths,
        "frequencies": frequencies,
        "times": times,
        "tfr_config_sha256": run["config"]["sha256"],
    }


def export_frozen_behaviour(repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read selected frozen RDS columns through base R without changing R state."""

    bdat_path = repo_root / BDAT_RDS
    bad_path = repo_root / BAD_IMAGERY_RDS
    with tempfile.TemporaryDirectory(prefix="demi_features_behaviour_", dir="/tmp") as temp:
        behaviour_csv = Path(temp) / "bdat2.csv"
        bad_csv = Path(temp) / "bad_imagery.csv"
        selected = ",".join(f'"{column}"' for column in BEHAVIOURAL_COLUMNS)
        expression = (
            f'x <- readRDS("{bdat_path}"); '
            f'x <- x[, c({selected})]; '
            f'write.csv(x, "{behaviour_csv}", row.names=FALSE, na=""); '
            f'y <- readRDS("{bad_path}"); '
            f'write.csv(y[, "figure_file", drop=FALSE], "{bad_csv}", row.names=FALSE, na="")'
        )
        environment = os.environ.copy()
        environment["R_PROFILE_USER"] = "/dev/null"
        subprocess.run(
            ["Rscript", "-e", expression], cwd=repo_root, env=environment, check=True
        )
        behaviour = pd.read_csv(behaviour_csv)
        bad_imagery = pd.read_csv(bad_csv)
    return behaviour, bad_imagery


def artifact_trial_fields(artifacts: pd.DataFrame) -> pd.DataFrame:
    """Return one row per trial with family-specific diagnostic flags."""

    parts = []
    for family in ("red_on", "response_onset", "response_end"):
        part = artifacts.loc[artifacts["family"].eq(family), ["canonical_event_key", *ARTIFACT_FLAGS]].copy()
        if len(part) != EXPECTED_TRIAL_COUNT or part["canonical_event_key"].duplicated().any():
            raise RuntimeError(f"artifact_family_identity_mismatch:{family}")
        part[f"artifact_{family}_any_flag"] = part[list(ARTIFACT_FLAGS)].any(axis=1)
        part = part.rename(
            columns={column: f"artifact_{family}_{column}" for column in ARTIFACT_FLAGS}
        )
        parts.append(part)
    result = parts[0]
    for part in parts[1:]:
        result = result.merge(part, on="canonical_event_key", how="inner", validate="one_to_one")
    return result


def link_frozen_behaviour(
    tfr: Mapping[str, Any],
    behaviour: pd.DataFrame,
    bad_imagery: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Create and validate the authoritative accepted trial-lineage table."""

    keys = ["participant", "session_num", "block_num", "trial_num"]
    duplicate_behaviour = int(behaviour.duplicated(keys, keep=False).sum())
    if len(behaviour) != 9_444 or duplicate_behaviour:
        raise RuntimeError(f"frozen_behaviour_authority_identity_mismatch:{len(behaviour)}:{duplicate_behaviour}")
    onset = tfr["metadata"]["response_onset"]
    end = tfr["metadata"]["response_end"]
    keep = [
        "tfr_row_index", "canonical_order_index", "canonical_event_key", "behavioural_id",
        "eeg_source_id", "source_recording_filename", "source_file_role", "task_filename",
        "offset_session", "offset_block", "offset_trial", "physical", "condition_semantics",
        "strict_clean_eligibility", "duration_warning_flag", "continuous_qc_warning",
        "interpolation_count", "event_policy_version", "ledger_version", "continuous_v2_run_id",
        "continuous_v2_manifest_id", "continuous_source_sha256",
    ]
    trial = onset[keep].copy()
    trial["response_onset_seconds"] = onset["anchor_time_seconds"].to_numpy(dtype=float)
    trial["response_end_seconds"] = end["anchor_time_seconds"].to_numpy(dtype=float)
    trial["response_duration_seconds"] = trial["response_end_seconds"] - trial["response_onset_seconds"]
    linked = trial.merge(
        behaviour,
        how="left",
        left_on=["behavioural_id", "offset_session", "offset_block", "offset_trial"],
        right_on=keys,
        validate="one_to_one",
        indicator=True,
    )
    unmatched = int(linked["_merge"].ne("both").sum())
    if unmatched or linked["canonical_event_key"].duplicated().any():
        raise RuntimeError(f"behavioural_lineage_failure:unmatched={unmatched}")
    linked = linked.drop(columns=["_merge"])

    physical_missing_mt_clip = int(linked.loc[linked["physical"], "mt_clip"].isna().sum())
    bad_figures = set(bad_imagery["figure_file"].dropna().astype(str))
    imagery_bad_overlap = int(
        linked.loc[~linked["physical"], "figure_file"].astype(str).isin(bad_figures).sum()
    )
    if physical_missing_mt_clip or imagery_bad_overlap:
        raise RuntimeError(
            f"frozen_behaviour_inclusion_failure:physical_missing={physical_missing_mt_clip}:imagery_bad={imagery_bad_overlap}"
        )

    linked = apply_task_hand_mapping(linked, config)
    recording_map = {item["source_filename"]: item["recording_stem"] for item in tfr["recordings"]}
    linked["recording_stem"] = linked["source_recording_filename"].map(recording_map)
    if linked["recording_stem"].isna().any():
        raise RuntimeError("recording_stem_linkage_failure")
    linked["source_shard"] = linked["recording_stem"]
    linked["behavioural_row_key"] = (
        linked["participant"].astype(str)
        + ":" + linked["session_num"].astype(str)
        + ":" + linked["block_num"].astype(str)
        + ":" + linked["trial_num"].astype(str)
    )
    linked["performed_condition"] = linked["condition"].astype(str)
    linked["familiarity_repetition_condition"] = linked["rep"].astype(str)
    linked["response_duration_lt_0p5"] = linked["response_duration_seconds"].lt(0.5)
    linked["response_duration_lt_1p0"] = linked["response_duration_seconds"].lt(1.0)
    linked["scope_all_accepted"] = True
    linked["scope_primary_blocks_1_5"] = linked["offset_block"].le(5)
    linked["scope_imagery_final_overt_bridge"] = (
        linked["group"].eq("imagery") & linked["physical"] & linked["offset_block"].eq(6)
    )
    linked["behavioural_inclusion_verified"] = True
    linked["overt_tracing_filter_provenance"] = np.where(
        linked["physical"], "frozen_bdat2_post_mt_clip_filter", "not_applicable_imagery"
    )
    linked["imagery_cleanup_provenance"] = np.where(
        linked["physical"], "not_applicable_overt", "frozen_bad_imagery_list_applied_before_bdat2"
    )
    linked = linked.merge(
        artifact_trial_fields(tfr["artifacts"]),
        on="canonical_event_key",
        how="left",
        validate="one_to_one",
    )
    if linked.filter(like="artifact_").isna().any().any():
        raise RuntimeError("artifact_trial_linkage_failure")
    linked = linked.sort_values("canonical_order_index", kind="stable").reset_index(drop=True)

    facts = {
        "accepted_eeg_trials": len(linked),
        "matched_frozen_behaviour_rows": len(linked),
        "unmatched_accepted_keys": unmatched,
        "duplicate_accepted_keys": int(linked["canonical_event_key"].duplicated().sum()),
        "duplicate_frozen_behaviour_keys": duplicate_behaviour,
        "frozen_behaviour_rows": len(behaviour),
        "physical_missing_mt_clip": physical_missing_mt_clip,
        "imagery_bad_figure_overlap": imagery_bad_overlap,
        "strict_clean_trials": int(linked["strict_clean_eligibility"].sum()),
        "duration_warning_trials": int(linked["duration_warning_flag"].sum()),
        "duration_lt_0p5_trials": int(linked["response_duration_lt_0p5"].sum()),
        "duration_lt_1p0_trials": int(linked["response_duration_lt_1p0"].sum()),
        "participants": int(linked["behavioural_id"].nunique()),
        "recordings": int(linked["recording_stem"].nunique()),
        "file49_trials": int(linked["eeg_source_id"].eq(49).sum()),
        "participant_4_trials": int(linked["behavioural_id"].eq(4).sum()),
        "participant_60_trials": int(linked["behavioural_id"].eq(60).sum()),
        "participant_70_trials": int(linked["behavioural_id"].eq(70).sum()),
        "participant_70_original_handedness": sorted(
            linked.loc[linked["behavioural_id"].eq(70), "original_handedness"].unique().tolist()
        ),
        "participant_70_analysis_hand": sorted(
            linked.loc[linked["behavioural_id"].eq(70), "analysis_hand"].unique().tolist()
        ),
        "participant_70_mapping_source": sorted(
            linked.loc[linked["behavioural_id"].eq(70), "mapping_source"].unique().tolist()
        ),
        "file54_1_present": bool(linked["source_recording_filename"].str.contains("54_1", regex=False).any()),
        "id86_present": bool(linked["behavioural_id"].eq(86).any() or linked["eeg_source_id"].eq(86).any()),
    }
    expected = {
        "accepted_eeg_trials": EXPECTED_TRIAL_COUNT,
        "matched_frozen_behaviour_rows": EXPECTED_TRIAL_COUNT,
        "unmatched_accepted_keys": 0,
        "duplicate_accepted_keys": 0,
        "duplicate_frozen_behaviour_keys": 0,
        "physical_missing_mt_clip": 0,
        "imagery_bad_figure_overlap": 0,
        "strict_clean_trials": EXPECTED_STRICT_CLEAN_COUNT,
        "duration_warning_trials": EXPECTED_DURATION_WARNING_COUNT,
        "duration_lt_0p5_trials": EXPECTED_SHORT_0P5_COUNT,
        "duration_lt_1p0_trials": EXPECTED_SHORT_1P0_COUNT,
        "participants": EXPECTED_PARTICIPANT_COUNT,
        "recordings": EXPECTED_RECORDING_COUNT,
        "file49_trials": EXPECTED_FILE49_COUNT,
        "participant_4_trials": EXPECTED_P4_TRIALS,
        "participant_60_trials": EXPECTED_P60_TRIALS,
        "participant_70_trials": 103,
        "file54_1_present": False,
        "id86_present": False,
    }
    if any(facts[key] != value for key, value in expected.items()):
        raise RuntimeError(f"accepted_lineage_count_mismatch:{facts}")
    if (
        facts["participant_70_original_handedness"] != ["a"]
        or facts["participant_70_analysis_hand"] != ["right"]
        or facts["participant_70_mapping_source"] != ["owner_recollection"]
    ):
        raise RuntimeError("participant_70_mapping_provenance_mismatch")
    return linked, facts


def authority_seed(
    repo_root: Path,
    config_sha256: str,
    tfr: Mapping[str, Any],
    linked: pd.DataFrame,
) -> tuple[dict[str, Any], dict[str, str]]:
    """Return complete feature authority seed and tracked code hashes."""

    code_paths = [
        repo_root / "analysis/eeg_mne/feature_construction.py",
        repo_root / "analysis/eeg_mne/17_construct_trial_level_features.py",
    ]
    code_files = {path.relative_to(repo_root).as_posix(): sha256_file(path) for path in code_paths}
    seed = {
        "schema_version": FEATURE_SCHEMA_VERSION,
        "namespace": FEATURE_NAMESPACE,
        "feature_config_sha256": config_sha256,
        "source_tfr_run_manifest_sha256": tfr["run_sha256"],
        "frozen_behaviour_rds_sha256": sha256_file(repo_root / BDAT_RDS),
        "frozen_bad_imagery_rds_sha256": sha256_file(repo_root / BAD_IMAGERY_RDS),
        "historical_offset_authority_sha256": sha256_file(repo_root / OFFSET_AUTHORITY),
        "linked_canonical_key_sha256": content_fingerprint(linked["canonical_event_key"].tolist()),
        "linked_task_hand_provenance_sha256": stable_frame_hash(
            linked[
                [
                    "canonical_event_key", "original_handedness", "analysis_hand",
                    "motor_laterality_mapping", "mapping_source", "mapping_override",
                ]
            ]
        ),
        "code_files": code_files,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "mne": mne.__version__,
        },
    }
    return seed, code_files


def recording_fingerprint(
    authority_fingerprint: str,
    recording: Mapping[str, Any],
    trials: pd.DataFrame,
) -> str:
    """Return the deterministic per-recording feature cache identity."""

    return content_fingerprint(
        {
            "authority_fingerprint": authority_fingerprint,
            "recording": recording["recording_stem"],
            "source_tfr_recording_fingerprint": recording["fingerprint"],
            "source_tfr_manifest_sha256": recording["manifest_sha256"],
            "trial_metadata_sha256": stable_frame_hash(trials),
        }
    )


def validate_recording_cache(
    directory: Path,
    expected_fingerprint: str,
    *,
    expected_trials: int,
) -> dict[str, Any] | None:
    """Hash/reopen one current recording shard or return None if not reusable."""

    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        return None
    manifest = read_json(manifest_path)
    if manifest.get("status") != "complete" or manifest.get("fingerprint") != expected_fingerprint:
        return None
    for label, expected_rows in (
        ("channel_features", expected_trials * EXPECTED_CHANNEL_COUNT * 11),
        ("roi_features", expected_trials * 13),
    ):
        descriptor = manifest["tables"][label]
        path = directory / descriptor["relative_path"]
        validate_descriptor(path, descriptor)
        if int(descriptor["row_count"]) != expected_rows:
            raise RuntimeError(f"cached_recording_row_count_mismatch:{directory}:{label}")
        values = pd.read_parquet(path, columns=["value_db", "value_log_power"])
        if not np.isfinite(values.to_numpy(dtype=float)).all():
            raise RuntimeError(f"cached_recording_nonfinite_features:{directory}:{label}")
    reused = dict(manifest)
    reused["invocation_cache_action"] = "reused_after_hash_reopen_validation"
    return reused


def archive_stale_recording(directory: Path, output_root: Path, reason: str) -> list[str]:
    """Move a stale generated shard into ignored history before replacement."""

    if not directory.exists():
        return []
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    target = output_root / "history" / stamp / directory.name
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(directory, target)
    atomic_write_json(target / "archive_reason.json", {"archived_at": utc_now(), "reason": reason})
    return [target.relative_to(output_root).as_posix()]


def archive_incomplete_recording_temporaries(output_root: Path, recording_stem: str) -> list[str]:
    """Preserve abandoned atomic-build directories before a new attempt."""

    archived = []
    for path in sorted((output_root / "recordings").glob(f".{recording_stem}.tmp-*")):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        target = output_root / "history" / stamp / f"{recording_stem}_incomplete"
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(path, target)
        atomic_write_json(
            target / "archive_reason.json",
            {"archived_at": utc_now(), "reason": "incomplete_atomic_build_recovered"},
        )
        archived.append(target.relative_to(output_root).as_posix())
    return archived


def build_recording(
    *,
    output_root: Path,
    recording: Mapping[str, Any],
    trials: pd.DataFrame,
    config: Mapping[str, Any],
    config_sha256: str,
    tfr: Mapping[str, Any],
    channels: pd.DataFrame,
    definitions: Sequence[Any],
    fingerprint: str,
    force: bool,
) -> dict[str, Any]:
    """Construct, validate, and atomically publish one recording feature shard."""

    final_dir = output_root / "recordings" / recording["recording_stem"]
    cached = validate_recording_cache(final_dir, fingerprint, expected_trials=len(trials))
    if cached is not None and not force:
        return cached
    archived = archive_incomplete_recording_temporaries(
        output_root, recording["recording_stem"]
    )
    archived.extend(archive_stale_recording(
        final_dir,
        output_root,
        "explicit_force_rebuild" if force else "cache_fingerprint_or_validation_mismatch",
    ))
    temp_dir = final_dir.parent / f".{final_dir.name}.tmp-{os.getpid()}-{time.time_ns()}"
    temp_dir.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    print_progress(f"construct recording={recording['recording_stem']} trials={len(trials)}")

    source_manifest = recording["manifest"]
    for label in ("response_onset_raw", "response_onset_db", "response_end_raw", "response_end_db"):
        path = recording["array_paths"][label]
        descriptor = source_manifest["arrays"][label]
        if sha256_file(path) != descriptor["sha256"]:
            raise RuntimeError(f"source_tfr_array_hash_mismatch:{path}")

    source_row_metadata = pd.read_parquet(recording["manifest_path"].parent / "row_metadata.parquet")
    if source_row_metadata["canonical_event_key"].tolist() != trials["canonical_event_key"].tolist():
        raise RuntimeError(f"recording_trial_order_mismatch:{recording['recording_stem']}")
    trials = trials.copy()
    trials["tfr_row_index"] = source_row_metadata["tfr_row_index"].to_numpy(dtype=int)
    trials["feature_config_sha256"] = config_sha256
    trials["source_tfr_run_manifest_sha256"] = tfr["run_sha256"]
    trials["source_tfr_config_sha256"] = tfr["tfr_config_sha256"]
    trials["source_tfr_recording_fingerprint"] = recording["fingerprint"]

    arrays = {
        family: {
            kind: np.load(
                recording["array_paths"][f"{family}_{kind}"], mmap_mode="r", allow_pickle=False
            )
            for kind in ("raw", "db")
        }
        for family in ("response_onset", "response_end")
    }
    expected_shape = (len(trials), EXPECTED_CHANNEL_COUNT, 37, 201)
    if any(array.shape != expected_shape or array.dtype != np.float32 for family in arrays.values() for array in family.values()):
        raise RuntimeError(f"source_tfr_array_contract_mismatch:{recording['recording_stem']}")

    channel_parts = []
    for definition in definitions:
        frequency_indices = inclusive_indices(
            tfr["frequencies"], definition.minimum_hz, definition.maximum_hz,
            expected_count=definition.expected_frequency_count,
            label=f"{definition.band}_frequency",
        )
        time_indices = inclusive_indices(
            tfr["times"], definition.minimum_seconds, definition.maximum_seconds,
            expected_count=definition.expected_time_sample_count,
            label=f"{definition.name}_time",
        )
        family = definition.source_family
        feature_trials = trials.copy()
        feature_trials["artifact_red_on_any_flag"] = feature_trials["artifact_red_on_any_flag"].astype(bool)
        feature_trials["artifact_source_any_flag"] = feature_trials[f"artifact_{family}_any_flag"].astype(bool)
        feature_trials["artifact_source_robust_logp2p_flag"] = feature_trials[
            f"artifact_{family}_diagnostic_robust_logp2p_z_gt_6"
        ].astype(bool)
        feature_trials["artifact_source_jump_flag"] = feature_trials[
            f"artifact_{family}_diagnostic_jump_gt_50uv_per_native_sample"
        ].astype(bool)
        feature_trials["artifact_source_flat_flag"] = feature_trials[
            f"artifact_{family}_diagnostic_flat_p2p_lt_1uv"
        ].astype(bool)
        values_db, values_log = summarize_transformed_cells(
            arrays[family]["db"], arrays[family]["raw"], frequency_indices, time_indices
        )
        channel_parts.append(
            build_channel_rows(
                feature_trials,
                channels,
                definition,
                values_db,
                values_log,
                source_raw_sha256=source_manifest["arrays"][f"{family}_raw"]["sha256"],
                source_db_sha256=source_manifest["arrays"][f"{family}_db"]["sha256"],
            )
        )
    channel_rows = pd.concat(channel_parts, ignore_index=True)
    channel_rows = channel_rows.sort_values(
        ["canonical_order_index", "feature_definition", "band", "channel_index"], kind="stable"
    ).reset_index(drop=True)
    if len(channel_rows) != len(trials) * EXPECTED_CHANNEL_COUNT * 11:
        raise RuntimeError("recording_channel_feature_count_mismatch")
    roi_rows = derive_roi_rows(channel_rows, config)
    roi_equality = roi_equality_validation(channel_rows, roi_rows)
    if (
        roi_equality["maximum_absolute_db_difference"] > 1e-12
        or roi_equality["maximum_absolute_log_power_difference"] > 1e-12
    ):
        raise RuntimeError("roi_derivation_not_exact_before_persistence")

    channel_path = temp_dir / "channel_features.parquet"
    roi_path = temp_dir / "roi_features.parquet"
    atomic_write_parquet(channel_path, channel_rows)
    atomic_write_parquet(roi_path, roi_rows)
    reopened_channel = pd.read_parquet(channel_path)
    reopened_roi = pd.read_parquet(roi_path)
    persisted_equality = roi_equality_validation(reopened_channel, reopened_roi)
    if (
        persisted_equality["maximum_absolute_db_difference"] > 1e-12
        or persisted_equality["maximum_absolute_log_power_difference"] > 1e-12
    ):
        raise RuntimeError("roi_derivation_not_exact_after_persistence")
    manifest = {
        "schema_version": FEATURE_SCHEMA_VERSION,
        "feature_namespace": FEATURE_NAMESPACE,
        "status": "complete",
        "created_at": utc_now(),
        "fingerprint": fingerprint,
        "recording": {
            "recording_stem": recording["recording_stem"],
            "source_filename": recording["source_filename"],
            "trial_count": len(trials),
        },
        "source_tfr": {
            "recording_fingerprint": recording["fingerprint"],
            "recording_manifest_sha256": recording["manifest_sha256"],
            "raw_db_array_sha256": {
                label: source_manifest["arrays"][label]["sha256"]
                for label in (
                    "response_onset_raw", "response_onset_db", "response_end_raw", "response_end_db"
                )
            },
        },
        "tables": {
            "channel_features": frame_descriptor(channel_path, reopened_channel, root=temp_dir),
            "roi_features": frame_descriptor(roi_path, reopened_roi, root=temp_dir),
        },
        "validation": {
            "channel_values_finite": bool(
                np.isfinite(reopened_channel[["value_db", "value_log_power"]]).all().all()
            ),
            "roi_values_finite": bool(
                np.isfinite(reopened_roi[["value_db", "value_log_power"]]).all().all()
            ),
            "roi_equality": persisted_equality,
            "zero_trial_loss": reopened_channel["canonical_event_key"].nunique() == len(trials),
        },
        "runtime_seconds": time.perf_counter() - started,
        "archived_previous_paths": archived,
    }
    atomic_write_json(temp_dir / "manifest.json", manifest)
    os.replace(temp_dir, final_dir)
    manifest["invocation_cache_action"] = "constructed_and_reopened"
    return manifest


def participant_summary(linked: pd.DataFrame) -> pd.DataFrame:
    """Return neutral participant trial-information counts."""

    return (
        linked.groupby(
            [
                "behavioural_id", "group", "original_handedness", "analysis_hand",
                "motor_laterality_mapping", "mapping_source", "mapping_override",
            ],
            observed=True,
            dropna=False,
        )
        .agg(
            trials=("canonical_event_key", "size"),
            recordings=("recording_stem", "nunique"),
            blocks=("offset_block", "nunique"),
            overt_trials=("physical", "sum"),
            primary_blocks_1_5_trials=("scope_primary_blocks_1_5", "sum"),
            final_overt_bridge_trials=("scope_imagery_final_overt_bridge", "sum"),
            strict_clean_trials=("strict_clean_eligibility", "sum"),
            duration_warning_trials=("duration_warning_flag", "sum"),
            duration_lt_0p5_trials=("response_duration_lt_0p5", "sum"),
            duration_lt_1p0_trials=("response_duration_lt_1p0", "sum"),
            accuracy_rating_nonmissing=("accuracy_rating", "count"),
            objective_error_nonmissing=("error", "count"),
        )
        .reset_index()
    )


def descriptive_feature_summary(channel_path: Path, roi_path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return neutral finite/range facts without condition or accuracy comparisons."""

    channel = pd.read_parquet(
        channel_path,
        columns=["feature_definition", "band", "value_db", "value_log_power"],
    )
    roi = pd.read_parquet(
        roi_path,
        columns=["roi_name", "feature_definition", "band", "value_db", "value_log_power"],
    )
    rows = []
    for surface, frame, groups in (
        ("channel", channel, ["feature_definition", "band"]),
        ("roi", roi, ["roi_name", "feature_definition", "band"]),
    ):
        for identity, subset in frame.groupby(groups, observed=True, sort=True):
            identity = identity if isinstance(identity, tuple) else (identity,)
            row = {"surface": surface, **dict(zip(groups, identity, strict=True))}
            row.update(
                {
                    "rows": len(subset),
                    "value_db_nonfinite": int((~np.isfinite(subset["value_db"])).sum()),
                    "value_log_power_nonfinite": int((~np.isfinite(subset["value_log_power"])).sum()),
                    "value_db_min": float(subset["value_db"].min()),
                    "value_db_max": float(subset["value_db"].max()),
                    "value_log_power_min": float(subset["value_log_power"].min()),
                    "value_log_power_max": float(subset["value_log_power"].max()),
                }
            )
            rows.append(row)
    correlation = float(np.corrcoef(channel["value_db"], channel["value_log_power"])[0, 1])
    correspondence = {
        "channel_rows": len(channel),
        "finite_pairs": int(np.isfinite(channel[["value_db", "value_log_power"]]).all(axis=1).sum()),
        "overall_descriptive_pearson_correlation": correlation,
        "role": "representation_correspondence_diagnostic_collapsed_across_condition_and_accuracy",
    }
    return pd.DataFrame(rows), correspondence


def write_summary_markdown(summary: Mapping[str, Any]) -> str:
    """Return concise neutral production summary Markdown."""

    return f"""# DEMI EEG features v1 summary

Status: **{summary['status']}**  
Completed: {summary['completed_at']}

- Accepted trials: {summary['trial_count']:,}
- Participants/recordings: {summary['participant_count']} / {summary['recording_count']}
- Channel feature rows: {summary['channel_row_count']:,}
- ROI feature rows: {summary['roi_row_count']:,}
- Strict-clean / duration-warning trials: {summary['strict_clean_count']:,} / {summary['duration_warning_count']}
- Durations below 0.5 / 1.0 s: {summary['duration_lt_0p5_count']} / {summary['duration_lt_1p0_count']}
- Behavioral linkage failures: 0
- Participant 70: original `a`; right task hand; owner-recollection override
- All channel and ROI feature values finite: yes
- ROI equality reproduced from channel rows: yes
- Source TFR mutation: none

The surface contains fixed theta, alpha, and beta summaries for the accepted
response-onset and response-end windows. Trial-matched dB is primary and
cell-level unnormalized log power is retained as the required sensitivity.
All 30 physical channels are preserved; ROI rows are equal-weight derivatives
of the channel table. No model, contrast, significance test, CSD derivative,
participant averaging, or scientific interpretation was produced.
"""


def complete_output_reusable(
    output_root: Path,
    authority_fingerprint: str,
    source_before: Sequence[Mapping[str, Any]],
) -> dict[str, Any] | None:
    """Fully hash/reopen an unchanged complete output without rewriting it."""

    path = output_root / RUN_MANIFEST
    if not path.is_file():
        return None
    manifest = read_json(path)
    if manifest.get("status") != "complete" or manifest.get("authority_fingerprint") != authority_fingerprint:
        return None
    for descriptor in manifest["outputs"].values():
        validate_descriptor(output_root / descriptor["relative_path"], descriptor)
    for recording in manifest["recordings"]:
        directory = output_root / "recordings" / recording["recording_stem"]
        cached = validate_recording_cache(
            directory,
            recording["fingerprint"],
            expected_trials=int(recording["trial_count"]),
        )
        if cached is None:
            raise RuntimeError(f"complete_recording_cache_not_reusable:{directory}")
    immutability = read_json(output_root / SOURCE_IMMUTABILITY_JSON)
    if immutability["before"] != list(source_before) or immutability["after"] != list(source_before):
        raise RuntimeError("current_source_tfr_snapshot_differs_from_completed_run")
    result = dict(manifest)
    result["invocation_cache_action"] = "complete_output_reused_after_hash_reopen_validation"
    return result


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute preflight, recording construction/reuse, aggregate validation, and publication."""

    repo_root = repo_root_from_script()
    os.chdir(repo_root)
    output_root = repo_root / OUTPUT_ROOT
    require_ignored_output_root(repo_root, output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    config, config_sha256 = load_and_validate_config(repo_root / CONFIG_PATH)
    tfr = load_tfr_authority(repo_root)
    source_before = source_snapshot(tfr["source_paths"], root=repo_root)
    behaviour, bad_imagery = export_frozen_behaviour(repo_root)
    linked, lineage_facts = link_frozen_behaviour(tfr, behaviour, bad_imagery, config)
    seed, code_files = authority_seed(repo_root, config_sha256, tfr, linked)
    authority_fingerprint = content_fingerprint(seed)

    if args.recording is None:
        reusable = complete_output_reusable(output_root, authority_fingerprint, source_before)
        if reusable is not None:
            print_progress(
                f"PASS features_v1 unchanged reuse trials={reusable['summary']['trial_count']} "
                f"channel_rows={reusable['summary']['channel_row_count']}"
            )
            return reusable
        if args.verify_current:
            raise RuntimeError("no_unchanged_complete_feature_output_to_verify")

    definitions = feature_definitions(config)
    channels = standard_1005_channel_table()
    feature_definitions_frame = definition_table(config, tfr["frequencies"], tfr["times"])
    roi_definitions_frame = roi_definition_table(config)
    selected_recordings = tfr["recordings"]
    if args.recording:
        selected_recordings = [item for item in selected_recordings if item["recording_stem"] == args.recording]
        if len(selected_recordings) != 1:
            raise RuntimeError(f"recording_selection_not_exact:{args.recording}")

    run_state = {
        "schema_version": FEATURE_SCHEMA_VERSION,
        "feature_namespace": FEATURE_NAMESPACE,
        "status": "running",
        "started_at": utc_now(),
        "authority_fingerprint": authority_fingerprint,
        "selected_recordings": [item["recording_stem"] for item in selected_recordings],
        "completed_recordings": 0,
    }
    atomic_write_json(output_root / RUN_STATE, run_state)
    manifests = []
    started = time.perf_counter()
    for index, recording in enumerate(selected_recordings, start=1):
        trials = linked.loc[linked["recording_stem"].eq(recording["recording_stem"])].copy()
        fingerprint = recording_fingerprint(authority_fingerprint, recording, trials)
        manifest = build_recording(
            output_root=output_root,
            recording=recording,
            trials=trials,
            config=config,
            config_sha256=config_sha256,
            tfr=tfr,
            channels=channels,
            definitions=definitions,
            fingerprint=fingerprint,
            force=args.force,
        )
        manifests.append(manifest)
        run_state["completed_recordings"] = index
        run_state["last_recording"] = recording["recording_stem"]
        run_state["elapsed_seconds"] = time.perf_counter() - started
        atomic_write_json(output_root / RUN_STATE, run_state)
        print_progress(
            f"recording {index}/{len(selected_recordings)} {recording['recording_stem']} "
            f"cache={manifest['invocation_cache_action']} runtime={manifest['runtime_seconds']:.2f}s"
        )

    source_after_recordings = source_snapshot(tfr["source_paths"], root=repo_root)
    compare_source_snapshots(source_before, source_after_recordings)
    if args.recording:
        run_state["status"] = "partial_recording_complete"
        run_state["completed_at"] = utc_now()
        atomic_write_json(output_root / RUN_STATE, run_state)
        return {
            "status": "partial_recording_complete",
            "recording": args.recording,
            "manifest": manifests[0],
        }

    recording_order = [item["recording_stem"] for item in tfr["recordings"]]
    channel_shards = [output_root / "recordings" / name / "channel_features.parquet" for name in recording_order]
    roi_shards = [output_root / "recordings" / name / "roi_features.parquet" for name in recording_order]
    channel_path = output_root / "channel_features.parquet"
    roi_path = output_root / "roi_features.parquet"
    channel_descriptor = atomic_combine_parquet(channel_shards, channel_path)
    roi_descriptor = atomic_combine_parquet(roi_shards, roi_path)
    if channel_descriptor["row_count"] != EXPECTED_CHANNEL_ROWS or roi_descriptor["row_count"] != EXPECTED_ROI_ROWS:
        raise RuntimeError("aggregate_feature_row_count_mismatch")

    config_snapshot = output_root / "configuration/feature_config_v1.yaml"
    atomic_write_text(config_snapshot, (repo_root / CONFIG_PATH).read_text(encoding="utf-8"))
    if sha256_file(config_snapshot) != config_sha256:
        raise RuntimeError("feature_config_snapshot_hash_mismatch")

    views: dict[str, tuple[Path, pd.DataFrame]] = {}
    lineage_columns = [
        column for column in linked.columns
        if not column.startswith("artifact_") or column.endswith("_any_flag")
    ]
    views["behavioural_lineage"] = (
        output_root / "metadata/behavioural_lineage.parquet", linked[lineage_columns].copy()
    )
    views["feature_definitions"] = (
        output_root / "definitions/feature_definitions.parquet", feature_definitions_frame
    )
    views["roi_definitions"] = (
        output_root / "definitions/roi_definitions.parquet", roi_definitions_frame
    )
    views["channels"] = (output_root / "definitions/channels_standard_1005.parquet", channels)
    views["analysis_scopes"] = (
        output_root / "views/analysis_scope_indices.parquet",
        linked[
            [
                "canonical_order_index", "canonical_event_key", "scope_all_accepted",
                "scope_primary_blocks_1_5", "scope_imagery_final_overt_bridge",
            ]
        ].copy(),
    )
    views["strict_clean"] = (
        output_root / "views/strict_clean_row_indices.parquet",
        linked.loc[linked["strict_clean_eligibility"], ["canonical_order_index", "canonical_event_key"]].copy(),
    )
    views["duration_warning"] = (
        output_root / "views/duration_warning_rows.parquet",
        linked.loc[linked["duration_warning_flag"]].copy(),
    )
    views["short_duration"] = (
        output_root / "views/short_duration_rows.parquet",
        linked.loc[linked["response_duration_lt_1p0"]].copy(),
    )
    views["artifact_qc"] = (
        output_root / "views/artifact_qc_metadata.parquet", tfr["artifacts"].copy()
    )
    participant = participant_summary(linked)
    views["participant_summary"] = (
        output_root / "summaries/participant_trial_counts.parquet", participant
    )
    design = (
        linked.groupby(["group", "performed_condition", "offset_block"], observed=True)
        .agg(
            trials=("canonical_event_key", "size"),
            participants=("behavioural_id", "nunique"),
            recordings=("recording_stem", "nunique"),
            strict_clean_trials=("strict_clean_eligibility", "sum"),
            duration_warning_trials=("duration_warning_flag", "sum"),
        )
        .reset_index()
    )
    views["design_counts"] = (output_root / "summaries/design_trial_counts.parquet", design)
    output_descriptors: dict[str, dict[str, Any]] = {}
    for label, (path, frame) in views.items():
        atomic_write_parquet(path, frame)
        output_descriptors[label] = frame_descriptor(path, frame, root=output_root)

    descriptive, representation_correspondence = descriptive_feature_summary(channel_path, roi_path)
    descriptive_path = output_root / "summaries/descriptive_feature_validation.parquet"
    atomic_write_parquet(descriptive_path, descriptive)
    output_descriptors["descriptive_feature_validation"] = frame_descriptor(
        descriptive_path, descriptive, root=output_root
    )

    channel_descriptor["relative_path"] = channel_path.relative_to(output_root).as_posix()
    roi_descriptor["relative_path"] = roi_path.relative_to(output_root).as_posix()
    output_descriptors["channel_features"] = channel_descriptor
    output_descriptors["roi_features"] = roi_descriptor
    output_descriptors["feature_config_snapshot"] = file_descriptor(config_snapshot, root=output_root)

    source_after = source_snapshot(tfr["source_paths"], root=repo_root)
    compare_source_snapshots(source_before, source_after)
    source_fingerprint = content_fingerprint(source_before)
    immutability = {
        "status": "unchanged",
        "source": "accepted_tfr_v1_arrays",
        "array_count": len(source_before),
        "total_bytes": sum(int(item["size_bytes"]) for item in source_before),
        "size_mtime_fingerprint": source_fingerprint,
        "before": source_before,
        "after": source_after,
        "source_tfr_run_manifest_sha256": tfr["run_sha256"],
    }
    immutability_path = output_root / SOURCE_IMMUTABILITY_JSON
    atomic_write_json(immutability_path, immutability)
    output_descriptors["source_tfr_immutability"] = file_descriptor(
        immutability_path, root=output_root
    )

    validation = {
        "status": "pass",
        "lineage": lineage_facts,
        "channel_row_count": pq.read_metadata(channel_path).num_rows,
        "roi_row_count": pq.read_metadata(roi_path).num_rows,
        "all_channel_values_finite": int(descriptive.loc[descriptive["surface"].eq("channel"), ["value_db_nonfinite", "value_log_power_nonfinite"]].sum().sum()) == 0,
        "all_roi_values_finite": int(descriptive.loc[descriptive["surface"].eq("roi"), ["value_db_nonfinite", "value_log_power_nonfinite"]].sum().sum()) == 0,
        "recording_roi_equality": {
            item["recording"]["recording_stem"]: item["validation"]["roi_equality"]
            for item in manifests
        },
        "all_roi_values_reproduce_from_channel_rows": all(
            item["validation"]["roi_equality"]["maximum_absolute_db_difference"] <= 1e-12
            and item["validation"]["roi_equality"]["maximum_absolute_log_power_difference"] <= 1e-12
            for item in manifests
        ),
        "representation_correspondence": representation_correspondence,
        "participant_4_retained": bool(participant["behavioural_id"].eq(4).any()),
        "participant_60_retained": bool(participant["behavioural_id"].eq(60).any()),
        "participant_70_task_hand_mapping": {
            "original_handedness": "a",
            "analysis_hand": "right",
            "motor_laterality_mapping": "right_hand_task",
            "mapping_source": "owner_recollection",
            "mapping_override": True,
        },
        "no_31_40_hz_feature_input": bool(feature_definitions_frame["frequency_maximum_hz"].le(30.0).all()),
        "no_low_high_beta_split": set(feature_definitions_frame["band"]) == {"theta", "alpha", "beta"},
        "no_participant_averaging_authoritative_output": True,
        "no_eligibility_change": linked["canonical_event_key"].nunique() == EXPECTED_TRIAL_COUNT,
        "source_tfr_unchanged": True,
        "forbidden_outputs_absent": [
            "models", "contrasts", "p_values", "posterior_intervals", "gamm", "cluster_tests", "csd"
        ],
    }
    if not all(
        validation[key]
        for key in (
            "all_channel_values_finite", "all_roi_values_finite",
            "all_roi_values_reproduce_from_channel_rows", "participant_4_retained",
            "participant_60_retained", "no_31_40_hz_feature_input",
            "no_low_high_beta_split", "no_participant_averaging_authoritative_output",
            "no_eligibility_change", "source_tfr_unchanged",
        )
    ):
        raise RuntimeError(f"aggregate_feature_validation_failed:{validation}")
    validation_path = output_root / VALIDATION_JSON
    atomic_write_json(validation_path, validation)
    output_descriptors["validation"] = file_descriptor(validation_path, root=output_root)

    completed_at = utc_now()
    summary = {
        "status": "complete",
        "completed_at": completed_at,
        "trial_count": EXPECTED_TRIAL_COUNT,
        "participant_count": EXPECTED_PARTICIPANT_COUNT,
        "recording_count": EXPECTED_RECORDING_COUNT,
        "channel_row_count": EXPECTED_CHANNEL_ROWS,
        "roi_row_count": EXPECTED_ROI_ROWS,
        "strict_clean_count": EXPECTED_STRICT_CLEAN_COUNT,
        "duration_warning_count": EXPECTED_DURATION_WARNING_COUNT,
        "duration_lt_0p5_count": EXPECTED_SHORT_0P5_COUNT,
        "duration_lt_1p0_count": EXPECTED_SHORT_1P0_COUNT,
        "runtime_seconds": time.perf_counter() - started,
        "storage_bytes": sum(descriptor["size_bytes"] for descriptor in output_descriptors.values()),
        "cache_actions": dict(Counter(item["invocation_cache_action"] for item in manifests)),
        "explicit_nonoperations": [
            "no TFR mutation or convolution", "no new baseline", "no 31-40 Hz feature rows",
            "no trial or participant loss", "no participant averaging", "no model design table",
            "no centering or standardization", "no inferential model or contrast", "no CSD",
        ],
    }
    summary_json_path = output_root / SUMMARY_JSON
    summary_md_path = output_root / SUMMARY_MD
    atomic_write_json(summary_json_path, summary)
    atomic_write_text(summary_md_path, write_summary_markdown(summary))
    output_descriptors["summary_json"] = file_descriptor(summary_json_path, root=output_root)
    output_descriptors["summary_markdown"] = file_descriptor(summary_md_path, root=output_root)

    run_manifest = {
        "schema_version": FEATURE_SCHEMA_VERSION,
        "feature_namespace": FEATURE_NAMESPACE,
        "status": "complete",
        "started_at": run_state["started_at"],
        "completed_at": completed_at,
        "authority_fingerprint": authority_fingerprint,
        "authority": seed,
        "code_files": code_files,
        "source_tfr": {
            "root": TFR_ROOT.as_posix(),
            "run_manifest_sha256": tfr["run_sha256"],
            "config_sha256": tfr["tfr_config_sha256"],
            "immutability_path": SOURCE_IMMUTABILITY_JSON,
        },
        "frozen_behaviour": {
            "path": BDAT_RDS.as_posix(),
            "sha256": sha256_file(repo_root / BDAT_RDS),
            "row_count": len(behaviour),
            "accepted_matched_rows": EXPECTED_TRIAL_COUNT,
        },
        "outputs": output_descriptors,
        "recordings": [
            {
                "recording_stem": item["recording"]["recording_stem"],
                "trial_count": item["recording"]["trial_count"],
                "fingerprint": item["fingerprint"],
                "manifest_path": (
                    Path("recordings") / item["recording"]["recording_stem"] / "manifest.json"
                ).as_posix(),
                "cache_action": item["invocation_cache_action"],
            }
            for item in manifests
        ],
        "summary": summary,
    }
    atomic_write_json(output_root / RUN_MANIFEST, run_manifest)
    run_state.update(
        {
            "status": "complete",
            "completed_at": completed_at,
            "elapsed_seconds": summary["runtime_seconds"],
            "channel_row_count": EXPECTED_CHANNEL_ROWS,
            "roi_row_count": EXPECTED_ROI_ROWS,
        }
    )
    atomic_write_json(output_root / RUN_STATE, run_state)
    print_progress(
        f"PASS features_v1 trials={EXPECTED_TRIAL_COUNT} channel_rows={EXPECTED_CHANNEL_ROWS} "
        f"roi_rows={EXPECTED_ROI_ROWS} runtime={summary['runtime_seconds']:.1f}s"
    )
    return run_manifest


def main() -> None:
    """Run stage 17 and persist a failed state without masking the exception."""

    args = parse_args()
    try:
        run(args)
    except Exception as error:
        try:
            repo_root = repo_root_from_script()
            output_root = repo_root / OUTPUT_ROOT
            if output_root.exists():
                atomic_write_json(
                    output_root / RUN_STATE,
                    {
                        "schema_version": FEATURE_SCHEMA_VERSION,
                        "feature_namespace": FEATURE_NAMESPACE,
                        "status": "failed",
                        "failed_at": utc_now(),
                        "error_type": type(error).__name__,
                        "error": str(error),
                    },
                )
        finally:
            raise


if __name__ == "__main__":
    main()
