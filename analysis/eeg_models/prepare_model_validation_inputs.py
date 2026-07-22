#!/usr/bin/env python3
"""Prepare isolated inputs for DEMI EEG synthetic model/prior validation.

This script exists to enforce a hard boundary between two authorized reads of
the accepted Stage 18 model tables. The prior-scale path projects only pooled
``value_db`` and returns median/MAD scale constants. The synthetic-design path
projects only explicitly allowlisted predictor, grouping, and key columns,
then verifies that no accepted EEG outcome column survived. It never combines
the two projections and never computes an EEG association.

Inputs are the accepted ``_Data/eeg/model_tables_v1`` views, validation,
manifest, and registries. The ``prepare`` command writes safe CSV design files,
scale/preflight JSON, and a source snapshot beneath a caller-provided staging
directory. The ``verify`` command compares the current accepted sources with
that snapshot and writes source-immutability evidence.

The script does not fit a model, generate a scientific estimand, inspect EEG by
condition/rating/familiarity/participant/phase, modify accepted production
files, or write outside the requested ignored validation staging directory.

Run from the repository root, normally through the Stage 19 R driver:

``python analysis/eeg_models/prepare_model_validation_inputs.py prepare --output-root PATH``
``python analysis/eeg_models/prepare_model_validation_inputs.py verify --before PATH --evidence PATH``
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any, Iterable

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


REPO_ROOT = Path(__file__).resolve().parents[2]
STAGE18_ROOT = REPO_ROOT / "_Data/eeg/model_tables_v1"
STAGE18_AUTHORITY_FINGERPRINT = (
    "685cdf04bb9d638b69e7b216f656ae4088659c7570bc3b451f8daafc80a8f43b"
)
MODEL_VIEWS = {
    "H1": "views/h1_theta_primary.parquet",
    "H2": "views/h2_alpha_primary.parquet",
    "H3": "views/h3_beta_primary.parquet",
}
EXPECTED_ROWS = {"H1": 7_307, "H2": 7_307, "H3": 14_614}
DESIGN_COLUMNS = {
    "H1": (
        "canonical_event_key",
        "participant_id",
        "condition_c",
        "accuracy_rating_within",
        "accuracy_rating_between",
        "familiarity_c",
    ),
    "H2": (
        "canonical_event_key",
        "participant_id",
        "condition_c",
        "accuracy_rating_within",
        "accuracy_rating_between",
        "familiarity_c",
    ),
    "H3": (
        "canonical_event_key",
        "participant_id",
        "condition_c",
        "phase_c",
        "familiarity_c",
    ),
}
KNOWN_OBSERVED_EEG_OUTCOMES = frozenset({"value_db", "value_log_power"})
PRIMARY_ESTIMANDS = ("T1", "T2", "T3", "A1", "A2", "A3", "B1", "B2", "B3")
SUPPLEMENTARY_ESTIMANDS = ("B4", "B5", "B6")


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for ``path`` without modifying it."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_read(path: Path) -> dict[str, Any]:
    """Read and return a JSON object from ``path``."""

    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def json_write(path: Path, value: Any) -> None:
    """Write deterministic, human-readable JSON and create parent folders."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def observed_outcome_columns(columns: Iterable[str]) -> tuple[str, ...]:
    """Return column names that could contain copied accepted EEG outcomes."""

    found: list[str] = []
    for column in columns:
        lowered = column.lower()
        if (
            column in KNOWN_OBSERVED_EEG_OUTCOMES
            or "value_db" in lowered
            or "value_log_power" in lowered
            or lowered.startswith("observed_eeg")
            or lowered.startswith("accepted_eeg")
        ):
            found.append(column)
    return tuple(found)


def assert_synthetic_fit_input_safe(frame: pd.DataFrame) -> None:
    """Reject a proposed synthetic fit frame containing observed EEG outcomes.

    The check is deliberately independent of whether a ``value_synthetic``
    column exists so callers can validate a design immediately after loading
    and again after injecting a synthetic response.
    """

    forbidden = observed_outcome_columns(frame.columns)
    if forbidden:
        raise ValueError(f"synthetic_fit_input_contains_observed_eeg:{','.join(forbidden)}")


def pooled_scale_from_value_only(path: Path) -> dict[str, Any]:
    """Calculate pooled median and robust scale from a ``value_db`` projection.

    No predictor, grouping, key, or other EEG outcome column is requested from
    Arrow. ``s_y`` uses raw median absolute deviation multiplied once by
    1.4826, matching the accepted contract rather than R's already-scaled
    default ``mad`` behavior.
    """

    projected_columns = ("value_db",)
    table = pq.read_table(path, columns=list(projected_columns))
    if tuple(table.column_names) != projected_columns:
        raise RuntimeError("prior_scale_projection_not_value_db_only")
    values = table.column("value_db").to_numpy(zero_copy_only=False)
    if values.size == 0 or not np.isfinite(values).all():
        raise RuntimeError("prior_scale_values_missing_or_nonfinite")
    median = float(np.median(values))
    raw_mad = float(np.median(np.abs(values - median)))
    scale = float(1.4826 * raw_mad)
    if not np.isfinite(scale) or scale <= 0:
        raise RuntimeError("prior_scale_nonpositive_or_nonfinite")
    try:
        source_label = str(path.relative_to(REPO_ROOT))
    except ValueError:
        source_label = str(path)
    return {
        "columns_read": list(projected_columns),
        "pooled_rows": int(values.size),
        "m_y": median,
        "mad_raw": raw_mad,
        "s_y": scale,
        "intercept_scale": float(2.5 * scale),
        "source_relative_path": source_label,
        "source_sha256": sha256_file(path),
    }


def load_synthetic_design(path: Path, model_id: str) -> pd.DataFrame:
    """Load one allowlisted predictor/group/key design without EEG outcomes."""

    columns = DESIGN_COLUMNS[model_id]
    frame = pq.read_table(path, columns=list(columns)).to_pandas()
    if tuple(frame.columns) != columns:
        raise RuntimeError(f"{model_id}_design_projection_changed")
    assert_synthetic_fit_input_safe(frame)
    if len(frame) != EXPECTED_ROWS[model_id]:
        raise RuntimeError(f"{model_id}_design_row_count_changed:{len(frame)}")
    if frame.isna().any().any():
        raise RuntimeError(f"{model_id}_design_contains_missing_values")
    return frame


def _verify_manifest_outputs(manifest: dict[str, Any]) -> None:
    """Hash every declared Stage 18 output and require its recorded identity."""

    for name, descriptor in manifest["outputs"].items():
        path = STAGE18_ROOT / descriptor["relative_path"]
        if not path.is_file():
            raise RuntimeError(f"stage18_output_missing:{name}:{path}")
        if path.stat().st_size != int(descriptor["size_bytes"]):
            raise RuntimeError(f"stage18_output_size_changed:{name}")
        if sha256_file(path) != descriptor["sha256"]:
            raise RuntimeError(f"stage18_output_hash_changed:{name}")


def _verify_stage18_sources(manifest: dict[str, Any]) -> None:
    """Require the accepted feature, TFR, and behavioural source fingerprints."""

    sources = manifest["sources"]
    checks = {
        REPO_ROOT / "_Data/eeg/features_v1/feature_run_manifest.json": sources[
            "feature_run_manifest_sha256"
        ],
        REPO_ROOT / "_Data/eeg/features_v1/roi_features.parquet": sources[
            "roi_features_sha256"
        ],
        REPO_ROOT / "_Data/eeg/features_v1/metadata/behavioural_lineage.parquet": sources[
            "behavioural_lineage_sha256"
        ],
        REPO_ROOT / "_Data/eeg/features_v1/validation/feature_validation.json": sources[
            "feature_validation_sha256"
        ],
        REPO_ROOT / "_Data/eeg/tfr_v1/tfr_run_manifest.json": sources[
            "tfr_run_manifest_sha256"
        ],
        REPO_ROOT / "_Private/old_rds/bdat2.rds": sources["frozen_behaviour_sha256"],
    }
    for path, expected in checks.items():
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"stage18_source_fingerprint_changed:{path}")

    immutability = json_read(STAGE18_ROOT / sources["immutability_path"])
    for descriptor in immutability["tfr_after"]:
        path = REPO_ROOT / descriptor["relative_path"]
        stat = path.stat()
        if (
            stat.st_size != int(descriptor["size_bytes"])
            or stat.st_mtime_ns != int(descriptor["mtime_ns"])
        ):
            raise RuntimeError(f"accepted_tfr_inventory_changed:{path}")


def stage18_preflight() -> dict[str, Any]:
    """Perform the complete read-only contract preflight required by Stage 19."""

    validation = json_read(STAGE18_ROOT / "validation/model_table_validation.json")
    state = json_read(STAGE18_ROOT / "model_table_run_state.json")
    manifest = json_read(STAGE18_ROOT / "model_table_run_manifest.json")
    if validation["status"] != "pass" or state["status"] != "complete":
        raise RuntimeError("stage18_not_complete_and_validated")
    if manifest["authority_fingerprint"] != STAGE18_AUTHORITY_FINGERPRINT:
        raise RuntimeError("stage18_authority_fingerprint_changed")
    if int(validation["general_rows"]) != 43_990:
        raise RuntimeError("stage18_general_row_count_changed")
    if int(validation["participants"]) != 81 or int(validation["recordings"]) != 81:
        raise RuntimeError("stage18_participant_or_recording_count_changed")
    if tuple(validation["nine_primary_estimands"]) != PRIMARY_ESTIMANDS:
        raise RuntimeError("stage18_primary_estimands_changed")
    if tuple(validation["supplementary_estimands"]) != SUPPLEMENTARY_ESTIMANDS:
        raise RuntimeError("stage18_supplementary_estimands_changed")

    _verify_manifest_outputs(manifest)
    _verify_stage18_sources(manifest)

    registry = pq.read_table(
        STAGE18_ROOT / "registries/estimand_registry.parquet",
        columns=["estimand_id", "classification", "contains_calculated_estimate"],
    ).to_pandas()
    primary = tuple(registry.loc[registry["classification"].eq("primary"), "estimand_id"])
    supplementary = tuple(
        registry.loc[registry["classification"].eq("supplementary"), "estimand_id"]
    )
    if primary != PRIMARY_ESTIMANDS or supplementary != SUPPLEMENTARY_ESTIMANDS:
        raise RuntimeError("production_estimand_registry_changed")
    if registry["contains_calculated_estimate"].any():
        raise RuntimeError("production_registry_contains_scientific_estimate")

    design_facts: dict[str, Any] = {}
    for model_id, relative in MODEL_VIEWS.items():
        path = STAGE18_ROOT / relative
        design = load_synthetic_design(path, model_id)
        fact = {
            "rows": int(len(design)),
            "participants": int(design["participant_id"].nunique()),
            "canonical_trials": int(design["canonical_event_key"].nunique()),
        }
        if model_id == "H3":
            per_trial = design.groupby("canonical_event_key", sort=False).size()
            if set(per_trial.unique()) != {2} or set(design["phase_c"].unique()) != {-0.5, 0.5}:
                raise RuntimeError("H3_two_phase_contract_changed")
            fact["rows_per_trial"] = 2
        design_facts[model_id] = fact

    return {
        "status": "pass",
        "stage18_authority_fingerprint": STAGE18_AUTHORITY_FINGERPRINT,
        "general_rows": 43_990,
        "participants": 81,
        "recordings": 81,
        "primary_estimands": list(PRIMARY_ESTIMANDS),
        "supplementary_estimands": list(SUPPLEMENTARY_ESTIMANDS),
        "design_facts": design_facts,
        "source_fingerprints_match": True,
        "accepted_outputs_opened_read_only": True,
    }


def accepted_source_snapshot() -> dict[str, Any]:
    """Snapshot accepted Stage 18 outputs and upstream authority identities."""

    manifest_path = STAGE18_ROOT / "model_table_run_manifest.json"
    manifest = json_read(manifest_path)
    files: list[dict[str, Any]] = []
    candidates = [manifest_path]
    candidates.extend(
        STAGE18_ROOT / descriptor["relative_path"]
        for descriptor in manifest["outputs"].values()
    )
    source_paths = [
        REPO_ROOT / "_Data/eeg/features_v1/feature_run_manifest.json",
        REPO_ROOT / "_Data/eeg/features_v1/roi_features.parquet",
        REPO_ROOT / "_Data/eeg/features_v1/metadata/behavioural_lineage.parquet",
        REPO_ROOT / "_Data/eeg/features_v1/validation/feature_validation.json",
        REPO_ROOT / "_Data/eeg/tfr_v1/tfr_run_manifest.json",
        REPO_ROOT / "_Private/old_rds/bdat2.rds",
    ]
    candidates.extend(source_paths)
    for path in sorted(set(candidates)):
        stat = path.stat()
        files.append(
            {
                "relative_path": str(path.relative_to(REPO_ROOT)),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": sha256_file(path),
            }
        )
    immutability = json_read(STAGE18_ROOT / "source_immutability.json")
    return {
        "stage18_authority_fingerprint": manifest["authority_fingerprint"],
        "hashed_files": files,
        "tfr_arrays": immutability["tfr_after"],
        "tfr_array_count": immutability["tfr_array_count"],
        "tfr_total_bytes": immutability["tfr_total_bytes"],
    }


def require_ignored_output_root(output_root: Path) -> None:
    """Require a narrow ignored output path inside ``_Data/eeg``."""

    resolved = output_root.resolve()
    safe_parent = (REPO_ROOT / "_Data/eeg").resolve()
    if resolved == safe_parent or safe_parent not in resolved.parents:
        raise RuntimeError(f"unsafe_model_validation_output_root:{resolved}")
    relative = resolved.relative_to(REPO_ROOT)
    check = subprocess.run(
        ["git", "check-ignore", "--quiet", str(relative)],
        cwd=REPO_ROOT,
        check=False,
    )
    if check.returncode != 0:
        raise RuntimeError(f"model_validation_output_not_ignored:{relative}")


def prepare(output_root: Path) -> None:
    """Write isolated design/scale inputs and source evidence to staging."""

    require_ignored_output_root(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    preflight = stage18_preflight()
    json_write(output_root / "validation/preflight.json", preflight)

    scales: dict[str, Any] = {
        "scale_definition": "m_y=median(value_db); s_y=1.4826*median(abs(value_db-m_y))",
        "path_boundary": "Each Arrow projection contained value_db only.",
        "models": {},
    }
    input_root = output_root / "inputs"
    input_root.mkdir(parents=True, exist_ok=True)
    for model_id, relative in MODEL_VIEWS.items():
        source = STAGE18_ROOT / relative
        scales["models"][model_id] = pooled_scale_from_value_only(source)
        design = load_synthetic_design(source, model_id)
        design.to_csv(
            input_root / f"{model_id}_accepted_design_no_outcomes.csv",
            index=False,
            float_format="%.17g",
        )
        json_write(
            input_root / f"{model_id}_design_metadata.json",
            {
                "model_id": model_id,
                "columns_read": list(DESIGN_COLUMNS[model_id]),
                "rows": int(len(design)),
                "observed_outcome_columns": [],
                "source_relative_path": relative,
                "source_sha256": sha256_file(source),
            },
        )
    json_write(output_root / "prior_scale_snapshot.json", scales)
    json_write(output_root / "source_snapshot_before.json", accepted_source_snapshot())


def verify(before_path: Path, evidence_path: Path) -> None:
    """Compare current accepted sources with ``before_path`` and write evidence."""

    stage18_preflight()
    before = json_read(before_path)
    after = accepted_source_snapshot()
    unchanged = before == after
    evidence = {
        "status": "pass" if unchanged else "fail",
        "accepted_sources_unchanged": unchanged,
        "before": before,
        "after": after,
    }
    json_write(evidence_path, evidence)
    if not unchanged:
        raise RuntimeError("accepted_stage18_or_upstream_sources_changed")


def parse_args() -> argparse.Namespace:
    """Parse the bounded prepare/verify command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--before", type=Path, required=True)
    verify_parser.add_argument("--evidence", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    """Run one authorized input-preparation or immutability-verification action."""

    args = parse_args()
    if args.command == "prepare":
        prepare(args.output_root.resolve())
    else:
        verify(args.before.resolve(), args.evidence.resolve())


if __name__ == "__main__":
    main()
