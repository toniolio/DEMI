"""Pure provenance and geometry guards for DEMI montage coordinates.

The tracked ``_Data/eeg/BESA-81.csv`` table is an exact copy of Robert
Oostenveld's published ``besa_81.txt`` unit-sphere template. That establishes
its upstream source, row meanings, unitless geometry, orientation, and
template fiducials. It does not establish that DEMI electrodes were digitized
at those locations, a participant head radius, or an acquisition-specific
transform.

This module validates the tracked sidecar contract and coordinate table. It
also provides an explicit guard that refuses to treat the unit-sphere table as
metric MNE head coordinates. Comparison scripts may transform and scale a copy
for documented audit purposes, but active interpolation or CSD code must not
silently pass this guard.

The module performs no raw-data I/O, signal preprocessing, interpolation,
topographic calculation, CSD, or file writing. Public tests use only synthetic
tables and dictionaries.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml


CONTRACT_VERSION = "montage_coordinate_contract_v1"
REQUIRED_TOP_LEVEL_KEYS = {
    "contract_version",
    "coordinate_file",
    "coordinate_file_sha256",
    "canonical_coordinate_rows_sha256",
    "expected_row_count",
    "unit_norm_tolerance",
    "source",
    "geometry",
    "historical_use",
    "active_pipeline",
}
REQUIRED_SOURCE_KEYS = {
    "name",
    "url",
    "source_page",
    "reviewed_match_status",
    "reviewed_match_date",
    "source_description",
    "acquisition_digitization_source",
}
REQUIRED_GEOMETRY_KEYS = {
    "coordinate_frame",
    "units",
    "origin",
    "orientation",
    "fiducials",
    "head_radius_m",
    "individual_head_shape_available",
    "individual_electrode_digitization_available",
}
EXPECTED_ORIENTATION = {
    "x_positive": "right",
    "y_positive": "anterior_toward_nasion",
    "z_positive": "superior_toward_cz",
}
EXPECTED_FIDUCIAL_KEYS = {
    "nasion",
    "left_preauricular",
    "right_preauricular",
}


def normalize_channel_name(value: object) -> str:
    """Normalize a channel label by trimming and uppercasing only."""

    return str(value or "").strip().upper()


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a coordinate or contract file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_coordinate_rows_sha256(coordinates: pd.DataFrame) -> str:
    """Hash ordered labels and coordinates rounded to the source precision."""

    required = {"chan", "x", "y", "z"}
    missing = required.difference(coordinates.columns)
    if missing:
        raise ValueError(f"coordinate table missing columns: {sorted(missing)}")
    payload = "".join(
        f"{row.chan},{float(row.x):.6f},{float(row.y):.6f},{float(row.z):.6f}\n"
        for row in coordinates.itertuples(index=False)
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_montage_contract(contract: Mapping[str, Any]) -> None:
    """Fail on missing or internally inconsistent required metadata.

    Args:
        contract: Parsed montage contract mapping.

    Returns:
        None when every required provenance and geometry field is coherent.

    Side effects:
        None.
    """

    missing = REQUIRED_TOP_LEVEL_KEYS.difference(contract)
    if missing:
        raise ValueError(f"montage contract missing required keys: {sorted(missing)}")
    if contract["contract_version"] != CONTRACT_VERSION:
        raise ValueError("unexpected montage contract version")
    source = contract["source"]
    geometry = contract["geometry"]
    active = contract["active_pipeline"]
    if not isinstance(source, Mapping):
        raise ValueError("montage source metadata must be a mapping")
    if not isinstance(geometry, Mapping):
        raise ValueError("montage geometry metadata must be a mapping")
    if not isinstance(active, Mapping):
        raise ValueError("active_pipeline metadata must be a mapping")
    source_missing = REQUIRED_SOURCE_KEYS.difference(source)
    geometry_missing = REQUIRED_GEOMETRY_KEYS.difference(geometry)
    if source_missing:
        raise ValueError(f"montage source missing required keys: {sorted(source_missing)}")
    if geometry_missing:
        raise ValueError(f"montage geometry missing required keys: {sorted(geometry_missing)}")
    if bool(source["acquisition_digitization_source"]):
        raise ValueError("historical unit-sphere source cannot claim acquisition digitization")
    if geometry["coordinate_frame"] != "besa_unit_sphere_template":
        raise ValueError("unexpected or unsupported historical coordinate frame")
    if geometry["units"] != "unitless":
        raise ValueError("BESA template coordinates must remain unitless")
    if geometry["orientation"] != EXPECTED_ORIENTATION:
        raise ValueError("montage orientation metadata is missing or inconsistent")
    fiducials = geometry["fiducials"]
    if not isinstance(fiducials, Mapping) or set(fiducials) != EXPECTED_FIDUCIAL_KEYS:
        raise ValueError("montage fiducial metadata is missing or inconsistent")
    if geometry["head_radius_m"] is not None:
        raise ValueError("historical unit-sphere contract must not invent head_radius_m")
    if bool(geometry["individual_head_shape_available"]):
        raise ValueError("historical template must not claim individual head shape")
    if bool(geometry["individual_electrode_digitization_available"]):
        raise ValueError("historical template must not claim individual digitization")
    if active.get("active_montage") != "standard_1005":
        raise ValueError("approved active montage must be standard_1005")
    if active.get("decision_status") != "approved":
        raise ValueError("active montage decision status is inconsistent")
    if active.get("historical_besa_role") != "historical_gam_and_visualization_provenance_only":
        raise ValueError("historical BESA role is inconsistent with the approved decision")
    if active.get("csd_status") != "deferred":
        raise ValueError("CSD status is inconsistent with the approved decision")
    if int(contract["expected_row_count"]) <= 0:
        raise ValueError("expected_row_count must be positive")
    if float(contract["unit_norm_tolerance"]) <= 0:
        raise ValueError("unit_norm_tolerance must be positive")


def load_montage_contract(path: Path) -> dict[str, Any]:
    """Read and validate the tracked YAML montage contract."""

    if not path.is_file():
        raise RuntimeError(f"missing montage contract: {path}")
    parsed = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError("montage contract must parse to a mapping")
    validate_montage_contract(parsed)
    return parsed


def load_coordinate_table(path: Path) -> pd.DataFrame:
    """Read a coordinate CSV and add strict normalized labels and norms."""

    if not path.is_file():
        raise RuntimeError(f"missing coordinate table: {path}")
    coordinates = pd.read_csv(path)
    required = {"chan", "x", "y", "z"}
    missing = required.difference(coordinates.columns)
    if missing:
        raise ValueError(f"coordinate table missing columns: {sorted(missing)}")
    coordinates = coordinates.copy()
    coordinates["normalized_channel_name"] = coordinates["chan"].map(normalize_channel_name)
    if coordinates["normalized_channel_name"].eq("").any():
        raise ValueError("coordinate table contains a blank channel label")
    duplicates = coordinates.loc[
        coordinates["normalized_channel_name"].duplicated(keep=False),
        "normalized_channel_name",
    ].unique()
    if len(duplicates):
        raise ValueError(f"coordinate table has duplicate normalized labels: {sorted(duplicates)}")
    for column in ("x", "y", "z"):
        coordinates[column] = pd.to_numeric(coordinates[column], errors="raise")
    xyz = coordinates[["x", "y", "z"]].to_numpy(dtype=float)
    if not np.all(np.isfinite(xyz)):
        raise ValueError("coordinate table contains non-finite coordinates")
    coordinates["coordinate_norm"] = np.linalg.norm(xyz, axis=1)
    return coordinates


def validate_coordinate_table(
    coordinates: pd.DataFrame,
    contract: Mapping[str, Any],
    coordinate_path: Path | None = None,
) -> None:
    """Validate table shape, unit sphere, fiducials, and optional file hashes."""

    validate_montage_contract(contract)
    if len(coordinates) != int(contract["expected_row_count"]):
        raise ValueError(
            f"coordinate row count {len(coordinates)} does not match contract "
            f"{contract['expected_row_count']}"
        )
    tolerance = float(contract["unit_norm_tolerance"])
    if not np.allclose(coordinates["coordinate_norm"], 1.0, atol=tolerance, rtol=0.0):
        raise ValueError("coordinates are inconsistent with the declared unit sphere")
    labels = set(coordinates["normalized_channel_name"])
    for role, label in contract["geometry"]["fiducials"].items():
        if normalize_channel_name(label) not in labels:
            raise ValueError(f"required {role} fiducial {label!r} is missing")
    canonical_digest = canonical_coordinate_rows_sha256(coordinates)
    if canonical_digest != contract["canonical_coordinate_rows_sha256"]:
        raise ValueError("canonical coordinate-row digest does not match contract")
    if coordinate_path is not None:
        if sha256_file(coordinate_path) != contract["coordinate_file_sha256"]:
            raise ValueError("coordinate file SHA-256 does not match contract")


def require_metric_head_coordinates(contract: Mapping[str, Any]) -> None:
    """Guard active MNE spatial operations against unitless template geometry.

    This function succeeds only for an explicit metric MNE head-coordinate
    contract with a positive head radius. The reviewed historical BESA
    contract is expected to fail. That failure is the safety behavior: callers
    must make an explicit conversion/strategy decision rather than silently
    passing unitless values to interpolation or CSD.
    """

    geometry = contract.get("geometry")
    if not isinstance(geometry, Mapping):
        raise ValueError("missing geometry metadata")
    if geometry.get("coordinate_frame") != "head":
        raise ValueError("active spatial operation requires MNE head coordinates")
    if geometry.get("units") != "m":
        raise ValueError("active spatial operation requires coordinates in meters")
    head_radius = geometry.get("head_radius_m")
    if head_radius is None or float(head_radius) <= 0:
        raise ValueError("active spatial operation requires a positive head_radius_m")
