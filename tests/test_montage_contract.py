"""Synthetic tests for DEMI montage provenance and coordinate guards."""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_ANALYSIS_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_ANALYSIS_DIR))

from montage_contract import (  # noqa: E402
    canonical_coordinate_rows_sha256,
    load_coordinate_table,
    load_montage_contract,
    normalize_channel_name,
    require_metric_head_coordinates,
    validate_coordinate_table,
    validate_montage_contract,
)


CONTRACT_PATH = EEG_ANALYSIS_DIR / "montage_coordinate_contract_v1.yaml"
COORDINATE_PATH = REPO_ROOT / "_Data" / "eeg" / "BESA-81.csv"


def valid_contract() -> dict:
    """Return a mutable copy of the tracked reviewed contract."""

    return copy.deepcopy(load_montage_contract(CONTRACT_PATH))


def test_name_normalization_is_case_only_and_stable() -> None:
    """Expected EDF/BESA case variants normalize without speculative aliases."""

    assert normalize_channel_name(" Fp1 ") == "FP1"
    assert normalize_channel_name("FCz") == "FCZ"
    assert normalize_channel_name("EMG-L") == "EMG-L"


def test_tracked_coordinate_table_and_contract_validate() -> None:
    """The tracked 112-row unit-sphere source is hash- and metadata-guarded."""

    contract = valid_contract()
    coordinates = load_coordinate_table(COORDINATE_PATH)
    validate_coordinate_table(coordinates, contract, COORDINATE_PATH)
    assert len(coordinates) == 112
    assert coordinates["coordinate_norm"].between(0.99999, 1.00001).all()
    assert canonical_coordinate_rows_sha256(coordinates) == (
        contract["canonical_coordinate_rows_sha256"]
    )


@pytest.mark.parametrize(
    "missing_key",
    ["source", "geometry", "coordinate_file_sha256", "active_pipeline"],
)
def test_missing_required_contract_metadata_fails(missing_key: str) -> None:
    """Provenance/frame/unit metadata cannot be silently omitted."""

    contract = valid_contract()
    del contract[missing_key]
    with pytest.raises(ValueError, match="missing|required"):
        validate_montage_contract(contract)


def test_inconsistent_orientation_and_fiducials_fail() -> None:
    """Axis direction and all three template fiducials are required."""

    orientation = valid_contract()
    orientation["geometry"]["orientation"]["x_positive"] = "left"
    with pytest.raises(ValueError, match="orientation"):
        validate_montage_contract(orientation)

    fiducials = valid_contract()
    del fiducials["geometry"]["fiducials"]["nasion"]
    with pytest.raises(ValueError, match="fiducial"):
        validate_montage_contract(fiducials)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("coordinate_frame", "head", "coordinate frame"),
        ("units", "m", "unitless"),
        ("head_radius_m", 0.095, "head_radius"),
    ],
)
def test_historical_frame_unit_and_scale_cannot_be_relabelled(
    field: str,
    value: object,
    message: str,
) -> None:
    """The provenance contract refuses an unsupported metric reinterpretation."""

    contract = valid_contract()
    contract["geometry"][field] = value
    with pytest.raises(ValueError, match=message):
        validate_montage_contract(contract)


def test_unit_sphere_contract_fails_metric_head_coordinate_guard() -> None:
    """Historical coordinates cannot silently reach interpolation or CSD."""

    with pytest.raises(ValueError, match="MNE head coordinates"):
        require_metric_head_coordinates(valid_contract())


def test_explicit_metric_head_contract_can_pass_guard() -> None:
    """The guard accepts only explicit head-frame meters and a positive scale."""

    contract = valid_contract()
    contract["geometry"]["coordinate_frame"] = "head"
    contract["geometry"]["units"] = "m"
    contract["geometry"]["head_radius_m"] = 0.095
    require_metric_head_coordinates(contract)


def test_coordinate_parser_surfaces_duplicate_normalized_names(tmp_path: Path) -> None:
    """Case variants that collapse to one channel name are rejected."""

    path = tmp_path / "coordinates.csv"
    pd.DataFrame(
        [
            {"chan": "Fp1", "x": 1, "y": 0, "z": 0},
            {"chan": "FP1", "x": 0, "y": 1, "z": 0},
        ]
    ).to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate normalized labels"):
        load_coordinate_table(path)


def test_missing_required_fiducial_in_coordinates_fails() -> None:
    """A coordinate surface cannot validate with an absent fiducial row."""

    contract = valid_contract()
    coordinates = load_coordinate_table(COORDINATE_PATH)
    coordinates = coordinates[coordinates["normalized_channel_name"].ne("NAS")].copy()
    contract["expected_row_count"] = len(coordinates)
    contract["canonical_coordinate_rows_sha256"] = canonical_coordinate_rows_sha256(
        coordinates
    )
    with pytest.raises(ValueError, match="nasion fiducial"):
        validate_coordinate_table(coordinates, contract)
