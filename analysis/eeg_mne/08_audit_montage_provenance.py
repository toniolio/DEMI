"""Audit DEMI historical montage provenance and compare MNE templates.

This script validates the tracked BESA unit-sphere table against
``montage_coordinate_contract_v1.yaml`` and performs a documented comparison
of the 32 DEMI scalp labels with MNE ``standard_1005``. It exists to answer
coordinate provenance/frame/unit/fiducial questions and records the later
approved active MNE montage decision.

Inputs:

* ``_Data/eeg/BESA-81.csv``;
* ``analysis/eeg_mne/montage_coordinate_contract_v1.yaml``;
* the installed MNE ``standard_1005`` template.

Local outputs below
``_Data/eeg/mne_preprocessing/channel_qc_v1/montage_audit``:

* ``historical_coordinate_inventory.csv``;
* ``besa_standard_1005_comparison.csv``;
* ``montage_audit_summary.json``;
* ``montage_audit_summary.md``;
* ``montage_audit_run_manifest.json``.

The comparison first constructs a BESA-template montage with the table's
NAS/LPA/RPA fiducials, converts it to an MNE head-oriented copy, and then fits
one uniform scale to ``standard_1005`` across the 32 labels. The fitted scale
is a comparison device only. It is not a recovered head radius, acquisition
digitization, or approved conversion for interpolation/CSD.

Safety boundaries:

* no raw EDF is opened;
* no montage is applied to EEG data;
* no interpolation, topographic signal calculation, CSD, filtering,
  referencing, ICA, bad-channel assignment, or epoch construction occurs;
* the active montage decision is reported but no montage is applied to signal.

Run from the repository root:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/08_audit_montage_provenance.py
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

# MNE can initialize matplotlib during import. Keep caches out of the user
# configuration and repository even though this audit writes no figures.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "demi_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "demi_cache"))

import mne
import numpy as np
import pandas as pd
from mne.transforms import apply_trans

from montage_contract import (
    load_coordinate_table,
    load_montage_contract,
    normalize_channel_name,
    sha256_file,
    validate_coordinate_table,
)


CONTRACT_PATH = Path("analysis") / "eeg_mne" / "montage_coordinate_contract_v1.yaml"
COORDINATE_PATH = Path("_Data") / "eeg" / "BESA-81.csv"
OUTPUT_DIR = (
    Path("_Data")
    / "eeg"
    / "mne_preprocessing"
    / "channel_qc_v1"
    / "montage_audit"
)

INVENTORY_FILENAME = "historical_coordinate_inventory.csv"
COMPARISON_FILENAME = "besa_standard_1005_comparison.csv"
SUMMARY_JSON_FILENAME = "montage_audit_summary.json"
SUMMARY_MARKDOWN_FILENAME = "montage_audit_summary.md"
RUN_MANIFEST_FILENAME = "montage_audit_run_manifest.json"

EXPECTED_EEG_CHANNELS = (
    "FP1", "FP2", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ",
    "FC4", "FT8", "T7", "C3", "CZ", "C4", "T8", "M1", "TP7", "CP3",
    "CPZ", "CP4", "TP8", "M2", "P7", "P3", "PZ", "P4", "P8", "O1",
    "OZ", "O2",
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script path."""

    return Path(__file__).resolve().parents[2]


def coordinate_role(normalized_name: str, contract: dict[str, Any]) -> str:
    """Classify one historical coordinate row for the audit inventory."""

    fiducials = {
        normalize_channel_name(value)
        for value in contract["geometry"]["fiducials"].values()
    }
    aliases = {
        normalize_channel_name(value)
        for value in contract["geometry"].get("fiducial_aliases", {}).values()
    }
    if normalized_name in fiducials:
        return "template_fiducial"
    if normalized_name in aliases:
        return "template_fiducial_alias"
    if normalized_name in EXPECTED_EEG_CHANNELS:
        return "demi_expected_eeg"
    return "other_template_location"


def transformed_standard_1005_positions() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Return standard_1005 channel and fiducial positions in MNE head frame."""

    montage = mne.channels.make_standard_montage("standard_1005")
    positions = montage.get_positions()
    transform = mne.channels.compute_native_head_t(montage)
    channels = {
        normalize_channel_name(name): apply_trans(transform, xyz)
        for name, xyz in positions["ch_pos"].items()
    }
    fiducials = {
        name: apply_trans(transform, positions[name])
        for name in ("nasion", "lpa", "rpa")
    }
    return channels, fiducials


def transformed_besa_positions(
    coordinates: pd.DataFrame,
    contract: dict[str, Any],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Return an audit-only head-oriented copy of unit-sphere coordinates."""

    by_name = {
        row.normalized_channel_name: np.array([row.x, row.y, row.z], dtype=float)
        for row in coordinates.itertuples(index=False)
    }
    fiducial_labels = contract["geometry"]["fiducials"]
    montage = mne.channels.make_dig_montage(
        ch_pos={name: by_name[name] for name in EXPECTED_EEG_CHANNELS},
        nasion=by_name[normalize_channel_name(fiducial_labels["nasion"])],
        lpa=by_name[normalize_channel_name(fiducial_labels["left_preauricular"])],
        rpa=by_name[normalize_channel_name(fiducial_labels["right_preauricular"])],
        coord_frame="unknown",
    )
    transform = mne.channels.compute_native_head_t(montage)
    channels = {
        name: apply_trans(transform, by_name[name]) for name in EXPECTED_EEG_CHANNELS
    }
    fiducials = {
        "nasion": apply_trans(
            transform, by_name[normalize_channel_name(fiducial_labels["nasion"])]
        ),
        "lpa": apply_trans(
            transform,
            by_name[normalize_channel_name(fiducial_labels["left_preauricular"])],
        ),
        "rpa": apply_trans(
            transform,
            by_name[normalize_channel_name(fiducial_labels["right_preauricular"])],
        ),
    }
    return channels, fiducials


def build_comparison(
    coordinates: pd.DataFrame,
    contract: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Compare head-oriented unit-sphere and standard_1005 coordinates."""

    besa, besa_fiducials = transformed_besa_positions(coordinates, contract)
    standard, standard_fiducials = transformed_standard_1005_positions()
    missing = sorted(set(EXPECTED_EEG_CHANNELS).difference(standard))
    if missing:
        raise RuntimeError(f"standard_1005 missing expected DEMI labels: {missing}")

    besa_matrix = np.array([besa[name] for name in EXPECTED_EEG_CHANNELS])
    standard_matrix = np.array([standard[name] for name in EXPECTED_EEG_CHANNELS])
    # One least-squares uniform scale makes the coordinate-shape differences
    # inspectable in millimetres. No acquisition head size is inferred.
    fitted_scale = float(
        np.sum(besa_matrix * standard_matrix) / np.sum(besa_matrix * besa_matrix)
    )
    scaled_besa = besa_matrix * fitted_scale
    distances_mm = np.linalg.norm(scaled_besa - standard_matrix, axis=1) * 1000.0

    rows = []
    for index, name in enumerate(EXPECTED_EEG_CHANNELS):
        rows.append(
            {
                "normalized_channel_name": name,
                "besa_head_oriented_unit_x": besa_matrix[index, 0],
                "besa_head_oriented_unit_y": besa_matrix[index, 1],
                "besa_head_oriented_unit_z": besa_matrix[index, 2],
                "comparison_uniform_scale_m_per_unit": fitted_scale,
                "besa_scaled_comparison_x_m": scaled_besa[index, 0],
                "besa_scaled_comparison_y_m": scaled_besa[index, 1],
                "besa_scaled_comparison_z_m": scaled_besa[index, 2],
                "standard_1005_head_x_m": standard_matrix[index, 0],
                "standard_1005_head_y_m": standard_matrix[index, 1],
                "standard_1005_head_z_m": standard_matrix[index, 2],
                "euclidean_distance_mm": distances_mm[index],
                "comparison_only_not_approved_montage_conversion": True,
            }
        )

    fiducial_distances = {
        name: float(
            np.linalg.norm(
                besa_fiducials[name] * fitted_scale - standard_fiducials[name]
            )
            * 1000.0
        )
        for name in ("nasion", "lpa", "rpa")
    }
    summary = {
        "comparison_label_count": len(EXPECTED_EEG_CHANNELS),
        "comparison_uniform_scale_m_per_unit": fitted_scale,
        "comparison_uniform_scale_is_not_recovered_head_radius": True,
        "distance_mm_mean": float(np.mean(distances_mm)),
        "distance_mm_median": float(np.median(distances_mm)),
        "distance_mm_minimum": float(np.min(distances_mm)),
        "distance_mm_maximum": float(np.max(distances_mm)),
        "distance_mm_maximum_channel": EXPECTED_EEG_CHANNELS[int(np.argmax(distances_mm))],
        "template_fiducial_distance_mm_after_same_scale": fiducial_distances,
    }
    return pd.DataFrame(rows), summary


def git_text(repo_root: Path, arguments: list[str]) -> str:
    """Return Git provenance without changing repository state."""

    result = subprocess.run(
        ["git", *arguments],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def write_markdown(path: Path, summary: dict[str, Any]) -> None:
    """Write concise direct answers and comparison results."""

    comparison = summary["standard_1005_comparison"]
    lines = [
        "# Montage provenance audit",
        "",
        "## Direct evidence",
        "",
        "- The tracked 112-row CSV exactly matches Robert Oostenveld's published `besa_81.txt` labels and six-decimal Cartesian coordinates.",
        "- The source describes these as coordinates according to BESA on a unit sphere.",
        "- Orientation is x right, y anterior/nasion, z superior/Cz.",
        "- NAS, LPA, and RPA template fiducials are present; Nz, T9, and T10 duplicate them as aliases.",
        "- Historical DEMI R code used this surface for spherical plotting/spatial coordinates. Austin's historical Python preprocessing used MNE `standard_1005` instead.",
        "",
        "## Limits",
        "",
        "No participant electrode digitization, participant head shape, physical head radius, acquisition-specific transform, or documented cap-to-template measurement was found. The template fiducials orient a unit sphere but do not recover individual acquisition geometry.",
        "",
        "## standard_1005 comparison",
        "",
        f"All {comparison['comparison_label_count']} DEMI scalp labels are present in `standard_1005`. After an audit-only fiducial orientation and one best-fitting uniform scale ({comparison['comparison_uniform_scale_m_per_unit']:.6f} m/unit), coordinate differences have mean {comparison['distance_mm_mean']:.2f} mm, median {comparison['distance_mm_median']:.2f} mm, and maximum {comparison['distance_mm_maximum']:.2f} mm at {comparison['distance_mm_maximum_channel']}. The fitted scale is not an inferred DEMI head radius.",
        "",
        "## Use boundary",
        "",
        "- Historical/template visualization and channel-name validation: supported.",
        "- Current quantitative topography/interpolation: use the approved MNE `standard_1005` template.",
        "- Interpolation and CSD: the unitless historical table is not approved as-is.",
        "- Active montage: `standard_1005`, approved 2026-07-12. The historical BESA template remains historical GAM/visualization provenance only; CSD is deferred.",
        "",
        "No EEG data were opened or modified by this audit.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Validate provenance, compare coordinate surfaces, and write local evidence."""

    repo_root = repo_root_from_script()
    contract_path = repo_root / CONTRACT_PATH
    coordinate_path = repo_root / COORDINATE_PATH
    output_dir = repo_root / OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    contract = load_montage_contract(contract_path)
    coordinates = load_coordinate_table(coordinate_path)
    validate_coordinate_table(coordinates, contract, coordinate_path)
    inventory = coordinates.copy()
    inventory["coordinate_role"] = inventory["normalized_channel_name"].map(
        lambda name: coordinate_role(name, contract)
    )
    inventory["coordinate_frame"] = contract["geometry"]["coordinate_frame"]
    inventory["coordinate_units"] = contract["geometry"]["units"]
    inventory["acquisition_digitization"] = False
    inventory["safe_as_metric_head_coordinates"] = False

    comparison, comparison_summary = build_comparison(coordinates, contract)
    source = contract["source"]
    summary = {
        "contract_version": contract["contract_version"],
        "coordinate_source_name": source["name"],
        "coordinate_source_url": source["url"],
        "reviewed_source_match_status": source["reviewed_match_status"],
        "coordinate_row_count": len(coordinates),
        "coordinate_norm_minimum": float(coordinates["coordinate_norm"].min()),
        "coordinate_norm_maximum": float(coordinates["coordinate_norm"].max()),
        "coordinate_frame": contract["geometry"]["coordinate_frame"],
        "coordinate_units": contract["geometry"]["units"],
        "orientation": contract["geometry"]["orientation"],
        "fiducials": contract["geometry"]["fiducials"],
        "fiducial_aliases": contract["geometry"]["fiducial_aliases"],
        "individual_digitization_available": False,
        "individual_head_shape_available": False,
        "physical_head_radius_available": False,
        "standard_1005_comparison": comparison_summary,
        "active_montage": contract["active_pipeline"]["active_montage"],
        "active_montage_decision_status": contract["active_pipeline"]["decision_status"],
        "active_montage_decision_date": str(contract["active_pipeline"]["decision_date"]),
        "historical_besa_role": contract["active_pipeline"]["historical_besa_role"],
        "safe_use": {
            "historical_spherical_visualization": True,
            "channel_name_coverage": True,
            "current_quantitative_topography_with_standard_1005": True,
            "interpolation_as_is": False,
            "csd_as_is": False,
        },
    }

    inventory.to_csv(output_dir / INVENTORY_FILENAME, index=False)
    comparison.to_csv(output_dir / COMPARISON_FILENAME, index=False)
    (output_dir / SUMMARY_JSON_FILENAME).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_markdown(output_dir / SUMMARY_MARKDOWN_FILENAME, summary)

    script_path = Path(__file__).resolve()
    helper_path = script_path.parent / "montage_contract.py"
    manifest = {
        "run_timestamp_local": datetime.now().astimezone().isoformat(),
        "script_path": script_path.relative_to(repo_root).as_posix(),
        "script_sha256": sha256_file(script_path),
        "helper_sha256": sha256_file(helper_path),
        "contract_sha256": sha256_file(contract_path),
        "coordinate_file_sha256": sha256_file(coordinate_path),
        "git_commit": git_text(repo_root, ["rev-parse", "HEAD"]),
        "git_status_short_at_run": git_text(repo_root, ["status", "--short"]),
        "python_version": platform.python_version(),
        "mne_version": mne.__version__,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "outputs": sorted(path.name for path in output_dir.iterdir() if path.is_file()),
        "safety": {
            "raw_edf_opened": False,
            "montage_applied_to_signal": False,
            "interpolation": False,
            "csd": False,
            "active_montage_selected": True,
            "active_montage_applied_to_signal": False,
        },
    }
    (output_dir / RUN_MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Validated {len(inventory)} historical coordinate rows and compared "
        f"{len(comparison)} DEMI channels in {output_dir.relative_to(repo_root)}"
    )


if __name__ == "__main__":
    main()
