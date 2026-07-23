"""Tested contracts and rendering helpers for the participant TFR atlas.

The atlas is a bounded descriptive visual-QC product.  It reads only the
accepted stage-16 response-onset and response-end *dB* arrays, through
read-only NumPy memory maps, and averages their already-normalized trial cells
within an owner's assigned condition.  It presents all 30 physical scalp
channels using the recovered historical BESA-derived layout coordinates.

The module provides no TFR computation, baseline calculation, CSD, montage
editing, handedness mirroring, participant/trial eligibility decision,
statistic, contrast, model, or group/ROI visualisation.  BESA coordinates are
used solely as stable page-layout coordinates, not as head geometry or
individual electrode positions.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping, Sequence
import uuid

import numpy as np
import pandas as pd
import yaml


ATLAS_NAMESPACE = "participant_sensor_tfr_atlas_v1"
ATLAS_SCHEMA_VERSION = 1
EXPECTED_PARTICIPANTS = 81
EXPECTED_OVERT_PARTICIPANTS = 41
EXPECTED_IMAGERY_PARTICIPANTS = 40
EXPECTED_ASSIGNED_TRIALS = 8_060
EXPECTED_BRIDGE_EXCLUDED = 738
PRIMARY_CHANNELS = (
    "FP1", "FP2", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ",
    "FC4", "FT8", "T7", "C3", "CZ", "C4", "T8", "TP7", "CP3", "CPZ",
    "CP4", "TP8", "P7", "P3", "PZ", "P4", "P8", "O1", "OZ", "O2",
)


@dataclass(frozen=True)
class LayoutPosition:
    """One normalized historical presentation centre for a physical channel."""

    channel: str
    latitude_degrees: float
    longitude_degrees: float
    x_fraction: float
    y_fraction: float


def sha256_file(path: Path) -> str:
    """Return a SHA-256 digest without changing the requested file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_json_hash(value: Any) -> str:
    """Return a deterministic SHA-256 identity for JSON-serializable input."""

    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_and_validate_config(path: Path) -> tuple[dict[str, Any], str]:
    """Load the tracked atlas configuration and enforce its closed contract."""

    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if config.get("schema_version") != ATLAS_SCHEMA_VERSION:
        raise RuntimeError("atlas_config_schema_mismatch")
    if config.get("pipeline_version") != ATLAS_NAMESPACE:
        raise RuntimeError("atlas_config_namespace_mismatch")
    safety = config["safety"]
    required_false = (
        "rerun_convolution", "create_baseline", "apply_csd", "mirror_handedness",
        "change_trial_eligibility", "fit_models_or_statistics",
        "create_group_or_roi_plots", "write_png_previews",
    )
    if not safety.get("read_only_tfr_arrays") or not safety.get("memory_map_arrays"):
        raise RuntimeError("atlas_read_only_memory_map_contract_mismatch")
    if any(safety.get(name) is not False for name in required_false):
        raise RuntimeError("atlas_forbidden_operation_enabled")
    if tuple(safety["accepted_array_representations"]) != (
        "response_onset_db", "response_end_db"
    ):
        raise RuntimeError("atlas_array_representation_mismatch")
    accepted = config["accepted_surface"]
    expected = {
        "participant_count": EXPECTED_PARTICIPANTS,
        "overt_participant_count": EXPECTED_OVERT_PARTICIPANTS,
        "imagery_participant_count": EXPECTED_IMAGERY_PARTICIPANTS,
        "assigned_condition_trial_count": EXPECTED_ASSIGNED_TRIALS,
        "imagery_final_overt_bridge_excluded_from_means": EXPECTED_BRIDGE_EXCLUDED,
        "channel_count": len(PRIMARY_CHANNELS),
        "frequency_count": 37,
        "time_sample_count": 201,
    }
    if any(int(accepted[name]) != value for name, value in expected.items()):
        raise RuntimeError("atlas_accepted_surface_count_mismatch")
    display = config["display"]
    if (
        display["palette"] != "plasma"
        or float(display["db_minimum"]) != -8.0
        or float(display["db_maximum"]) != 8.0
        or display["output_format"] != "pdf"
    ):
        raise RuntimeError("atlas_display_contract_mismatch")
    return config, sha256_file(path)


def load_layout_positions(path: Path, channels: Sequence[str]) -> dict[str, LayoutPosition]:
    """Load historical BESA layout coordinates without treating them as montage data."""

    table = pd.read_csv(path)
    required = {"lat", "long", "chan"}
    if set(table.columns) != required:
        raise RuntimeError("historical_layout_columns_mismatch")
    table = table.copy()
    table["channel"] = table["chan"].astype(str).str.upper()
    expected = list(channels)
    if len(table) != len(expected) or set(table["channel"]) != set(expected):
        raise RuntimeError("historical_layout_channel_set_mismatch")
    radians = np.deg2rad(table["long"].to_numpy(dtype=float))
    # This is the exact historical radial projection recovered in the audit.
    x = table["lat"].to_numpy(dtype=float) * np.cos(radians)
    y = table["lat"].to_numpy(dtype=float) * np.sin(radians)
    x_fraction = 0.06 + 0.88 * (x - x.min()) / (x.max() - x.min())
    y_fraction = 0.10 + 0.78 * (y - y.min()) / (y.max() - y.min())
    result = {
        row.channel: LayoutPosition(
            channel=row.channel,
            latitude_degrees=float(row.lat),
            longitude_degrees=float(row.long),
            x_fraction=float(x_fraction[index]),
            y_fraction=float(y_fraction[index]),
        )
        for index, row in enumerate(table.itertuples(index=False))
    }
    return {channel: result[channel] for channel in expected}


def select_assigned_condition_trials(lineage: pd.DataFrame, config: Mapping[str, Any]) -> pd.DataFrame:
    """Select display trials by assigned participant group, retaining bridge provenance.

    Physical-group participants retain their accepted overt-movement trials.
    Imagery-group participants retain their accepted motor-imagery trials only.
    The latter's accepted final-overt bridge rows are marked and excluded from
    the atlas mean, never removed from the accepted source surface.
    """

    required = {
        "behavioural_id", "canonical_event_key", "recording_stem", "tfr_row_index",
        "group", "performed_condition", "condition_semantics", "physical",
        "scope_all_accepted", "scope_imagery_final_overt_bridge",
    }
    missing = required.difference(lineage.columns)
    if missing:
        raise RuntimeError(f"atlas_lineage_columns_missing:{sorted(missing)}")
    if not lineage["scope_all_accepted"].all() or lineage["canonical_event_key"].duplicated().any():
        raise RuntimeError("atlas_lineage_accepted_identity_mismatch")
    selection = config["selection"]
    overt = lineage["group"].eq(selection["overt_group"])
    imagery = lineage["group"].eq(selection["imagery_group"])
    if not (overt | imagery).all():
        raise RuntimeError("atlas_unknown_participant_group")
    selected = (overt & lineage["condition_semantics"].eq(selection["overt_condition_semantics"])) | (
        imagery & lineage["condition_semantics"].eq(selection["imagery_condition_semantics"])
    )
    bridge = imagery & lineage["scope_imagery_final_overt_bridge"].astype(bool)
    if int(bridge.sum()) != EXPECTED_BRIDGE_EXCLUDED or bool(selected.loc[bridge].any()):
        raise RuntimeError("atlas_bridge_selection_mismatch")
    result = lineage.copy()
    result["atlas_pdf_group"] = np.where(overt, "overt", "imagery")
    result["atlas_assigned_condition"] = np.where(overt, "overt_movement", "motor_imagery")
    result["atlas_selected_for_mean"] = selected
    result["atlas_bridge_excluded_from_mean"] = bridge
    return result


def participant_page_table(selection: pd.DataFrame) -> pd.DataFrame:
    """Create a deterministic participant-to-PDF/page table and validate counts."""

    selected = selection.loc[selection["atlas_selected_for_mean"]].copy()
    grouped = selected.groupby("behavioural_id", sort=True, observed=True)
    rows: list[dict[str, Any]] = []
    for pdf_group, part in selected.groupby("atlas_pdf_group", sort=True, observed=True):
        for page, (participant_id, frame) in enumerate(part.groupby("behavioural_id", sort=True), start=1):
            if frame["recording_stem"].nunique() != 1:
                raise RuntimeError(f"atlas_participant_recording_mismatch:{participant_id}")
            rows.append(
                {
                    "behavioural_id": int(participant_id), "pdf_group": str(pdf_group),
                    "page_number": page, "recording_stem": frame["recording_stem"].iloc[0],
                    "assigned_condition": frame["atlas_assigned_condition"].iloc[0],
                    "contributing_trial_count": int(len(frame)),
                }
            )
    table = pd.DataFrame(rows).sort_values(["pdf_group", "page_number"]).reset_index(drop=True)
    if len(table) != EXPECTED_PARTICIPANTS or table["behavioural_id"].duplicated().any():
        raise RuntimeError("atlas_participant_page_identity_mismatch")
    if len(grouped) != EXPECTED_PARTICIPANTS or int(selected.shape[0]) != EXPECTED_ASSIGNED_TRIALS:
        raise RuntimeError("atlas_assigned_trial_total_mismatch")
    if (table["pdf_group"].eq("overt")).sum() != EXPECTED_OVERT_PARTICIPANTS:
        raise RuntimeError("atlas_overt_page_count_mismatch")
    if (table["pdf_group"].eq("imagery")).sum() != EXPECTED_IMAGERY_PARTICIPANTS:
        raise RuntimeError("atlas_imagery_page_count_mismatch")
    return table


def mean_memory_mapped_trials(array_path: Path, row_indices: Sequence[int], expected_shape: tuple[int, int, int]) -> np.ndarray:
    """Average selected trial cells one at a time from a read-only NumPy memory map."""

    values = np.load(array_path, mmap_mode="r")
    if values.ndim != 4 or tuple(values.shape[1:]) != expected_shape:
        raise RuntimeError(f"atlas_array_shape_mismatch:{array_path}:{values.shape}")
    rows = np.asarray(row_indices, dtype=int)
    if rows.size == 0 or rows.min() < 0 or rows.max() >= values.shape[0]:
        raise RuntimeError(f"atlas_row_indices_invalid:{array_path}")
    accumulator = np.zeros(expected_shape, dtype=np.float64)
    # Avoid fancy indexing, which would duplicate the selected trial cube.
    for row in rows:
        accumulator += values[int(row)]
    return accumulator / float(rows.size)


def saturation_accounting(onset: np.ndarray, end: np.ndarray, minimum: float, maximum: float) -> dict[str, Any]:
    """Count out-of-scale rendered mean cells before image clipping is applied."""

    if onset.shape != end.shape or onset.ndim != 3:
        raise RuntimeError("atlas_rendered_surface_shape_mismatch")
    values = np.concatenate((onset.ravel(), end.ravel()))
    below = int(np.count_nonzero(values < minimum))
    above = int(np.count_nonzero(values > maximum))
    count = int(values.size)
    return {
        "rendered_source_cell_count": count,
        "below_minimum_count": below,
        "above_maximum_count": above,
        "saturated_cell_count": below + above,
        "saturated_percent": 100.0 * (below + above) / count,
        "minimum_rendered_value_db": float(values.min()),
        "maximum_rendered_value_db": float(values.max()),
        "display_minimum_db": float(minimum),
        "display_maximum_db": float(maximum),
    }


def source_snapshot(paths: Iterable[Path], *, root: Path) -> list[dict[str, Any]]:
    """Capture source path/size/mtime evidence for a before/after immutability check."""

    unique = sorted({Path(path) for path in paths}, key=lambda item: item.as_posix())
    return [
        {
            "relative_path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for path in unique
    ]


def compare_source_snapshots(before: Sequence[Mapping[str, Any]], after: Sequence[Mapping[str, Any]]) -> None:
    """Fail closed if an atlas source path, size, or modification time changes."""

    if list(before) != list(after):
        raise RuntimeError("atlas_source_size_or_mtime_changed")


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Write and parse-check JSON beside its destination before atomic replacement."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        json.loads(temporary.read_text(encoding="utf-8"))
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    """Write a CSV beside its destination, reopen it, then atomically publish it."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        frame.to_csv(temporary, index=False)
        if len(pd.read_csv(temporary)) != len(frame):
            raise RuntimeError(f"atlas_csv_reopen_mismatch:{path}")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_text(path: Path, text: str) -> None:
    """Write UTF-8 text atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def pdf_page_count(path: Path) -> int:
    """Reopen a PDF with the system PDF inspector and return its page count."""

    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"atlas_pdf_page_count_not_found:{path}")


def output_descriptor(path: Path, *, root: Path) -> dict[str, Any]:
    """Return a content descriptor for one atlas artifact."""

    return {
        "relative_path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }

