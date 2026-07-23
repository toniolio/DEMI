"""Focused contract tests for the read-only participant sensor TFR atlas."""

from __future__ import annotations

from pathlib import Path
import importlib.util
import sys

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from participant_tfr_atlas import (  # noqa: E402
    EXPECTED_ASSIGNED_TRIALS,
    EXPECTED_BRIDGE_EXCLUDED,
    PRIMARY_CHANNELS,
    atomic_write_csv,
    atomic_write_json,
    compare_source_snapshots,
    load_and_validate_config,
    load_layout_positions,
    mean_memory_mapped_trials,
    participant_page_table,
    saturation_accounting,
    select_assigned_condition_trials,
    source_snapshot,
)
import participant_tfr_atlas as atlas  # noqa: E402


CONFIG = EEG_DIR / "participant_tfr_atlas_config_v1.yaml"
LAYOUT = REPO_ROOT / "_Scripts/_misc/sensor_latlong_chan_map.csv"


def test_layout_is_complete_physical_and_not_handedness_mirrored() -> None:
    """All 30 physical channels have the historical presentation coordinate once."""

    positions = load_layout_positions(LAYOUT, PRIMARY_CHANNELS)
    assert tuple(positions) == PRIMARY_CHANNELS
    assert len(positions) == 30
    # Historical radial coordinates place C3 left and C4 right. This atlas never
    # applies an analysis-hand-dependent coordinate transform.
    assert positions["C3"].x_fraction < positions["C4"].x_fraction
    assert positions["FP1"].y_fraction > positions["OZ"].y_fraction


def tiny_lineage() -> pd.DataFrame:
    """Return a minimal assigned-condition/bridge fixture with complete columns."""

    rows = []
    for participant, group, semantics, bridge in (
        (1, "physical", "overt_movement", False),
        (2, "imagery", "motor_imagery", False),
        (2, "imagery", "overt_movement", True),
    ):
        rows.append(
            {
                "behavioural_id": participant,
                "canonical_event_key": f"key-{len(rows)}",
                "recording_stem": f"demi_{participant:02d}_data",
                "tfr_row_index": len(rows),
                "group": group,
                "performed_condition": "imagery" if semantics == "motor_imagery" else "physical",
                "condition_semantics": semantics,
                "physical": semantics == "overt_movement",
                "scope_all_accepted": True,
                "scope_imagery_final_overt_bridge": bridge,
            }
        )
    return pd.DataFrame(rows)


def test_assigned_condition_selection_excludes_imagery_bridge_only(monkeypatch: pytest.MonkeyPatch) -> None:
    """The bridge is excluded only from atlas means, not asserted unavailable."""

    config = load_and_validate_config(CONFIG)[0]
    monkeypatch.setattr(atlas, "EXPECTED_BRIDGE_EXCLUDED", 1)
    selected = select_assigned_condition_trials(tiny_lineage(), config)
    assert selected["atlas_selected_for_mean"].tolist() == [True, True, False]
    assert selected["atlas_bridge_excluded_from_mean"].tolist() == [False, False, True]
    assert selected["scope_all_accepted"].all()


def test_page_assignment_is_one_page_per_participant(monkeypatch: pytest.MonkeyPatch) -> None:
    """Page assignment is stable within each PDF and one-to-one by participant."""

    config = load_and_validate_config(CONFIG)[0]
    monkeypatch.setattr(atlas, "EXPECTED_PARTICIPANTS", 2)
    monkeypatch.setattr(atlas, "EXPECTED_OVERT_PARTICIPANTS", 1)
    monkeypatch.setattr(atlas, "EXPECTED_IMAGERY_PARTICIPANTS", 1)
    monkeypatch.setattr(atlas, "EXPECTED_ASSIGNED_TRIALS", 2)
    monkeypatch.setattr(atlas, "EXPECTED_BRIDGE_EXCLUDED", 1)
    selected = select_assigned_condition_trials(tiny_lineage(), config)
    pages = participant_page_table(selected)
    assert pages[["behavioural_id", "pdf_group", "page_number"]].to_dict(orient="records") == [
        {"behavioural_id": 2, "pdf_group": "imagery", "page_number": 1},
        {"behavioural_id": 1, "pdf_group": "overt", "page_number": 1},
    ]


def test_memory_mapped_mean_uses_only_requested_rows(tmp_path: Path) -> None:
    """The helper opens an array read-only and averages selected trial cells exactly."""

    source = tmp_path / "db.npy"
    values = np.arange(3 * 2 * 2 * 2, dtype=np.float32).reshape(3, 2, 2, 2)
    np.save(source, values)
    observed = mean_memory_mapped_trials(source, [0, 2], (2, 2, 2))
    assert np.array_equal(observed, values[[0, 2]].mean(axis=0))
    with pytest.raises(RuntimeError, match="atlas_row_indices_invalid"):
        mean_memory_mapped_trials(source, [], (2, 2, 2))


def test_scale_and_saturation_accounting_are_rendered_cell_based() -> None:
    """Out-of-range values are counted before fixed -8/+8 image clipping."""

    onset = np.array([[[[-9.0]]]]).reshape(1, 1, 1)
    end = np.array([[[8.1]]])
    summary = saturation_accounting(onset, end, -8.0, 8.0)
    assert summary["below_minimum_count"] == 1
    assert summary["above_maximum_count"] == 1
    assert summary["saturated_cell_count"] == 2
    assert summary["rendered_source_cell_count"] == 2
    assert summary["saturated_percent"] == 100.0


def test_source_snapshot_detects_mutation(tmp_path: Path) -> None:
    """Source immutability validation fails after a source size/mtime change."""

    source = tmp_path / "source.txt"
    source.write_text("before", encoding="utf-8")
    before = source_snapshot([source], root=tmp_path)
    source.write_text("after-and-different", encoding="utf-8")
    after = source_snapshot([source], root=tmp_path)
    with pytest.raises(RuntimeError, match="atlas_source_size_or_mtime_changed"):
        compare_source_snapshots(before, after)


def test_atomic_writers_publish_reopenable_files(tmp_path: Path) -> None:
    """The small manifest/summary writers use atomic replacement and reopen checks."""

    json_path = tmp_path / "state.json"
    csv_path = tmp_path / "pages.csv"
    atomic_write_json(json_path, {"status": "complete", "trial_count": EXPECTED_ASSIGNED_TRIALS})
    atomic_write_csv(csv_path, pd.DataFrame({"participant": [4, 60], "page": [1, 2]}))
    assert json_path.read_text(encoding="utf-8").startswith("{")
    assert pd.read_csv(csv_path)["participant"].tolist() == [4, 60]
    assert not list(tmp_path.glob(".*.tmp-*"))


def test_completed_output_reuses_without_rewrite(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A matching complete atlas reopens and returns without touching artifact mtimes."""

    spec = importlib.util.spec_from_file_location(
        "atlas_driver", EEG_DIR / "20_construct_participant_sensor_tfr_atlas.py"
    )
    assert spec and spec.loader
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    monkeypatch.setattr(driver, "EXPECTED_PARTICIPANTS", 2)
    monkeypatch.setattr(driver, "EXPECTED_OVERT_PARTICIPANTS", 1)
    monkeypatch.setattr(driver, "EXPECTED_IMAGERY_PARTICIPANTS", 1)
    root = tmp_path / "atlas"
    root.mkdir()
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt

    pdfs = {}
    for group in ("overt", "imagery"):
        path = root / f"{group}.pdf"
        with PdfPages(path) as writer:
            figure = plt.figure(figsize=(1, 1))
            writer.savefig(figure)
            plt.close(figure)
        pdfs[group] = atlas.output_descriptor(path, root=root)
    pages = pd.DataFrame({"behavioural_id": [1, 2], "pdf_group": ["overt", "imagery"]})
    atomic_write_csv(root / "participant_page_manifest.csv", pages)
    atomic_write_json(root / "validation/atlas_validation.json", {"status": "pass"})
    outputs = {
        "pages": atlas.output_descriptor(root / "participant_page_manifest.csv", root=root),
        "validation": atlas.output_descriptor(root / "validation/atlas_validation.json", root=root),
        "overt": pdfs["overt"], "imagery": pdfs["imagery"],
    }
    authority = {"example": "matching"}
    atomic_write_json(root / driver.RUN_MANIFEST, {
        "status": "complete", "authority": authority, "outputs": outputs,
        "pdfs": {"overt": pdfs["overt"], "imagery": pdfs["imagery"]},
    })
    before = (root / "overt.pdf").stat().st_mtime_ns
    reused = driver.validate_existing(root, authority)
    assert reused is not None
    assert (root / "overt.pdf").stat().st_mtime_ns == before


def test_configuration_closes_scale_palette_and_bridge_contract() -> None:
    """The tracked config records the owner-selected display and selection rules."""

    config, _ = load_and_validate_config(CONFIG)
    assert config["display"]["palette"] == "plasma"
    assert (config["display"]["db_minimum"], config["display"]["db_maximum"]) == (-8.0, 8.0)
    assert config["accepted_surface"]["assigned_condition_trial_count"] == EXPECTED_ASSIGNED_TRIALS
    assert config["accepted_surface"]["imagery_final_overt_bridge_excluded_from_means"] == EXPECTED_BRIDGE_EXCLUDED
