"""Focused guardrail tests for config-driven raw EEG QC.

These tests exercise the fail-closed helper surface in
``analysis/eeg_mne/04_preprocess_continuous_raw.py`` without opening private
raw EDF files. Synthetic configs, small DataFrames, and temporary paths keep
the tests public and independent of local data.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "analysis" / "eeg_mne" / "04_preprocess_continuous_raw.py"


def load_script_module():
    """Import script 04 by path because its filename is not a valid module name."""

    spec = importlib.util.spec_from_file_location("raw_qc_script", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


raw_qc = load_script_module()


def minimal_valid_config(output_root: str = "_Data/eeg/mne_preprocessing/raw_qc") -> dict[str, Any]:
    """Return the smallest config shape needed for fail-closed validation."""

    return {
        "outputs": {
            "output_root": output_root,
            "local_only": True,
            "write_cleaned_signal_derivatives": False,
            "write_cleaned_fif": False,
            "write_preprocessed_edf": False,
            "write_epochs": False,
        },
        "channels": {
            "type_mapping": {
                "eeg": list(raw_qc.EXPECTED_EEG_CHANNELS),
                "eog": list(raw_qc.EXPECTED_EOG_CHANNELS),
                "emg": list(raw_qc.EXPECTED_EMG_CHANNELS),
                "stim": ["Trigger"],
            },
            "future_mne_type_overrides": {
                "HEO": "eog",
                "VEO": "eog",
                "EMG-L": "emg",
                "EMG-A": "emg",
                "Trigger": "stim",
            },
            "trigger_policy": {
                "use_for_event_repair": False,
                "use_for_epoching": False,
            },
        },
        "montage_validation": {
            "apply_active_montage_to_signal": False,
        },
        "raw_qc_outputs": {
            "enabled": True,
            "signal_mutation_allowed": False,
            "event_context": {
                "make_event_policy_decisions": False,
            },
        },
        "disabled_signal_mutating_sections": {
            name: {"enabled": False}
            for name in (
                "filtering",
                "notch_filtering",
                "reference",
                "bad_channel_detection",
                "interpolation",
                "ica",
                "csd",
                "cleaned_fif_writing",
                "event_repair",
                "epoching",
            )
        },
    }


def assert_config_valid(config: dict[str, Any], tmp_path: Path) -> None:
    """Run validation against a temporary repository root."""

    raw_qc.validate_raw_qc_config(config, tmp_path / "raw_qc_config.yaml", tmp_path)


def set_nested(config: dict[str, Any], path: tuple[str, ...], value: Any) -> None:
    """Set a nested config value for mutation tests."""

    current: dict[str, Any] = config
    for key in path[:-1]:
        current = current[key]
    current[path[-1]] = value


def file_row(
    participant_id: int = 1,
    source_filename: str = "demi_01 Data.edf",
    file_role: str = "single",
    annotation_count: int | None = 700,
) -> dict[str, Any]:
    """Return a small file-row dictionary used by table helpers."""

    return {
        "participant_id": participant_id,
        "participant_id_padded": f"{participant_id:03d}",
        "source_filename": source_filename,
        "file_path": f"_Data/eeg/raw/{source_filename}",
        "manifest_present": True,
        "discovered_in_raw_dir": True,
        "file_exists": True,
        "file_role": file_role,
        "split_part": None,
        "manifest_annotation_count": annotation_count,
        "annotation_count_from_counts_csv": annotation_count,
    }


def header_for_channels(
    channel_names: list[str],
    channel_types_after: list[str] | None = None,
    annotation_count: int | None = 700,
) -> dict[str, Any]:
    """Return a synthetic successful MNE-header dictionary."""

    if channel_types_after is None:
        channel_types_after = [raw_qc.expected_channel_type(raw_qc.normalize_channel_name(name)) for name in channel_names]
    return {
        "read_status": "ok",
        "mne_read_error_type": "",
        "mne_read_error_message": "",
        "sampling_frequency_hz": 1000.0,
        "duration_seconds": 10.0,
        "channel_names": channel_names,
        "mne_channel_types_before_overrides": ["eeg"] * len(channel_names),
        "channel_types_after_overrides": channel_types_after,
        "applied_channel_type_overrides_json": "{}",
        "channel_type_override_status": "applied",
        "channel_type_override_error_type": "",
        "channel_type_override_error_message": "",
        "mne_annotation_count": annotation_count,
        "mne_annotation_description_counts_json": "{}",
    }


def test_minimal_valid_config_passes(tmp_path: Path) -> None:
    """A minimal config preserving the approved raw-QC boundary should pass."""

    assert_config_valid(minimal_valid_config(), tmp_path)


def test_missing_and_malformed_config_fail(tmp_path: Path) -> None:
    """Config loading should fail before validation when YAML is absent or bad."""

    with pytest.raises(RuntimeError, match="missing"):
        raw_qc.load_raw_qc_config(tmp_path / "missing.yaml")

    malformed = tmp_path / "malformed.yaml"
    malformed.write_text("outputs: [unterminated\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="failed to parse"):
        raw_qc.load_raw_qc_config(malformed)

    scalar = tmp_path / "scalar.yaml"
    scalar.write_text("not-a-mapping\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="mapping"):
        raw_qc.load_raw_qc_config(scalar)


@pytest.mark.parametrize(
    "section",
    [
        "filtering",
        "notch_filtering",
        "reference",
        "bad_channel_detection",
        "interpolation",
        "ica",
        "csd",
        "cleaned_fif_writing",
        "event_repair",
        "epoching",
    ],
)
def test_each_disabled_stage_fails_when_enabled(section: str, tmp_path: Path) -> None:
    """Every future signal/event stage should fail closed if enabled."""

    config = minimal_valid_config()
    config["disabled_signal_mutating_sections"][section]["enabled"] = True
    with pytest.raises(RuntimeError, match=section):
        assert_config_valid(config, tmp_path)


@pytest.mark.parametrize(
    "path",
    [
        ("outputs", "write_cleaned_signal_derivatives"),
        ("outputs", "write_cleaned_fif"),
        ("outputs", "write_preprocessed_edf"),
        ("outputs", "write_epochs"),
        ("raw_qc_outputs", "signal_mutation_allowed"),
        ("raw_qc_outputs", "event_context", "make_event_policy_decisions"),
        ("channels", "trigger_policy", "use_for_event_repair"),
        ("channels", "trigger_policy", "use_for_epoching"),
        ("montage_validation", "apply_active_montage_to_signal"),
    ],
)
def test_signal_derivative_and_event_policy_switches_fail_when_true(path: tuple[str, ...], tmp_path: Path) -> None:
    """Cleaned-output, signal-mutation, and event-policy switches stay false."""

    config = minimal_valid_config()
    set_nested(config, path, True)
    with pytest.raises(RuntimeError, match="must stay false"):
        assert_config_valid(config, tmp_path)


def test_output_root_must_stay_under_data(tmp_path: Path) -> None:
    """Only output roots below repository-local _Data/ are accepted."""

    assert_config_valid(minimal_valid_config("_Data/eeg/raw_qc_test"), tmp_path)

    outside_data = minimal_valid_config("outputs/raw_qc")
    with pytest.raises(RuntimeError, match="_Data"):
        assert_config_valid(outside_data, tmp_path)

    outside_repo = minimal_valid_config((tmp_path.parent / "raw_qc").as_posix())
    with pytest.raises(RuntimeError, match="inside the repository"):
        assert_config_valid(outside_repo, tmp_path)


def test_changed_channel_mapping_or_missing_override_fails(tmp_path: Path) -> None:
    """The approved DEMI channel map and required overrides are fixed."""

    changed_mapping = minimal_valid_config()
    changed_mapping["channels"]["type_mapping"]["eog"].remove("VEO")
    with pytest.raises(RuntimeError, match="channel mapping differs"):
        assert_config_valid(changed_mapping, tmp_path)

    missing_override = minimal_valid_config()
    del missing_override["channels"]["future_mne_type_overrides"]["EMG-A"]
    with pytest.raises(RuntimeError, match="EMG-A=emg"):
        assert_config_valid(missing_override, tmp_path)


def test_configured_mapping_normalizes_expected_channel_names() -> None:
    """Configured labels should normalize to the approved channel-type map."""

    mapping = raw_qc.configured_channel_type_by_normalized_name(minimal_valid_config())
    assert mapping["HEO"] == "eog"
    assert mapping["VEO"] == "eog"
    assert mapping["EMG-L"] == "emg"
    assert mapping["EMG-A"] == "emg"
    assert mapping["TRIGGER"] == "stim"


def test_channel_rows_surface_unknown_names() -> None:
    """Unknown channel labels should produce explicit unmapped rows."""

    header = header_for_channels(
        ["HEO", "VEO", "EMG-L", "EMG-A", "Trigger", "Mystery"],
        ["eog", "eog", "emg", "emg", "stim", "eeg"],
    )
    rows = raw_qc.build_channel_type_rows(file_row(), header)
    by_name = {row["normalized_channel_name"]: row for row in rows}

    assert by_name["HEO"]["expected_raw_qc_type"] == "eog"
    assert by_name["VEO"]["expected_raw_qc_type"] == "eog"
    assert by_name["EMG-L"]["expected_raw_qc_type"] == "emg"
    assert by_name["EMG-A"]["expected_raw_qc_type"] == "emg"
    assert by_name["TRIGGER"]["expected_raw_qc_type"] == "stim"
    assert by_name["MYSTERY"]["channel_type_assignment_status"] == "unmapped_channel_label"
    assert by_name["MYSTERY"]["channel_type_issue_codes"] == "unmapped_channel_label"


def test_besa_name_normalization_matches_case_variants_and_surfaces_missing_coordinates() -> None:
    """BESA matching should normalize case only and expose missing coordinates."""

    besa_channels = {
        channel: channel
        for channel in raw_qc.EXPECTED_EEG_CHANNELS
        if channel != "PZ"
    }
    besa_channels.update(
        {
            "FP1": "Fp1",
            "FZ": "Fz",
            "FCZ": "FCz",
            "CPZ": "CPz",
            "CZ": "Cz",
            "OZ": "Oz",
        }
    )
    coordinates = pd.DataFrame(
        [
            {"chan": label, "x": index, "y": index + 1, "z": index + 2}
            for index, label in enumerate(besa_channels.values())
        ]
    )
    coordinates["normalized_channel_name"] = coordinates["chan"].map(raw_qc.normalize_channel_name)
    besa_by_name = {
        row["normalized_channel_name"]: row
        for _, row in coordinates.iterrows()
    }

    header = header_for_channels(list(raw_qc.EXPECTED_EEG_CHANNELS))
    rows = raw_qc.build_montage_coordinate_rows(file_row(), header, besa_by_name, None)
    by_name = {row["normalized_channel_name"]: row for row in rows}

    for name in ["FP1", "FZ", "FCZ", "CPZ", "CZ", "OZ"]:
        assert by_name[name]["besa_coordinate_match"] is True
        assert by_name[name]["montage_issue_codes"] == ""

    assert by_name["PZ"]["besa_coordinate_match"] is False
    assert "expected_eeg_channel_without_besa_coordinate" in by_name["PZ"]["montage_issue_codes"]


@pytest.mark.parametrize(
    ("source_filename", "expected_stem"),
    [
        ("demi_01 Data.edf", "demi_01_Data"),
        ("demi_54_1 Data.edf", "demi_54_1_Data"),
        ("demi_05 concatenated Data.edf", "demi_05_concatenated_Data"),
    ],
)
def test_safe_figure_stem_is_stable_for_common_file_roles(source_filename: str, expected_stem: str) -> None:
    """Figure stems should be deterministic for normal, split, and concatenated files."""

    assert raw_qc.safe_figure_stem({"source_filename": source_filename}) == expected_stem


@pytest.mark.parametrize("participant_id", [54, 56, 65])
def test_split_file_ids_get_split_context_code(participant_id: int) -> None:
    """Known split-file IDs should carry factual split event-context codes."""

    codes = raw_qc.pending_event_context_codes(participant_id, "single", 700)
    assert "event_policy_pending_split_file_continuity" in codes


def test_id5_gets_concatenated_file_start_context_code() -> None:
    """ID 5 should carry the concatenated/file-start context code."""

    codes = raw_qc.pending_event_context_codes(5, "single", 700)
    assert "event_policy_pending_concatenated_file_start" in codes


@pytest.mark.parametrize("participant_id", [11, 14, 89, 94, 100])
def test_zero_no_offset_ids_get_context_code(participant_id: int) -> None:
    """Known zero/no-offset IDs should carry factual pending context codes."""

    codes = raw_qc.pending_event_context_codes(participant_id, "single", 700)
    assert "event_policy_pending_zero_or_no_offset_context" in codes


def test_zero_and_low_annotation_flags_remain_factual_context() -> None:
    """Zero/low annotation counts should be facts in outputs, not event decisions."""

    zero_file = file_row(participant_id=11, source_filename="demi_11 Data.edf", annotation_count=0)
    low_file = file_row(participant_id=89, source_filename="demi_89 Data.edf", annotation_count=50)
    montage_rows = []

    zero_summary = raw_qc.build_file_summary_row(zero_file, header_for_channels([], [], 0), montage_rows)
    low_summary = raw_qc.build_file_summary_row(low_file, header_for_channels([], [], 50), montage_rows)
    summary = pd.DataFrame([zero_summary, low_summary])
    event_context = raw_qc.build_event_context_flags(summary, {})

    zero_row = event_context[event_context["source_filename"].eq("demi_11 Data.edf")].iloc[0]
    low_row = event_context[event_context["source_filename"].eq("demi_89 Data.edf")].iloc[0]

    assert bool(zero_row["zero_annotation_fact"])
    assert "zero_annotation_fact_only" in zero_row["event_context_pending_codes"]
    assert "event_policy_pending_zero_or_no_offset_context" in zero_row["event_context_pending_codes"]
    assert bool(low_row["low_annotation_fact"])
    assert "low_annotation_fact_only" in low_row["event_context_pending_codes"]
    assert "event_policy_pending_zero_or_no_offset_context" in low_row["event_context_pending_codes"]


def test_valid_minimal_config_round_trips_through_yaml_loader(tmp_path: Path) -> None:
    """Synthetic valid configs should parse from YAML and validate."""

    config_path = tmp_path / "raw_qc_config.yaml"
    config_path.write_text(yaml.safe_dump(minimal_valid_config()), encoding="utf-8")
    loaded = raw_qc.load_raw_qc_config(config_path)
    raw_qc.validate_raw_qc_config(loaded, config_path, tmp_path)
