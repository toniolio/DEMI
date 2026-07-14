"""Synthetic atomic-storage, ledger, resume, and derivative validation tests."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_DIR))

from continuous_preprocessing.contracts import (  # noqa: E402
    ALL_RETAINED_CHANNELS,
    load_config,
)
from continuous_preprocessing.cohort import select_production_surface  # noqa: E402
from continuous_preprocessing.runner import (  # noqa: E402
    aggregate_continuous_run,
    compare_source_inventories,
    execute_members,
    reopen_current_surface_derivatives,
    snapshot_source_inventory,
)
from continuous_preprocessing.source import apply_channel_types, apply_montage  # noqa: E402
from continuous_preprocessing.stages import (  # noqa: E402
    apply_fir_filter,
    apply_reference,
    calculate_reference_sources,
    interpolate_global_bads,
)
from continuous_preprocessing.storage import (  # noqa: E402
    ResultStore,
    StageLedger,
    artifact_descriptor,
    atomic_write_json,
    deterministic_manifest_id,
    recording_provenance,
    validate_derivative_path,
    validate_output_root,
    validate_saved_raw,
    write_raw_fif_atomic,
)


CONFIG_PATH = EEG_DIR / "continuous_preprocessing_config_v1.yaml"


def prepared_raw() -> tuple[mne.io.RawArray, dict]:
    """Return a small analysis-filtered/referenced/interpolated synthetic Raw."""

    cfg = load_config(CONFIG_PATH)
    rng = np.random.default_rng(4)
    data = rng.normal(0, 1e-6, (len(ALL_RETAINED_CHANNELS), 2500))
    data[-1] = 0.0
    raw = mne.io.RawArray(
        data,
        mne.create_info(list(ALL_RETAINED_CHANNELS), 250.0, "eeg"),
        verbose="ERROR",
    )
    apply_channel_types(raw)
    apply_montage(raw)
    apply_fir_filter(
        raw,
        high_pass_hz=0.5,
        low_pass_hz=45.0,
        settings=cfg["analysis_filter"],
    )
    reference = calculate_reference_sources(["F7"], cfg)
    apply_reference(raw, reference)
    interpolate_global_bads(raw, ["F7"], cfg)
    return raw, cfg


def roots(tmp_path: Path) -> tuple[Path, Path]:
    """Create synthetic authorized and immutable-root stand-ins."""

    output = tmp_path / "continuous_validation_v1"
    raw = tmp_path / "raw"
    output.mkdir()
    raw.mkdir()
    return output, raw


def publish_fake_result(
    store: ResultStore, output_root: Path, raw_root: Path, status: str = "complete"
) -> Path:
    """Publish one minimal hash-valid terminal result for cache tests."""

    decision = store.prepare()
    assert decision["action"] == "process"
    evidence = store.work_dir / "evidence.json"
    atomic_write_json(evidence, {"ok": True}, output_root, raw_root)
    descriptor = artifact_descriptor(evidence, store.work_dir)
    manifest = {
        "status": status,
        "provenance": {"fingerprint": store.expected_provenance["fingerprint"]},
        "artifacts": [descriptor],
        "manifest_id": "synthetic",
    }
    atomic_write_json(store.work_dir / "manifest.json", manifest, output_root, raw_root)
    return store.publish()


def test_stage_ledger_transitions_and_not_reached_states(tmp_path: Path) -> None:
    """The ordered ledger flushes complete, stopped, and not-reached states."""

    output, raw_root = roots(tmp_path)
    ledger = StageLedger(output / "ledger.json", output, raw_root)
    ledger.start("source_validation_and_read")
    ledger.complete("source_validation_and_read", {"ok": True})
    ledger.start("channel_typing")
    ledger.stop("channel_typing", "synthetic_stop", "test")
    ledger.finalize("stopped")
    payload = json.loads((output / "ledger.json").read_text(encoding="utf-8"))
    assert payload["overall_status"] == "stopped"
    assert payload["stages"][0]["status"] == "complete"
    assert payload["stages"][1]["status"] == "stopped"
    assert all(row["status"] == "not_reached" for row in payload["stages"][2:])


def test_output_path_safety_prohibits_raw_edf_epoch_and_escape(tmp_path: Path) -> None:
    """Derivative writes cannot target raw, EDF, epoch-like, or escaped paths."""

    output, raw_root = roots(tmp_path)
    validate_derivative_path(output / "recordings" / "x" / "post_ica_raw.fif", output, raw_root)
    with pytest.raises(ValueError, match="raw EDF"):
        validate_derivative_path(raw_root / "derived_raw.fif", output, raw_root)
    with pytest.raises(ValueError, match="EDF"):
        validate_derivative_path(output / "derived.edf", output, raw_root)
    with pytest.raises(ValueError, match="Prohibited"):
        validate_derivative_path(output / "subject_epochs.fif", output, raw_root)
    with pytest.raises(ValueError, match="escapes"):
        validate_derivative_path(tmp_path / "elsewhere.json", output, raw_root)


def test_atomic_derivative_write_reopen_and_full_validation(tmp_path: Path) -> None:
    """A synthetic saved continuous derivative reopens with all required metadata."""

    output, raw_root = roots(tmp_path)
    raw, cfg = prepared_raw()
    path = output / "post_ica_raw.fif"
    written = write_raw_fif_atomic(
        raw,
        path,
        fif_format=cfg["derivatives"]["fif_format"],
        output_root=output,
        raw_root=raw_root,
    )
    validation = validate_saved_raw(
        path,
        expected_sample_count=raw.n_times,
        expected_sfreq=raw.info["sfreq"],
        expected_highpass=0.5,
        expected_lowpass=45.0,
        expected_sha256=written["sha256"],
        chunk_seconds=2.0,
    )
    assert validation["valid"]
    assert validation["m1_m2_retained"]
    assert validation["montage_available"]
    assert validation["custom_reference_applied"]
    assert validation["operational_bad_list"] == []
    assert validation["nonfinite_sample_count"] == 0
    assert not list(output.glob(".*tmp*"))


def test_unchanged_terminal_result_is_skipped_and_force_is_explicit(tmp_path: Path) -> None:
    """Hash-valid terminal results skip; explicit force preserves then recomputes."""

    output, raw_root = roots(tmp_path)
    provenance = {"fingerprint": "same"}
    store = ResultStore(output, raw_root, "demi_01_data", provenance)
    publish_fake_result(store, output, raw_root)
    unchanged = ResultStore(output, raw_root, "demi_01_data", provenance).prepare()
    assert unchanged["action"] == "skip"
    forced = ResultStore(output, raw_root, "demi_01_data", provenance, force=True).prepare()
    assert forced["action"] == "process"
    assert forced["reason"] == "force_recompute"
    assert any("__force__" in path.name for path in (output / "history").iterdir())


def test_source_config_code_hash_changes_invalidate_terminal_result(tmp_path: Path) -> None:
    """A changed provenance fingerprint archives a result as stale before recompute."""

    output, raw_root = roots(tmp_path)
    first = {"fingerprint": "first"}
    publish_fake_result(ResultStore(output, raw_root, "demi_01_data", first), output, raw_root)
    stale = ResultStore(output, raw_root, "demi_01_data", {"fingerprint": "changed"}).prepare()
    assert stale["action"] == "process"
    assert stale["previous_assessment"]["cache_status"] == "stale"
    assert any("__stale__" in path.name for path in (output / "history").iterdir())

    base_run = {
        "pipeline_version": "continuous_preprocessing_v1",
        "config_sha256": "config-a",
        "code": {"sha256": "code-a"},
        "environment_sha256": "env-a",
    }
    base = recording_provenance(base_run, "source-a")
    assert recording_provenance(base_run, "source-b")["fingerprint"] != base["fingerprint"]
    changed_config = copy.deepcopy(base_run)
    changed_config["config_sha256"] = "config-b"
    assert recording_provenance(changed_config, "source-a")["fingerprint"] != base["fingerprint"]
    changed_code = copy.deepcopy(base_run)
    changed_code["code"]["sha256"] = "code-b"
    assert recording_provenance(changed_code, "source-a")["fingerprint"] != base["fingerprint"]


def test_interrupted_incomplete_directory_is_never_valid_and_is_preserved(tmp_path: Path) -> None:
    """An interrupted work directory is moved to history, never treated as terminal."""

    output, raw_root = roots(tmp_path)
    provenance = {"fingerprint": "same"}
    first = ResultStore(output, raw_root, "demi_01_data", provenance)
    assert first.prepare()["action"] == "process"
    atomic_write_json(first.work_dir / "partial.json", {"stage": 3}, output, raw_root)
    resumed = ResultStore(output, raw_root, "demi_01_data", provenance).prepare()
    assert resumed["action"] == "process"
    assert any("__incomplete__" in path.name for path in (output / "history").iterdir())
    assert resumed["result_dir"].name.endswith(".incomplete")


def test_continuation_after_one_file_failure() -> None:
    """Sequential execution records one failure and still processes the next file."""

    seen: list[str] = []
    members = [
        {"source_filename": "one.edf", "source_path": "/missing/one.edf"},
        {"source_filename": "two.edf", "source_path": "/missing/two.edf"},
    ]

    def processor(member):
        seen.append(member["source_filename"])
        if member["source_filename"] == "one.edf":
            raise RuntimeError("synthetic failure")
        return {"source_filename": "two.edf", "status": "complete"}

    results = execute_members(members, processor)
    assert seen == ["one.edf", "two.edf"]
    assert [row["status"] for row in results] == ["failed", "complete"]


def test_deterministic_manifest_id_ignores_runtime_and_output_container_noise() -> None:
    """Stable scientific/provenance content produces one deterministic ID."""

    manifest = {
        "schema_version": 1,
        "pipeline_version": "continuous_preprocessing_v1",
        "source": {"source_sha256": "source"},
        "provenance": {"fingerprint": "fingerprint"},
        "status": "complete",
        "channel_contract": {"targets": 32},
        "detector": {"bads": ["F7"]},
        "reference": {"count": 29},
        "interpolation": {"count": 1},
        "ica": {"exclude": [0]},
        "created_at": "first",
        "runtime": {"total_elapsed_seconds": 1.0},
        "artifacts": [{"sha256": "first-output"}],
    }
    changed = copy.deepcopy(manifest)
    changed["created_at"] = "second"
    changed["runtime"]["total_elapsed_seconds"] = 999.0
    changed["artifacts"] = [{"sha256": "second-output"}]
    assert deterministic_manifest_id(manifest) == deterministic_manifest_id(changed)


def test_production_surface_uses_all_readable_inventory_rows_in_stable_order(
    tmp_path: Path,
) -> None:
    """Production selection keeps split parts separate and ignores event eligibility."""

    raw_root = tmp_path / "_Data" / "eeg" / "raw"
    raw_root.mkdir(parents=True)
    filenames = [
        "demi_54_2 Data.edf",
        "demi_02 Data.edf",
        "demi_54_1 Data.edf",
    ]
    for filename in filenames:
        (raw_root / filename).write_bytes(b"synthetic-edf")
    manifest_path = tmp_path / "raw_eeg_file_manifest.csv"
    pd.DataFrame(
        [
            {
                "source_filename": filename,
                "file_path": f"_Data/eeg/raw/{filename}",
                "file_role": "split_part" if "54_" in filename else "single",
                "split_part": (
                    int(filename.split("_")[2].split()[0]) if "54_" in filename else None
                ),
                "read_status": "ok",
                "participant_id": 54 if "54_" in filename else 2,
            }
            for filename in filenames
        ]
    ).to_csv(manifest_path, index=False)

    surface = select_production_surface(
        raw_root=raw_root,
        raw_manifest_path=manifest_path,
        expected_readable_count=3,
    )
    assert [row["source_filename"] for row in surface["members"]] == [
        "demi_02 Data.edf",
        "demi_54_1 Data.edf",
        "demi_54_2 Data.edf",
    ]
    assert [row["split_part"] for row in surface["members"]] == [None, 1, 2]
    assert surface["process_split_parts_independently"] is True
    assert surface["selection_uses_event_or_epoch_eligibility"] is False
    assert all(
        not row["successful_continuous_preprocessing_implies_epoch_eligibility"]
        for row in surface["members"]
    )


def test_output_root_guard_accepts_only_separate_versioned_namespaces(
    tmp_path: Path,
) -> None:
    """The validation and production roots are exact; arbitrary roots fail closed."""

    for relative in (
        "_Data/eeg/mne_preprocessing/continuous_validation_v1",
        "_Data/eeg/mne_preprocessing/continuous_v1",
    ):
        validate_output_root(tmp_path, tmp_path / relative, relative)
    with pytest.raises(ValueError, match="Only the authorized"):
        validate_output_root(
            tmp_path,
            tmp_path / "_Data/eeg/mne_preprocessing/other",
            "_Data/eeg/mne_preprocessing/other",
        )


def test_complete_source_inventory_snapshot_detects_content_and_mtime_change(
    tmp_path: Path,
) -> None:
    """Run-level raw immutability compares hash, size, path, and timestamp."""

    source = tmp_path / "demi_01 Data.edf"
    source.write_bytes(b"before")
    members = [{"order": 1, "source_filename": source.name, "source_path": source.as_posix()}]
    before = snapshot_source_inventory(members)
    unchanged = snapshot_source_inventory(members)
    assert compare_source_inventories(before, unchanged)["unchanged"] is True
    source.write_bytes(b"after-with-different-size")
    after = snapshot_source_inventory(members)
    comparison = compare_source_inventories(before, after)
    assert comparison["unchanged"] is False
    assert comparison["changed_files"][0]["source_filename"] == source.name


def test_execution_namespace_changes_recording_cache_fingerprint() -> None:
    """Validation and production execution contexts cannot share cache identity."""

    base_run = {
        "pipeline_version": "continuous_preprocessing_v1",
        "config_sha256": "config",
        "code": {"sha256": "code"},
        "environment_sha256": "environment",
        "execution_sha256": "validation",
    }
    validation = recording_provenance(base_run, "source")
    production_run = {**base_run, "execution_sha256": "production"}
    production = recording_provenance(production_run, "source")
    assert validation["fingerprint"] != production["fingerprint"]


def test_current_surface_reopen_verifier_rescans_saved_raw(tmp_path: Path) -> None:
    """Post-run verification reopens current signal FIFs independently."""

    output, raw_root = roots(tmp_path)
    result_dir = output / "recordings" / "demi_01_data"
    result_dir.mkdir(parents=True)
    raw, cfg = prepared_raw()
    written = write_raw_fif_atomic(
        raw,
        result_dir / "post_ica_raw.fif",
        fif_format="single",
        output_root=output,
        raw_root=raw_root,
    )
    manifest = {
        "status": "complete",
        "source": {
            "sample_count": int(raw.n_times),
            "sampling_frequency_hz": float(raw.info["sfreq"]),
        },
        "ica": {"application": {"final_exclusions": []}},
        "artifacts": [
            {
                "kind": "post_ica_continuous",
                "relative_path": "post_ica_raw.fif",
                "sha256": written["sha256"],
                "size_bytes": written["size_bytes"],
            }
        ],
    }
    atomic_write_json(result_dir / "manifest.json", manifest, output, raw_root)
    surface = {
        "surface_kind": "synthetic",
        "surface_size": 1,
        "members": [{"source_filename": "demi_01 Data.edf"}],
    }
    verification = reopen_current_surface_derivatives(
        output_root=output,
        raw_root=raw_root,
        surface=surface,
        config=cfg,
    )
    assert verification["continuous_raw_artifact_count"] == 1
    assert verification["all_current_fif_and_ica_artifacts_reopened_and_valid"]


def test_bounded_production_invocation_is_not_labeled_full_dataset_run(
    tmp_path: Path,
) -> None:
    """A finalized smoke remains explicitly partial until all members are attempted."""

    surface = {
        "surface_kind": "authorized_all_inventoried_readable_edfs_v1",
        "surface_size": 95,
        "members": [
            {"source_filename": f"demi_{index:03d} Data.edf"}
            for index in range(1, 96)
        ],
    }
    aggregate = aggregate_continuous_run(
        repo_root=tmp_path,
        output_root=tmp_path / "continuous_v1",
        surface=surface,
        invocation_results=[],
        run_provenance_data={},
        not_attempted=[row["source_filename"] for row in surface["members"]][1:],
        run_started_at="synthetic",
        run_elapsed_seconds=0.0,
        run_state="finalized",
        preflight={},
        source_inventory_before=None,
        source_inventory_after=None,
    )
    assert aggregate["validation"]["full_dataset_run_performed"] is False
    assert aggregate["validation"]["all_surface_members_have_terminal_manifests"] is False
    assert aggregate["current_surface_terminal_counts"]["missing_or_incomplete"] == 95
