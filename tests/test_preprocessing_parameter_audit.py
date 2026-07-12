"""Synthetic guardrail tests for the preprocessing-parameter evidence pass."""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import mne
import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_ANALYSIS_DIR = REPO_ROOT / "analysis" / "eeg_mne"
if str(EEG_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_ANALYSIS_DIR))

from preprocessing_parameter_audit import (  # noqa: E402
    analysis_filter_response,
    compare_line_noise_branches,
    load_audit_contract,
    make_uppercase_standard_1005,
    parse_python_assignments,
    parse_r_scalar_assignments,
    run_lof,
    validate_audit_output_path,
    welch_band_powers,
)


CONTRACT_PATH = EEG_ANALYSIS_DIR / "preprocessing_parameter_audit_contract_v1.yaml"


def test_tracked_contract_is_audit_only_and_preserves_60_hz() -> None:
    """The tracked contract cannot write signals/epochs or alias away 60 Hz."""

    contract = load_audit_contract(CONTRACT_PATH)
    assert contract["safety"]["write_signal_derivatives"] is False
    assert contract["safety"]["construct_epochs"] is False
    assert contract["prepared_continuous_input"]["resample_hz"] > 120
    assert contract["prepared_continuous_input"]["montage"] == "standard_1005"
    assert contract["autoreject"]["ransac_validation_in_this_audit"] is False


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("safety", "write_signal_derivatives", True),
        ("safety", "construct_epochs", True),
        ("safety", "mutate_source_edf", True),
        ("prepared_continuous_input", "montage", "standard_1020"),
        ("prepared_continuous_input", "apply_reference", True),
        ("prepared_continuous_input", "apply_notch", True),
        ("prepared_continuous_input", "apply_analysis_filter", True),
        ("prepared_continuous_input", "resample_hz", 100.0),
    ],
)
def test_contract_fails_closed_on_signal_or_input_scope_expansion(
    tmp_path: Path, section: str, key: str, value: object
) -> None:
    """Unsafe mutations of the bounded audit contract are rejected."""

    import yaml

    contract = copy.deepcopy(load_audit_contract(CONTRACT_PATH))
    contract[section][key] = value
    path = tmp_path / "unsafe.yaml"
    path.write_text(yaml.safe_dump(contract), encoding="utf-8")
    with pytest.raises(ValueError):
        load_audit_contract(path)


@pytest.mark.parametrize("name", ["subject_cleaned.fif", "subject.edf", "epochs.csv"])
def test_output_guard_rejects_signal_or_epoch_like_products(tmp_path: Path, name: str) -> None:
    """Audit helpers cannot be redirected into derivative or epoch products."""

    with pytest.raises(ValueError):
        validate_audit_output_path(tmp_path / name, tmp_path, (".csv", ".json", ".md"))


def test_output_guard_accepts_flat_local_evidence_path(tmp_path: Path) -> None:
    """A plain CSV evidence table directly under the audit root is allowed."""

    validate_audit_output_path(
        tmp_path / "global_bad_method_runs.csv", tmp_path, (".csv", ".json", ".md")
    )


def test_uppercase_standard_1005_matches_demi_labels() -> None:
    """The PyPREP-compatible montage preserves the approved MNE coordinates."""

    montage = make_uppercase_standard_1005()
    for channel in ("FP1", "FZ", "CZ", "CPZ", "M1", "M2", "OZ"):
        assert channel in montage.ch_names
    assert len(montage.ch_names) == len(set(montage.ch_names))


def test_literal_historical_setting_parsers_do_not_execute_code(tmp_path: Path) -> None:
    """Historical Python/R settings are recovered through bounded literal parsing."""

    python_path = tmp_path / "settings.py"
    python_path.write_text("enabled = True\nthreshold = 0.25\nignored = func()\n", encoding="utf-8")
    assert parse_python_assignments(python_path, ("enabled", "threshold")) == {
        "enabled": True,
        "threshold": 0.25,
    }

    r_path = tmp_path / "settings.R"
    r_path.write_text("drop_interpolated <- FALSE\nusing_csd <- TRUE\n", encoding="utf-8")
    assert parse_r_scalar_assignments(r_path, ("drop_interpolated", "using_csd")) == {
        "drop_interpolated": False,
        "using_csd": True,
    }


def synthetic_raw(duration_seconds: float = 30.0) -> mne.io.RawArray:
    """Return deterministic multichannel EEG with an explicit 60-Hz component."""

    sfreq = 250.0
    times = np.arange(int(sfreq * duration_seconds)) / sfreq
    channels = ["FP1", "FP2", "F3", "F4", "C3", "C4", "P3", "P4"]
    data = []
    for index, _ in enumerate(channels):
        phase = index * 0.1
        signal_10 = 10e-6 * np.sin(2 * np.pi * 10 * times + phase)
        signal_60 = (index + 1) * 1e-6 * np.sin(2 * np.pi * 60 * times)
        data.append(signal_10 + signal_60)
    raw = mne.io.RawArray(
        np.asarray(data),
        mne.create_info(channels, sfreq, ch_types="eeg"),
        verbose="ERROR",
    )
    raw.set_montage(make_uppercase_standard_1005(), on_missing="raise", verbose="ERROR")
    return raw


def test_lof_repeated_run_is_deterministic_on_same_input() -> None:
    """MNE LOF produces identical labels and scores on an unchanged input."""

    raw = synthetic_raw(8.0)
    first = run_lof(raw, n_neighbors=4, metric="euclidean", threshold=1.5)
    second = run_lof(raw, n_neighbors=4, metric="euclidean", threshold=1.5)
    assert first == second


def test_welch_band_power_detects_synthetic_line_component() -> None:
    """The audit PSD helper separates primary alpha and 60-Hz energy."""

    raw = synthetic_raw(20.0)
    powers = welch_band_powers(
        raw.get_data(),
        raw.info["sfreq"],
        bands={"alpha": (9.0, 12.0), "line_60": (59.0, 61.0)},
        n_fft=2048,
        n_overlap=1024,
    )
    assert np.all(powers["alpha"] > 0)
    assert powers["line_60"][-1] > powers["line_60"][0]


def test_line_noise_branches_preserve_alpha_and_remove_60_hz() -> None:
    """Synthetic comparison shows the intended branch distinction without writes."""

    raw = synthetic_raw(30.0)
    settings = copy.deepcopy(load_audit_contract(CONTRACT_PATH)["line_noise_filter_comparison"])
    settings["measurement_edge_exclusion_seconds"] = 2.0
    settings["welch_n_fft"] = 2048
    settings["welch_n_overlap"] = 1024
    rows, response = compare_line_noise_branches(raw, settings)
    assert len(rows) == len(raw.ch_names)
    assert np.median([row["notch_only_line_60_mean_psd"] for row in rows]) < np.median(
        [row["raw_line_60_mean_psd"] for row in rows]
    )
    assert np.median(
        [abs(row["analysis_b_vs_a_alpha_relative_change"]) for row in rows]
    ) < 0.01
    assert response["filter_length_seconds"] > 1.0
    assert response["response_db_60_hz"] < -40.0


def test_candidate_filter_response_covers_primary_bands() -> None:
    """The proposed FIR response is flat through theta, alpha, and beta."""

    response = analysis_filter_response(250.0, high_pass_hz=0.5, low_pass_hz=45.0)
    for key in ("response_db_4_hz", "response_db_8_hz", "response_db_12_hz", "response_db_30_hz"):
        assert abs(response[key]) < 0.1
    assert response["response_db_60_hz"] < -40.0
