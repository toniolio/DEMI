"""Focused safety contracts for Stage 19 EEG input isolation.

Tests use temporary synthetic Parquet files. They prove predictor-blind scale
projection, observed-outcome rejection, outcome-free design projection,
production output-path safety, accepted estimand classification, and optional
reopening of the ignored completed validation namespace. No scientific
estimate is calculated from accepted EEG.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = REPO_ROOT / "analysis/eeg_models"
sys.path.insert(0, str(MODEL_DIR))

from prepare_model_validation_inputs import (  # noqa: E402
    DESIGN_COLUMNS,
    PRIMARY_ESTIMANDS,
    SUPPLEMENTARY_ESTIMANDS,
    assert_synthetic_fit_input_safe,
    load_synthetic_design,
    observed_outcome_columns,
    pooled_scale_from_value_only,
    require_ignored_output_root,
)


def test_predictor_blind_scale_projection_and_definition(tmp_path: Path) -> None:
    """Scale calculation requests value_db only and uses one 1.4826 factor."""

    path = tmp_path / "scale.parquet"
    pd.DataFrame(
        {
            "value_db": [-2.0, -1.0, 0.0, 4.0],
            "condition_c": [-0.5, -0.5, 0.5, 0.5],
            "accuracy_rating_within": [99.0, 99.0, 99.0, 99.0],
        }
    ).to_parquet(path)
    result = pooled_scale_from_value_only(path)
    median = np.median([-2.0, -1.0, 0.0, 4.0])
    expected = 1.4826 * np.median(np.abs(np.array([-2.0, -1.0, 0.0, 4.0]) - median))
    assert result["columns_read"] == ["value_db"]
    assert result["m_y"] == median
    assert result["s_y"] == expected


def test_synthetic_design_projection_drops_observed_outcomes(tmp_path: Path) -> None:
    """Only the explicit design allowlist crosses the synthetic-data route."""

    rows = 7_307
    frame = pd.DataFrame(
        {
            "canonical_event_key": [f"key:{index}" for index in range(rows)],
            "participant_id": np.arange(rows) % 81,
            "condition_c": np.where(np.arange(rows) % 2, -0.5, 0.5),
            "accuracy_rating_within": np.zeros(rows),
            "accuracy_rating_between": np.zeros(rows),
            "familiarity_c": np.where(np.arange(rows) % 2, -0.5, 0.5),
            "value_db": np.arange(rows, dtype=float),
            "value_log_power": np.arange(rows, dtype=float),
        }
    )
    path = tmp_path / "h1.parquet"
    frame.to_parquet(path)
    design = load_synthetic_design(path, "H1")
    assert tuple(design.columns) == DESIGN_COLUMNS["H1"]
    assert not observed_outcome_columns(design.columns)
    assert_synthetic_fit_input_safe(design)


@pytest.mark.parametrize(
    "column", ["value_db", "value_log_power", "copied_value_db", "observed_eeg_theta"]
)
def test_synthetic_fit_rejects_observed_outcome_columns(column: str) -> None:
    """Known and copied outcome-like names fail closed before fitting."""

    with pytest.raises(ValueError, match="synthetic_fit_input_contains_observed_eeg"):
        assert_synthetic_fit_input_safe(pd.DataFrame({column: [0.0]}))


def test_estimand_classification_and_output_namespace_safety() -> None:
    """The accepted 9+3 registry and narrow ignored namespace are fixed."""

    assert PRIMARY_ESTIMANDS == ("T1", "T2", "T3", "A1", "A2", "A3", "B1", "B2", "B3")
    assert SUPPLEMENTARY_ESTIMANDS == ("B4", "B5", "B6")
    require_ignored_output_root(REPO_ROOT / "_Data/eeg/model_validation_v1")
    with pytest.raises(RuntimeError, match="unsafe_model_validation_output_root"):
        require_ignored_output_root(REPO_ROOT)


def test_completed_namespace_records_no_observed_fit_and_immutability() -> None:
    """If local production exists, reopen its core safety declarations."""

    root = REPO_ROOT / "_Data/eeg/model_validation_v1"
    validation_path = root / "validation/model_validation.json"
    if not validation_path.is_file():
        pytest.skip("local model_validation_v1 output absent")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    assert validation["status"] == "pass"
    assert validation["accepted_EEG_models_fitted"] is False
    assert validation["scientific_estimands_calculated"] is False
    assert validation["accepted_outcome_predictor_associations_opened"] is False
    assert validation["source_immutability"] == "pass"
    assert validation["atomic_directory_publication"] is True
    assert validation["validated_unchanged_reuse_supported"] is True

