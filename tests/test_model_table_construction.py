"""Focused contracts for stage-18 DEMI EEG model-ready tables.

The pure contract tests require only tracked code/configuration. Integration
tests reopen the ignored production namespace when it is locally available;
they skip cleanly in public clones without private data.
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
EEG_DIR = REPO_ROOT / "analysis" / "eeg_mne"
sys.path.insert(0, str(EEG_DIR))

from model_table_construction import (  # noqa: E402
    ACCURACY_BETWEEN_CENTER,
    ALL_MODEL_ROLES,
    EXPECTED_GENERAL_ROWS,
    EXPECTED_VIEW_ROWS,
    FAMILYWISE_SUPPORT_PERCENT,
    LOG_ERROR_BETWEEN_CENTER,
    PRIMARY_ESTIMAND_IDS,
    SUPPLEMENTARY_ESTIMAND_IDS,
    archive_existing_output,
    estimand_registry,
    factor_coding_registry,
    load_and_validate_model_table_config,
    make_staging_directory,
    model_role_registry,
    publish_staging_directory,
    require_ignored_output_root,
)


CONFIG_PATH = EEG_DIR / "model_table_config_v1.yaml"
PRODUCTION_ROOT = REPO_ROOT / "_Data/eeg/model_tables_v1"
PRODUCTION_MANIFEST = PRODUCTION_ROOT / "model_table_run_manifest.json"


def load_config() -> dict:
    """Return the validated tracked stage-18 configuration."""

    config, _ = load_and_validate_model_table_config(CONFIG_PATH)
    return config


def test_exact_five_role_mapping_and_factor_codes() -> None:
    """The tracked contract selects only the five accepted ROI rows."""

    config = load_config()
    roles = model_role_registry(config)
    assert tuple(roles["model_role"]) == ALL_MODEL_ROLES
    assert len(roles) == 5
    assert roles[["roi_name", "feature_definition", "band"]].to_dict("records") == [
        {
            "roi_name": "frontal_theta",
            "feature_definition": "response_onset_0_1",
            "band": "theta",
        },
        {
            "roi_name": "frontal_theta",
            "feature_definition": "response_onset_0_0p5",
            "band": "theta",
        },
        {
            "roi_name": "posterior_alpha",
            "feature_definition": "response_end_0_1",
            "band": "alpha",
        },
        {
            "roi_name": "contralateral_motor_beta",
            "feature_definition": "response_end_pre_0p5_0",
            "band": "beta",
        },
        {
            "roi_name": "contralateral_motor_beta",
            "feature_definition": "response_end_pmbr_0p5_1p5",
            "band": "beta",
        },
    ]
    assert roles.loc[roles["model_role"].eq("beta_pre"), "phase_c"].item() == -0.5
    assert roles.loc[roles["model_role"].eq("beta_pmbr"), "phase_c"].item() == 0.5
    assert roles.loc[~roles["model_role"].str.startswith("beta"), "phase_c"].isna().all()

    factors = factor_coding_registry(config)
    observed = {
        (row.factor, row.source_label): row.code
        for row in factors.itertuples(index=False)
    }
    assert observed == {
        ("condition", "physical"): -0.5,
        ("condition", "imagery"): 0.5,
        ("familiarity", "random"): -0.5,
        ("familiarity", "repeated"): 0.5,
        ("phase", "pre_end"): -0.5,
        ("phase", "pmbr"): 0.5,
    }


def test_predictor_centers_scope_and_no_unauthorized_transformation() -> None:
    """The rating/error derivation constants and safety boundary are exact."""

    config = load_config()
    accuracy = config["predictors"]["accuracy_rating"]
    error = config["predictors"]["objective_error"]
    assert float(accuracy["between_center"]) == ACCURACY_BETWEEN_CENTER
    assert accuracy["derivation_scope"] == "primary_blocks_1_5"
    assert accuracy["units"] == "rating_points"
    assert not any(
        accuracy[key]
        for key in ("z_scored", "rounded", "binned", "strict_clean_recalculation")
    )
    assert float(error["between_center"]) == LOG_ERROR_BETWEEN_CENTER
    assert error["derivation_scope"] == "finite_primary_overt"
    assert error["transform"] == "natural_log"
    assert error["imputation"] is False
    assert error["speed_error_composite"] is False
    assert not any(config["safety"].values())


def test_exact_nine_primary_estimands_and_no_equivalence_registry() -> None:
    """Only T1--B3 are primary; beta familiarity moderation is supplementary."""

    registry = estimand_registry()
    primary = registry.loc[registry["primary_family_member"]]
    supplementary = registry.loc[~registry["primary_family_member"]]
    assert tuple(primary["estimand_id"]) == PRIMARY_ESTIMAND_IDS
    assert tuple(supplementary["estimand_id"]) == SUPPLEMENTARY_ESTIMAND_IDS
    assert len(primary) == 9
    assert primary["future_familywise_support_interval_percent"].eq(
        FAMILYWISE_SUPPORT_PERCENT
    ).all()
    assert primary["future_estimation_interval_percent"].eq(95.0).all()
    assert primary["future_fallback_holm_family"].eq("primary_nine_two_sided").all()
    assert registry["two_sided"].all()
    assert not registry["rope_authorized"].any()
    assert not registry["smallest_meaningful_db_authorized"].any()
    assert not registry["equivalence_test_authorized"].any()
    assert not registry["contains_calculated_estimate"].any()
    assert registry.loc[registry["estimand_id"].isin(["T3", "A3", "B3"]), "commonality_language_contract"].str.contains(
        "does not establish equivalence", regex=False
    ).all()


def test_expected_view_count_contract() -> None:
    """Every deterministic view count is frozen in tracked configuration."""

    config = load_config()
    assert {key: int(value) for key, value in config["expected_views"].items()} == dict(
        EXPECTED_VIEW_ROWS
    )
    assert int(config["accepted_source"]["expected_general_rows"]) == EXPECTED_GENERAL_ROWS


def test_atomic_directory_publication_and_non_destructive_archive(tmp_path: Path) -> None:
    """A staged directory publishes atomically and old output is archived."""

    target = tmp_path / "model_tables_v1"
    staging = make_staging_directory(target)
    (staging / "proof.txt").write_text("complete\n", encoding="utf-8")
    publish_staging_directory(staging, target)
    assert (target / "proof.txt").read_text(encoding="utf-8") == "complete\n"
    archived = archive_existing_output(target)
    assert archived is not None
    assert not target.exists()
    assert (archived / "proof.txt").is_file()


def test_output_namespace_is_ignored_and_broad_targets_fail() -> None:
    """Production publication is allowed only below an ignored repo path."""

    require_ignored_output_root(REPO_ROOT, PRODUCTION_ROOT)
    with pytest.raises(RuntimeError, match="unsafe_feature_output_root"):
        require_ignored_output_root(REPO_ROOT, REPO_ROOT)


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="local production model table absent")
def test_completed_production_surface_and_registries() -> None:
    """Reopen the local stage-18 production output and require all contracts."""

    manifest = json.loads(PRODUCTION_MANIFEST.read_text(encoding="utf-8"))
    validation = json.loads(
        (PRODUCTION_ROOT / "validation/model_table_validation.json").read_text(encoding="utf-8")
    )
    table = pd.read_parquet(PRODUCTION_ROOT / "model_ready_general.parquet")
    estimands = pd.read_parquet(PRODUCTION_ROOT / "registries/estimand_registry.parquet")
    assert manifest["status"] == "complete"
    assert validation["status"] == "pass"
    assert len(table) == EXPECTED_GENERAL_ROWS
    assert table["canonical_event_key"].nunique() == 8_798
    assert table["model_row_key"].nunique() == EXPECTED_GENERAL_ROWS
    assert not table.duplicated(["canonical_event_key", "model_role"]).any()
    assert set(table["model_role"]) == set(ALL_MODEL_ROLES)
    assert np.isfinite(table[["value_db", "value_log_power"]]).all().all()
    assert tuple(estimands.loc[estimands["primary_family_member"], "estimand_id"]) == PRIMARY_ESTIMAND_IDS
    assert tuple(estimands.loc[~estimands["primary_family_member"], "estimand_id"]) == SUPPLEMENTARY_ESTIMAND_IDS
    assert validation["source_value_identity"]["value_db_bitwise_equal"] is True
    assert validation["source_value_identity"]["value_log_power_bitwise_equal"] is True
    assert validation["no_model_or_inferential_output"] is True


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="local production model table absent")
def test_completed_predictor_derivations_error_nullability_and_participants() -> None:
    """Reopen persisted predictors and verify means, keys, nullability, and provenance."""

    table = pd.read_parquet(PRODUCTION_ROOT / "model_ready_general.parquet")
    trials = table.drop_duplicates("canonical_event_key").copy()
    primary = trials["scope_primary_blocks_1_5"]
    accuracy_sums = trials.loc[primary].groupby("behavioural_id")["accuracy_rating_within"].sum()
    assert float(accuracy_sums.abs().max()) < 1e-10
    assert trials["accuracy_between_center"].eq(ACCURACY_BETWEEN_CENTER).all()
    assert trials["rating_derivation_scope"].eq("primary_blocks_1_5").all()
    strict_means = trials.loc[trials["strict_clean_eligibility"]].groupby("behavioural_id")[
        "accuracy_rating_person_mean"
    ].nunique()
    assert strict_means.eq(1).all()

    valid = primary & trials["performed_condition_label"].eq("overt") & trials["objective_error_finite"]
    derived = [
        "log_error_raw", "log_error_person_mean", "log_error_within",
        "log_error_between", "log_error_between_center",
    ]
    assert int(valid.sum()) == 3_703
    assert trials.loc[valid, derived].notna().all().all()
    assert trials.loc[~valid, derived].isna().all().all()
    assert trials.loc[valid, "log_error_between_center"].eq(LOG_ERROR_BETWEEN_CENTER).all()
    missing = trials.loc[
        trials["objective_error_field_present"] & ~trials["objective_error_finite"],
        "behavioural_row_key",
    ].tolist()
    assert missing == ["91:1:6:1", "97:1:4:13"]
    assert "log_vresp_over_error" not in table.columns
    assert set(trials.loc[trials["behavioural_id"].eq(4), "canonical_event_key"]) and len(
        trials.loc[trials["behavioural_id"].eq(4)]
    ) == 48
    assert len(trials.loc[trials["behavioural_id"].eq(60)]) == 34
    p70 = trials.loc[trials["behavioural_id"].eq(70)]
    assert len(p70) == 103
    assert set(p70["original_handedness"]) == {"a"}
    assert set(p70["analysis_hand"]) == {"right"}
    assert set(p70["mapping_source"]) == {"owner_recollection"}
    assert set(p70["mapping_override"]) == {True}


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="local production model table absent")
def test_completed_views_counts_beta_pairs_and_rank() -> None:
    """Every persisted deterministic view has its exact accepted count."""

    for name, expected in EXPECTED_VIEW_ROWS.items():
        frame = pd.read_parquet(PRODUCTION_ROOT / f"views/{name}.parquet")
        assert len(frame) == expected, name
    beta = pd.read_parquet(PRODUCTION_ROOT / "views/h3_beta_primary.parquet")
    assert beta["canonical_event_key"].nunique() == 7_307
    assert beta.groupby("canonical_event_key")["phase"].nunique().eq(2).all()
    objective_theta = pd.read_parquet(PRODUCTION_ROOT / "views/objective_error_theta.parquet")
    objective_alpha = pd.read_parquet(PRODUCTION_ROOT / "views/objective_error_alpha.parquet")
    for frame in (objective_theta, objective_alpha):
        assert len(frame) == 3_703
        assert frame["behavioural_id"].nunique() == 41
        assert int(frame["strict_clean_eligibility"].sum()) == 3_699
    ranks = pd.read_parquet(PRODUCTION_ROOT / "validation/model_matrix_rank.parquet")
    assert ranks["rank"].tolist() == [6, 6, 8]
    assert ranks["full_column_rank"].all()
    assert not ranks["uses_eeg_outcome"].any()


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="local production model table absent")
def test_completed_source_immutability_and_no_forbidden_outputs() -> None:
    """Source evidence remains unchanged and no model product is published."""

    immutability = json.loads(
        (PRODUCTION_ROOT / "source_immutability.json").read_text(encoding="utf-8")
    )
    validation = json.loads(
        (PRODUCTION_ROOT / "validation/model_table_validation.json").read_text(encoding="utf-8")
    )
    assert immutability["status"] == "unchanged"
    assert immutability["feature_before"] == immutability["feature_after"]
    assert immutability["tfr_before"] == immutability["tfr_after"]
    assert immutability["behaviour_before"] == immutability["behaviour_after"]
    assert validation["no_eligibility_change"] is True
    assert validation["no_model_or_inferential_output"] is True
    forbidden = ["models", "fits", "priors", "posteriors", "contrasts", "p_values", "gamm", "csd"]
    paths = [path.relative_to(PRODUCTION_ROOT).as_posix().lower() for path in PRODUCTION_ROOT.rglob("*")]
    assert not any(token in path for token in forbidden for path in paths)


@pytest.mark.skipif(not PRODUCTION_MANIFEST.is_file(), reason="local production model table absent")
def test_completed_unchanged_reuse_does_not_rewrite_outputs() -> None:
    """Validation-only reuse hashes/reopens outputs without changing mtimes."""

    driver = EEG_DIR / "18_construct_model_ready_tables.py"
    tracked = [
        PRODUCTION_ROOT / "model_ready_general.parquet",
        PRODUCTION_ROOT / "model_table_run_manifest.json",
        PRODUCTION_ROOT / "validation/model_table_validation.json",
    ]
    before = {path: path.stat().st_mtime_ns for path in tracked}
    completed = subprocess.run(
        [sys.executable, str(driver), "--verify-current"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    after = {path: path.stat().st_mtime_ns for path in tracked}
    assert "unchanged reuse" in completed.stdout
    assert before == after
