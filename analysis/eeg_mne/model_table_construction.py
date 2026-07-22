"""Contracts and helpers for DEMI EEG model-ready tables.

This module implements the accepted stage-18 table contract. It selects five
predeclared rows per accepted trial directly from the persisted stage-17 ROI
feature table, derives the frozen within-/between-participant predictor
representations, assigns deterministic model roles and eligibility reasons,
constructs hypothesis-specific views, and validates predictor-only design
matrix rank.

Inputs:
    The tracked ``model_table_config_v1.yaml`` contract and read-only accepted
    stage-17 behavioural-lineage and ROI-feature Parquets.

Outputs:
    In-memory general/view tables and machine-readable registries used by the
    numbered stage-18 driver. Atomic file and directory helpers support safe
    publication below the ignored ``_Data/eeg/model_tables_v1`` namespace.

This module explicitly does not reconstruct EEG features from TFR arrays,
change trial eligibility, standardize or bin predictors, impute objective
error, instantiate priors, build scientific estimates, fit models, calculate
contrasts or intervals, inspect EEG associations, run sensitivities, apply
CSD, or interpret results.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml

from feature_construction import (
    atomic_write_json,
    atomic_write_parquet,
    atomic_write_text,
    file_descriptor,
    frame_descriptor,
    require_ignored_output_root,
    sha256_file,
)


MODEL_TABLE_NAMESPACE = "model_tables_v1"
MODEL_TABLE_SCHEMA_VERSION = 1
EXPECTED_TRIALS = 8_798
EXPECTED_PARTICIPANTS = 81
EXPECTED_RECORDINGS = 81
EXPECTED_PRIMARY_TRIALS = 7_307
EXPECTED_BRIDGE_TRIALS = 738
EXPECTED_STRICT_CLEAN = 8_789
EXPECTED_WARNINGS = 9
EXPECTED_GENERAL_ROWS = 43_990
EXPECTED_PRIMARY_OVERT = 3_704
EXPECTED_PRIMARY_FINITE_ERROR = 3_703
EXPECTED_ALL_OVERT = 5_195
EXPECTED_ALL_FINITE_ERROR = 5_193
EXPECTED_FILE49_TRIALS = 117
EXPECTED_P4_TRIALS = 48
EXPECTED_P60_TRIALS = 34
EXPECTED_P70_TRIALS = 103

ACCURACY_BETWEEN_CENTER = 5.70371756698893
LOG_ERROR_BETWEEN_CENTER = 4.467802672900884
PRIMARY_ESTIMAND_IDS = ("T1", "T2", "T3", "A1", "A2", "A3", "B1", "B2", "B3")
SUPPLEMENTARY_ESTIMAND_IDS = ("B4", "B5", "B6")
FAMILYWISE_SUPPORT_PERCENT = 99.4444

ACCEPTED_SOURCE_HASHES: Mapping[str, str] = {
    "feature_authority_fingerprint": "87152fb763a448354c416a0caf11bd2711a8e572181841f29c24ff33dcf06ae8",
    "feature_run_manifest_sha256": "487e0e88f2abe713654d6923a646f23ae886d41cd1a109251a53bc644f9ca066",
    "roi_features_sha256": "274e12a5fa5be04c159634288e1685d8607b27a69ac9b9f4553608553261e8ba",
    "behavioural_lineage_sha256": "5437b4ba1a6444276b8c4debfa00b1ef4539ff556ba05bd0fae2abd11fd8f6e1",
    "feature_validation_sha256": "0fedd56e4c0f5c955eadfd6daf247001df0df23e0b9f3510870de3fb0fe28bff",
    "tfr_run_manifest_sha256": "f127e80dcbc3f0292a26a81280ed79fb06e242c7a45d16292db64f77544ec48a",
    "frozen_behaviour_sha256": "3bf347da7b5007b65e0a61c135c55aa176c771cbbef49c221128eec9fc5bb831",
}

PRIMARY_MODEL_ROLES = ("theta_primary", "alpha_primary", "beta_pre", "beta_pmbr")
ALL_MODEL_ROLES = (
    "theta_primary",
    "theta_window_sensitivity",
    "alpha_primary",
    "beta_pre",
    "beta_pmbr",
)

EXPECTED_VIEW_ROWS: Mapping[str, int] = {
    "h1_theta_primary": 7_307,
    "h1_theta_window_sensitivity": 7_307,
    "h2_alpha_primary": 7_307,
    "h3_beta_primary": 14_614,
    "primary_core_union": 29_228,
    "primary_core_plus_h1_window": 36_535,
    "h1_theta_primary_strict_clean": 7_303,
    "h1_theta_window_sensitivity_strict_clean": 7_303,
    "h2_alpha_primary_strict_clean": 7_303,
    "h3_beta_primary_strict_clean": 14_606,
    "primary_core_union_strict_clean": 29_212,
    "primary_core_plus_h1_window_strict_clean": 36_515,
    "objective_error_theta": 3_703,
    "objective_error_alpha": 3_703,
    "bridge_descriptive": 2_952,
    "bridge_descriptive_plus_h1_window": 3_690,
}

EXPECTED_VIEW_STRICT_ROWS: Mapping[str, int] = {
    "h1_theta_primary": 7_303,
    "h1_theta_window_sensitivity": 7_303,
    "h2_alpha_primary": 7_303,
    "h3_beta_primary": 14_606,
    "primary_core_union": 29_212,
    "primary_core_plus_h1_window": 36_515,
    "h1_theta_primary_strict_clean": 7_303,
    "h1_theta_window_sensitivity_strict_clean": 7_303,
    "h2_alpha_primary_strict_clean": 7_303,
    "h3_beta_primary_strict_clean": 14_606,
    "primary_core_union_strict_clean": 29_212,
    "primary_core_plus_h1_window_strict_clean": 36_515,
    "objective_error_theta": 3_699,
    "objective_error_alpha": 3_699,
    "bridge_descriptive": 2_932,
    "bridge_descriptive_plus_h1_window": 3_665,
}


def utc_now() -> str:
    """Return a stable ISO-8601 UTC timestamp for manifests and archives."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def load_and_validate_model_table_config(path: Path) -> tuple[dict[str, Any], str]:
    """Load the declarative stage-18 contract and require every frozen value.

    Args:
        path: Tracked YAML configuration path.

    Returns:
        The parsed configuration and its SHA-256 digest.

    Side effects:
        Reads one YAML file. It does not write outputs.
    """

    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if config.get("schema_version") != MODEL_TABLE_SCHEMA_VERSION:
        raise RuntimeError("model_table_config_schema_mismatch")
    if config.get("pipeline_version") != MODEL_TABLE_NAMESPACE:
        raise RuntimeError("model_table_config_namespace_mismatch")
    accepted = config["accepted_source"]
    expected = {
        "accepted_trials": EXPECTED_TRIALS,
        "participants": EXPECTED_PARTICIPANTS,
        "recordings": EXPECTED_RECORDINGS,
        "primary_blocks_1_5_trials": EXPECTED_PRIMARY_TRIALS,
        "bridge_trials": EXPECTED_BRIDGE_TRIALS,
        "strict_clean_trials": EXPECTED_STRICT_CLEAN,
        "duration_warning_trials": EXPECTED_WARNINGS,
        "expected_general_rows": EXPECTED_GENERAL_ROWS,
    }
    if any(int(accepted[key]) != value for key, value in expected.items()):
        raise RuntimeError("model_table_accepted_count_contract_mismatch")
    if any(str(accepted[key]) != value for key, value in ACCEPTED_SOURCE_HASHES.items()):
        raise RuntimeError("model_table_accepted_source_hash_contract_mismatch")
    accuracy = config["predictors"]["accuracy_rating"]
    error = config["predictors"]["objective_error"]
    if (
        accuracy["derivation_scope"] != "primary_blocks_1_5"
        or not math.isclose(float(accuracy["between_center"]), ACCURACY_BETWEEN_CENTER, abs_tol=0.0)
        or accuracy["units"] != "rating_points"
        or any(bool(accuracy[key]) for key in ("z_scored", "rounded", "binned", "strict_clean_recalculation"))
    ):
        raise RuntimeError("accuracy_rating_predictor_contract_mismatch")
    if (
        error["derivation_scope"] != "finite_primary_overt"
        or error["transform"] != "natural_log"
        or not math.isclose(float(error["between_center"]), LOG_ERROR_BETWEEN_CENTER, abs_tol=0.0)
        or bool(error["imputation"])
        or bool(error["speed_error_composite"])
    ):
        raise RuntimeError("objective_error_predictor_contract_mismatch")
    roles = config["model_roles"]
    if tuple(sorted(roles, key=lambda name: int(roles[name]["order"]))) != ALL_MODEL_ROLES:
        raise RuntimeError("model_role_order_contract_mismatch")
    if tuple(config["estimands"]["primary_ids"]) != PRIMARY_ESTIMAND_IDS:
        raise RuntimeError("primary_estimand_registry_contract_mismatch")
    if tuple(config["estimands"]["supplementary_ids"]) != SUPPLEMENTARY_ESTIMAND_IDS:
        raise RuntimeError("supplementary_estimand_registry_contract_mismatch")
    if (
        int(config["estimands"]["primary_family_size"]) != 9
        or float(config["estimands"]["familywise_support_interval_percent"])
        != FAMILYWISE_SUPPORT_PERCENT
        or config["estimands"]["fallback_correction"] != "Holm"
        or config["estimands"]["tests_two_sided"] is not True
        or config["estimands"]["rope_authorized"] is not False
        or config["estimands"]["equivalence_test_authorized"] is not False
    ):
        raise RuntimeError("estimand_reporting_contract_mismatch")
    if {key: int(value) for key, value in config["expected_views"].items()} != dict(EXPECTED_VIEW_ROWS):
        raise RuntimeError("model_view_count_contract_mismatch")
    if any(bool(value) for value in config["safety"].values()):
        raise RuntimeError("model_table_config_enables_unauthorized_operation")
    storage = config["storage"]
    if (
        storage["output_root"] != "_Data/eeg/model_tables_v1"
        or storage["general_table"] != "model_ready_general.parquet"
        or storage["atomic_directory_publication"] is not True
        or storage["unchanged_validated_reuse_without_rewrite"] is not True
    ):
        raise RuntimeError("model_table_storage_contract_mismatch")
    return config, sha256_file(path)


def model_role_registry(config: Mapping[str, Any]) -> pd.DataFrame:
    """Return the exact five-role source-selection registry."""

    rows: list[dict[str, Any]] = []
    phase_codes = config["factor_coding"]["phase"]
    for name, item in sorted(config["model_roles"].items(), key=lambda pair: int(pair[1]["order"])):
        phase = item["phase"]
        rows.append(
            {
                "model_role": name,
                "role_order": int(item["order"]),
                "roi_name": item["roi_name"],
                "band": item["band"],
                "feature_definition": item["feature_definition"],
                "outcome_role": item["outcome_role"],
                "phase": phase,
                "phase_c": float(phase_codes[phase]["code"]) if phase is not None else np.nan,
                "primary_role": bool(item["primary_role"]),
                "source_is_persisted_roi_features": True,
            }
        )
    frame = pd.DataFrame(rows)
    if tuple(frame["model_role"]) != ALL_MODEL_ROLES or len(frame) != 5:
        raise RuntimeError("model_role_registry_construction_mismatch")
    return frame


def factor_coding_registry(config: Mapping[str, Any]) -> pd.DataFrame:
    """Return source labels, analysis labels, and frozen half-sum codes."""

    rows: list[dict[str, Any]] = []
    for factor, mapping in config["factor_coding"].items():
        for source_label, item in mapping.items():
            rows.append(
                {
                    "factor": factor,
                    "source_label": source_label,
                    "analysis_label": item.get("analysis_label", source_label),
                    "code": float(item["code"]),
                }
            )
    return pd.DataFrame(rows).sort_values(["factor", "code"], kind="stable").reset_index(drop=True)


def estimand_registry() -> pd.DataFrame:
    """Return the accepted nine-primary/three-supplementary estimand registry.

    The registry contains future reporting metadata only. It never contains a
    calculated EEG estimate, interval, test, or posterior quantity.
    """

    rows = [
        ("T1", "theta", "primary", "overt within-person accuracy-rating slope", "negative"),
        ("T2", "theta", "primary", "imagery within-person accuracy-rating slope", "negative_broad_direction"),
        ("T3", "theta", "primary", "imagery-minus-overt within-person accuracy-rating slope", "two_sided_no_direction"),
        ("A1", "alpha", "primary", "overt within-person accuracy-rating slope", "positive"),
        ("A2", "alpha", "primary", "imagery within-person accuracy-rating slope", "positive_broad_direction"),
        ("A3", "alpha", "primary", "imagery-minus-overt within-person accuracy-rating slope", "two_sided_no_direction"),
        ("B1", "beta", "primary", "overt PMBR-minus-pre-end contrast equally averaged over familiarity", "positive"),
        ("B2", "beta", "primary", "imagery PMBR-minus-pre-end contrast equally averaged over familiarity", "two_sided_no_direction"),
        ("B3", "beta", "primary", "imagery-minus-overt PMBR-minus-pre-end contrast difference", "two_sided_no_direction"),
        ("B4", "beta", "supplementary", "overt familiarity moderation of PMBR-minus-pre-end", "positive"),
        ("B5", "beta", "supplementary", "imagery familiarity moderation of PMBR-minus-pre-end", "two_sided_no_direction"),
        ("B6", "beta", "supplementary", "imagery-minus-overt difference in familiarity moderation", "two_sided_no_direction"),
    ]
    commonality = (
        "Direct condition differences are two-sided; an interval crossing zero does not establish "
        "equivalence. Analogous or parallel effects may be described when supported, but strict "
        "equivalence is not authorized."
    )
    records = []
    for identifier, family, classification, conceptual, direction in rows:
        primary = classification == "primary"
        records.append(
            {
                "estimand_id": identifier,
                "outcome_family": family,
                "classification": classification,
                "primary_family_member": primary,
                "primary_family_size": 9,
                "conceptual_contrast": conceptual,
                "historical_direction": direction,
                "two_sided": True,
                "future_estimation_interval_percent": 95.0 if primary else np.nan,
                "future_familywise_support_interval_percent": (
                    FAMILYWISE_SUPPORT_PERCENT if primary else np.nan
                ),
                "future_fallback_holm_family": "primary_nine_two_sided" if primary else pd.NA,
                "commonality_language_contract": commonality,
                "rope_authorized": False,
                "smallest_meaningful_db_authorized": False,
                "equivalence_test_authorized": False,
                "contains_calculated_estimate": False,
            }
        )
    frame = pd.DataFrame(records)
    primary_ids = tuple(frame.loc[frame["primary_family_member"], "estimand_id"])
    supplementary_ids = tuple(frame.loc[~frame["primary_family_member"], "estimand_id"])
    if primary_ids != PRIMARY_ESTIMAND_IDS or supplementary_ids != SUPPLEMENTARY_ESTIMAND_IDS:
        raise RuntimeError("estimand_registry_classification_mismatch")
    return frame


def controlled_reason_registry() -> pd.DataFrame:
    """Return the only authorized inclusion/exclusion reason codes."""

    rows = [
        ("primary", "eligible_primary_scope_role", "inclusion"),
        ("primary", "not_primary_scope", "exclusion"),
        ("primary", "sensitivity_role_not_primary", "exclusion"),
        ("objective_error", "eligible_finite_primary_overt_error", "inclusion"),
        ("objective_error", "wrong_model_role", "exclusion"),
        ("objective_error", "not_primary_scope", "exclusion"),
        ("objective_error", "not_overt_condition", "exclusion"),
        ("objective_error", "objective_error_missing_source_value", "exclusion"),
        ("bridge", "eligible_bridge_scope_role", "inclusion"),
        ("bridge", "not_bridge_scope", "exclusion"),
    ]
    return pd.DataFrame(rows, columns=["analysis_surface", "reason_code", "reason_type"])


def derive_trial_predictors(
    lineage: pd.DataFrame, config: Mapping[str, Any]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Derive frozen trial-level predictor and factor fields.

    Accuracy participant means use every primary blocks-1--5 trial, including
    warning rows. Objective-error means use finite primary-overt values only.
    Derived error fields remain null everywhere else.

    Args:
        lineage: Accepted stage-17 one-row-per-trial behavioural lineage.
        config: Validated stage-18 configuration.

    Returns:
        A one-row-per-trial derivation frame and exact validation facts.
    """

    required = {
        "canonical_event_key", "behavioural_row_key", "behavioural_id", "recording_stem",
        "group", "performed_condition", "familiarity_repetition_condition",
        "accuracy_rating", "error", "scope_primary_blocks_1_5",
        "scope_imagery_final_overt_bridge", "strict_clean_eligibility",
        "duration_warning_flag",
    }
    missing = sorted(required - set(lineage.columns))
    if missing:
        raise RuntimeError(f"predictor_lineage_columns_missing:{missing}")
    if len(lineage) != EXPECTED_TRIALS or lineage["canonical_event_key"].duplicated().any():
        raise RuntimeError("predictor_lineage_identity_mismatch")

    frame = lineage.copy()
    frame["accuracy_rating_raw"] = frame["accuracy_rating"].astype(np.int64)
    if frame["accuracy_rating_raw"].isna().any() or not frame["accuracy_rating_raw"].between(1, 10).all():
        raise RuntimeError("accuracy_rating_raw_domain_mismatch")
    primary = frame["scope_primary_blocks_1_5"]
    primary_means = (
        frame.loc[primary]
        .groupby("behavioural_id", sort=True)["accuracy_rating_raw"]
        .mean()
    )
    if len(primary_means) != EXPECTED_PARTICIPANTS:
        raise RuntimeError("accuracy_participant_mean_count_mismatch")
    observed_accuracy_center = float(primary_means.mean())
    if not math.isclose(observed_accuracy_center, ACCURACY_BETWEEN_CENTER, abs_tol=1e-12, rel_tol=0.0):
        raise RuntimeError(f"accuracy_between_center_mismatch:{observed_accuracy_center}")
    frame["accuracy_rating_person_mean"] = frame["behavioural_id"].map(primary_means)
    frame["accuracy_rating_within"] = (
        frame["accuracy_rating_raw"] - frame["accuracy_rating_person_mean"]
    )
    frame["accuracy_rating_between"] = (
        frame["accuracy_rating_person_mean"] - ACCURACY_BETWEEN_CENTER
    )
    frame["rating_derivation_scope"] = "primary_blocks_1_5"
    frame["accuracy_between_center"] = ACCURACY_BETWEEN_CENTER
    primary_within_sums = (
        frame.loc[primary].groupby("behavioural_id")["accuracy_rating_within"].sum()
    )
    max_accuracy_sum = float(primary_within_sums.abs().max())
    if max_accuracy_sum >= 1e-10:
        raise RuntimeError(f"accuracy_within_sum_not_zero:{max_accuracy_sum}")

    frame["objective_error_field_present"] = frame["performed_condition"].eq("physical")
    frame["objective_error_finite"] = (
        frame["objective_error_field_present"] & np.isfinite(frame["error"])
    )
    frame["objective_error_validity_status"] = np.select(
        [
            frame["objective_error_finite"],
            frame["objective_error_field_present"],
        ],
        ["finite_usable_source_value", "missing_source_value"],
        default="not_available_imagery",
    )
    valid_error = primary & frame["performed_condition"].eq("physical") & frame["objective_error_finite"]
    if not frame.loc[valid_error, "error"].gt(0).all():
        raise RuntimeError("finite_primary_overt_error_not_positive")
    frame["log_error_raw"] = np.nan
    frame.loc[valid_error, "log_error_raw"] = np.log(frame.loc[valid_error, "error"])
    error_means = (
        frame.loc[valid_error]
        .groupby("behavioural_id", sort=True)["log_error_raw"]
        .mean()
    )
    observed_error_center = float(error_means.mean())
    if len(error_means) != 41 or not math.isclose(
        observed_error_center, LOG_ERROR_BETWEEN_CENTER, abs_tol=1e-12, rel_tol=0.0
    ):
        raise RuntimeError(f"log_error_between_center_mismatch:{len(error_means)}:{observed_error_center}")
    frame["log_error_person_mean"] = np.nan
    frame.loc[valid_error, "log_error_person_mean"] = frame.loc[valid_error, "behavioural_id"].map(error_means)
    frame["log_error_within"] = np.nan
    frame.loc[valid_error, "log_error_within"] = (
        frame.loc[valid_error, "log_error_raw"]
        - frame.loc[valid_error, "log_error_person_mean"]
    )
    frame["log_error_between"] = np.nan
    frame.loc[valid_error, "log_error_between"] = (
        frame.loc[valid_error, "log_error_person_mean"] - LOG_ERROR_BETWEEN_CENTER
    )
    frame["error_derivation_scope"] = pd.NA
    frame.loc[valid_error, "error_derivation_scope"] = "finite_primary_overt"
    frame["log_error_between_center"] = np.nan
    frame.loc[valid_error, "log_error_between_center"] = LOG_ERROR_BETWEEN_CENTER
    max_error_sum = float(
        frame.loc[valid_error].groupby("behavioural_id")["log_error_within"].sum().abs().max()
    )
    if max_error_sum >= 1e-10:
        raise RuntimeError(f"log_error_within_sum_not_zero:{max_error_sum}")

    condition_map = config["factor_coding"]["condition"]
    familiarity_map = config["factor_coding"]["familiarity"]
    frame["assigned_condition_source"] = frame["group"]
    frame["assigned_condition"] = frame["group"].map(
        {label: item["analysis_label"] for label, item in condition_map.items()}
    )
    frame["performed_condition_source"] = frame["performed_condition"]
    frame["performed_condition_label"] = frame["performed_condition"].map(
        {label: item["analysis_label"] for label, item in condition_map.items()}
    )
    frame["condition_c"] = frame["performed_condition"].map(
        {label: float(item["code"]) for label, item in condition_map.items()}
    )
    frame["familiarity_source"] = frame["familiarity_repetition_condition"]
    frame["familiarity"] = frame["familiarity_repetition_condition"].map(
        {label: item["analysis_label"] for label, item in familiarity_map.items()}
    )
    frame["familiarity_c"] = frame["familiarity_repetition_condition"].map(
        {label: float(item["code"]) for label, item in familiarity_map.items()}
    )
    if frame[["assigned_condition", "performed_condition_label", "condition_c", "familiarity", "familiarity_c"]].isna().any().any():
        raise RuntimeError("factor_coding_source_domain_mismatch")

    missing_error_rows = frame.loc[
        frame["objective_error_field_present"] & ~frame["objective_error_finite"],
        "behavioural_row_key",
    ].tolist()
    if missing_error_rows != ["91:1:6:1", "97:1:4:13"]:
        raise RuntimeError(f"objective_error_missing_keys_mismatch:{missing_error_rows}")
    if frame.loc[valid_error & frame["behavioural_row_key"].eq("97:1:4:13")].any().any():
        raise RuntimeError("primary_missing_error_row_marked_valid")

    facts = {
        "primary_trials": int(primary.sum()),
        "accuracy_participants": len(primary_means),
        "accuracy_between_center": observed_accuracy_center,
        "maximum_absolute_accuracy_within_sum": max_accuracy_sum,
        "all_overt_rows": int(frame["objective_error_field_present"].sum()),
        "all_finite_error_rows": int(frame["objective_error_finite"].sum()),
        "primary_overt_rows": int((primary & frame["performed_condition"].eq("physical")).sum()),
        "primary_finite_error_rows": int(valid_error.sum()),
        "log_error_participants": len(error_means),
        "log_error_between_center": observed_error_center,
        "maximum_absolute_log_error_within_sum": max_error_sum,
        "missing_error_keys": missing_error_rows,
    }
    expected_facts = {
        "primary_trials": EXPECTED_PRIMARY_TRIALS,
        "all_overt_rows": EXPECTED_ALL_OVERT,
        "all_finite_error_rows": EXPECTED_ALL_FINITE_ERROR,
        "primary_overt_rows": EXPECTED_PRIMARY_OVERT,
        "primary_finite_error_rows": EXPECTED_PRIMARY_FINITE_ERROR,
    }
    if any(int(facts[key]) != value for key, value in expected_facts.items()):
        raise RuntimeError(f"predictor_surface_count_mismatch:{facts}")

    columns = [
        "canonical_event_key", "accuracy_rating_raw", "accuracy_rating_person_mean",
        "accuracy_rating_within", "accuracy_rating_between", "rating_derivation_scope",
        "accuracy_between_center", "error", "objective_error_field_present",
        "objective_error_finite", "objective_error_validity_status", "log_error_raw",
        "log_error_person_mean", "log_error_within", "log_error_between",
        "error_derivation_scope", "log_error_between_center", "assigned_condition_source",
        "assigned_condition", "performed_condition_source", "performed_condition_label",
        "condition_c", "familiarity_source", "familiarity", "familiarity_c",
    ]
    return frame[columns].copy(), facts


def select_required_roi_rows(
    roi_features: pd.DataFrame, role_registry: pd.DataFrame
) -> pd.DataFrame:
    """Select exactly five persisted ROI rows per accepted trial.

    EEG values are copied directly from ``roi_features`` without arithmetic or
    reconstruction. Source row identities and values are retained for later
    equality validation.
    """

    parts: list[pd.DataFrame] = []
    for role in role_registry.itertuples(index=False):
        selected = roi_features.loc[
            roi_features["roi_name"].eq(role.roi_name)
            & roi_features["feature_definition"].eq(role.feature_definition)
            & roi_features["band"].eq(role.band)
        ].copy()
        if (
            len(selected) != EXPECTED_TRIALS
            or selected["canonical_event_key"].nunique() != EXPECTED_TRIALS
            or selected.duplicated(["canonical_event_key", "roi_name", "feature_definition", "band"]).any()
            or not np.isfinite(selected[["value_db", "value_log_power"]]).all().all()
        ):
            raise RuntimeError(f"required_roi_role_source_mismatch:{role.model_role}")
        selected["model_role"] = role.model_role
        selected["role_order"] = int(role.role_order)
        selected["outcome_role"] = role.outcome_role
        selected["phase"] = role.phase
        selected["phase_c"] = role.phase_c
        selected["source_roi_row_identity"] = (
            selected["canonical_event_key"].astype(str)
            + ":" + selected["roi_name"].astype(str)
            + ":" + selected["feature_definition"].astype(str)
            + ":" + selected["band"].astype(str)
        )
        parts.append(selected)
    result = pd.concat(parts, ignore_index=True)
    result = result.sort_values(["canonical_order_index", "role_order"], kind="stable").reset_index(drop=True)
    if len(result) != EXPECTED_GENERAL_ROWS or result.duplicated(["canonical_event_key", "model_role"]).any():
        raise RuntimeError("required_roi_general_identity_mismatch")
    return result


def add_eligibility_contract(table: pd.DataFrame) -> pd.DataFrame:
    """Add controlled primary, objective-error, and bridge eligibility fields."""

    result = table.copy()
    primary_scope = result["scope_primary_blocks_1_5"].astype(bool)
    primary_role = result["model_role"].isin(PRIMARY_MODEL_ROLES)
    result["primary_eligible"] = primary_scope & primary_role
    result["primary_inclusion_reason"] = pd.NA
    result.loc[result["primary_eligible"], "primary_inclusion_reason"] = "eligible_primary_scope_role"
    result["primary_exclusion_reason"] = pd.NA
    result.loc[~primary_scope, "primary_exclusion_reason"] = "not_primary_scope"
    result.loc[primary_scope & ~primary_role, "primary_exclusion_reason"] = "sensitivity_role_not_primary"

    objective_role = result["model_role"].isin(["theta_primary", "alpha_primary"])
    objective_eligible = (
        objective_role
        & primary_scope
        & result["performed_condition_label"].eq("overt")
        & result["objective_error_finite"]
    )
    result["objective_error_eligible"] = objective_eligible
    result["objective_error_inclusion_reason"] = pd.NA
    result.loc[objective_eligible, "objective_error_inclusion_reason"] = "eligible_finite_primary_overt_error"
    result["objective_error_exclusion_reason"] = pd.NA
    result.loc[~objective_role, "objective_error_exclusion_reason"] = "wrong_model_role"
    result.loc[objective_role & ~primary_scope, "objective_error_exclusion_reason"] = "not_primary_scope"
    result.loc[
        objective_role & primary_scope & ~result["performed_condition_label"].eq("overt"),
        "objective_error_exclusion_reason",
    ] = "not_overt_condition"
    result.loc[
        objective_role & primary_scope & result["performed_condition_label"].eq("overt")
        & ~result["objective_error_finite"],
        "objective_error_exclusion_reason",
    ] = "objective_error_missing_source_value"

    bridge = result["scope_imagery_final_overt_bridge"].astype(bool)
    result["bridge_eligible"] = bridge
    result["bridge_inclusion_reason"] = pd.NA
    result.loc[bridge, "bridge_inclusion_reason"] = "eligible_bridge_scope_role"
    result["bridge_exclusion_reason"] = pd.NA
    result.loc[~bridge, "bridge_exclusion_reason"] = "not_bridge_scope"

    result["eligible_h1_primary"] = primary_scope & result["model_role"].eq("theta_primary")
    result["eligible_h1_window_sensitivity"] = primary_scope & result["model_role"].eq("theta_window_sensitivity")
    result["eligible_h2_primary"] = primary_scope & result["model_role"].eq("alpha_primary")
    result["eligible_h3_primary"] = primary_scope & result["model_role"].isin(["beta_pre", "beta_pmbr"])
    result["eligible_objective_error_secondary"] = objective_eligible
    result["eligible_bridge_descriptive"] = bridge
    return result


def build_general_table(
    lineage: pd.DataFrame,
    roi_features: pd.DataFrame,
    config: Mapping[str, Any],
    *,
    provenance: Mapping[str, str],
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Build the one-five-role-row authoritative general table in memory.

    Args:
        lineage: Accepted one-row-per-trial behavioural lineage.
        roi_features: Accepted persisted stage-17 ROI table.
        config: Validated stage-18 configuration.
        provenance: Frozen source/config/code hash strings to repeat per row.

    Returns:
        General table, predictor-derivation facts, and selected source ROI rows.
    """

    roles = model_role_registry(config)
    selected = select_required_roi_rows(roi_features, roles)
    predictors, predictor_facts = derive_trial_predictors(lineage, config)
    source_accuracy = selected["accuracy_rating"].to_numpy(copy=True)
    source_error = selected["error"].to_numpy(copy=True)
    selected = selected.drop(columns=["accuracy_rating", "error"])
    table = selected.merge(predictors, on="canonical_event_key", how="left", validate="many_to_one")
    if not np.array_equal(source_accuracy, table["accuracy_rating_raw"].to_numpy()):
        raise RuntimeError("roi_lineage_accuracy_identity_mismatch")
    source_error_finite = np.isfinite(source_error)
    table_error = table["error"].to_numpy(dtype=float)
    if not np.array_equal(source_error_finite, np.isfinite(table_error)) or not np.array_equal(
        source_error[source_error_finite], table_error[source_error_finite]
    ):
        raise RuntimeError("roi_lineage_error_identity_mismatch")

    table["model_table_contract_version"] = MODEL_TABLE_NAMESPACE
    table["model_row_key"] = table["canonical_event_key"].astype(str) + ":" + table["model_role"]
    table["participant_id"] = table["behavioural_id"]
    table["session"] = table["offset_session"]
    table["block"] = table["offset_block"]
    table["trial"] = table["offset_trial"]
    table["analysis_scope"] = np.select(
        [table["scope_primary_blocks_1_5"], table["scope_imagery_final_overt_bridge"]],
        ["primary_blocks_1_5", "imagery_final_overt_bridge"],
        default="accepted_nonprimary_other",
    )
    for key, value in provenance.items():
        table[key] = value
    table = add_eligibility_contract(table)
    table = table.sort_values(["canonical_order_index", "role_order"], kind="stable").reset_index(drop=True)
    if len(table) != EXPECTED_GENERAL_ROWS or table["model_row_key"].duplicated().any():
        raise RuntimeError("general_table_identity_mismatch")
    return table, predictor_facts, selected


def deterministic_views(table: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Return every accepted hypothesis/scope view with stable row order."""

    primary = table["scope_primary_blocks_1_5"]
    strict = table["strict_clean_eligibility"]
    role = table["model_role"]
    core = role.isin(PRIMARY_MODEL_ROLES)
    plus_window = role.isin(ALL_MODEL_ROLES)
    bridge = table["scope_imagery_final_overt_bridge"]
    views = {
        "h1_theta_primary": table.loc[primary & role.eq("theta_primary")],
        "h1_theta_window_sensitivity": table.loc[primary & role.eq("theta_window_sensitivity")],
        "h2_alpha_primary": table.loc[primary & role.eq("alpha_primary")],
        "h3_beta_primary": table.loc[primary & role.isin(["beta_pre", "beta_pmbr"])],
        "primary_core_union": table.loc[primary & core],
        "primary_core_plus_h1_window": table.loc[primary & plus_window],
        "h1_theta_primary_strict_clean": table.loc[primary & strict & role.eq("theta_primary")],
        "h1_theta_window_sensitivity_strict_clean": table.loc[primary & strict & role.eq("theta_window_sensitivity")],
        "h2_alpha_primary_strict_clean": table.loc[primary & strict & role.eq("alpha_primary")],
        "h3_beta_primary_strict_clean": table.loc[primary & strict & role.isin(["beta_pre", "beta_pmbr"])],
        "primary_core_union_strict_clean": table.loc[primary & strict & core],
        "primary_core_plus_h1_window_strict_clean": table.loc[primary & strict & plus_window],
        "objective_error_theta": table.loc[table["objective_error_eligible"] & role.eq("theta_primary")],
        "objective_error_alpha": table.loc[table["objective_error_eligible"] & role.eq("alpha_primary")],
        "bridge_descriptive": table.loc[bridge & core],
        "bridge_descriptive_plus_h1_window": table.loc[bridge & plus_window],
    }
    result = {name: frame.copy().reset_index(drop=True) for name, frame in views.items()}
    observed = {name: len(frame) for name, frame in result.items()}
    if observed != dict(EXPECTED_VIEW_ROWS):
        raise RuntimeError(f"deterministic_view_count_mismatch:{observed}")
    return result


def model_matrix_rank_validation(views: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Validate predictor-only fixed-effect ranks without using EEG outcomes."""

    rows: list[dict[str, Any]] = []
    for name, view_name in (("h1_theta", "h1_theta_primary"), ("h2_alpha", "h2_alpha_primary")):
        frame = views[view_name]
        condition = frame["condition_c"].to_numpy(dtype=float)
        within = frame["accuracy_rating_within"].to_numpy(dtype=float)
        between = frame["accuracy_rating_between"].to_numpy(dtype=float)
        familiarity = frame["familiarity_c"].to_numpy(dtype=float)
        matrix = np.column_stack(
            [np.ones(len(frame)), condition, within, condition * within, between, familiarity]
        )
        rank = int(np.linalg.matrix_rank(matrix))
        rows.append(
            {
                "model": name,
                "rows": len(frame),
                "columns": matrix.shape[1],
                "rank": rank,
                "full_column_rank": rank == matrix.shape[1],
                "fixed_effect_contract": "1 + condition_c * accuracy_rating_within + accuracy_rating_between + familiarity_c",
                "uses_eeg_outcome": False,
            }
        )
    beta = views["h3_beta_primary"]
    c = beta["condition_c"].to_numpy(dtype=float)
    p = beta["phase_c"].to_numpy(dtype=float)
    f = beta["familiarity_c"].to_numpy(dtype=float)
    beta_matrix = np.column_stack(
        [np.ones(len(beta)), c, p, f, c * p, c * f, p * f, c * p * f]
    )
    beta_rank = int(np.linalg.matrix_rank(beta_matrix))
    rows.append(
        {
            "model": "h3_beta",
            "rows": len(beta),
            "columns": beta_matrix.shape[1],
            "rank": beta_rank,
            "full_column_rank": beta_rank == beta_matrix.shape[1],
            "fixed_effect_contract": "1 + condition_c * phase_c * familiarity_c",
            "uses_eeg_outcome": False,
        }
    )
    result = pd.DataFrame(rows)
    if not result["full_column_rank"].all() or result["rank"].tolist() != [6, 6, 8]:
        raise RuntimeError(f"model_matrix_rank_validation_failed:{result.to_dict('records')}")
    return result


def source_value_identity(
    table: pd.DataFrame, selected_source_rows: pd.DataFrame
) -> dict[str, Any]:
    """Require persisted table EEG values to equal selected source ROI values."""

    identity = ["canonical_event_key", "model_role"]
    source = selected_source_rows[
        [*identity, "value_db", "value_log_power", "source_roi_row_identity"]
    ].rename(
        columns={
            "value_db": "source_value_db",
            "value_log_power": "source_value_log_power",
            "source_roi_row_identity": "expected_source_roi_row_identity",
        }
    )
    linked = table[
        [*identity, "value_db", "value_log_power", "source_roi_row_identity"]
    ].merge(source, on=identity, how="left", validate="one_to_one")
    db_equal = np.array_equal(linked["value_db"].to_numpy(), linked["source_value_db"].to_numpy())
    log_equal = np.array_equal(
        linked["value_log_power"].to_numpy(), linked["source_value_log_power"].to_numpy()
    )
    identity_equal = linked["source_roi_row_identity"].equals(
        linked["expected_source_roi_row_identity"]
    )
    if not (db_equal and log_equal and identity_equal):
        raise RuntimeError("persisted_roi_value_identity_failure")
    return {
        "rows_checked": len(linked),
        "value_db_bitwise_equal": db_equal,
        "value_log_power_bitwise_equal": log_equal,
        "source_roi_row_identity_equal": identity_equal,
    }


def validate_general_and_views(
    table: pd.DataFrame,
    views: Mapping[str, pd.DataFrame],
    selected_source_rows: pd.DataFrame,
    predictor_facts: Mapping[str, Any],
    estimands: pd.DataFrame,
    rank_validation: pd.DataFrame,
) -> dict[str, Any]:
    """Run the complete neutral stage-18 in-memory validation contract."""

    required_nonmissing = [
        "canonical_event_key", "model_row_key", "participant_id", "eeg_source_id",
        "recording_stem", "session", "block", "trial", "assigned_condition",
        "performed_condition_label", "condition_c", "familiarity", "familiarity_c",
        "accuracy_rating_raw", "accuracy_rating_person_mean", "accuracy_rating_within",
        "accuracy_rating_between", "roi_name", "band", "feature_definition",
        "window_minimum_seconds", "window_maximum_seconds", "value_db",
        "value_log_power", "original_handedness", "analysis_hand", "mapping_source",
    ]
    if table[required_nonmissing].isna().any().any():
        raise RuntimeError("general_table_required_value_missing")
    allowed_error = table["scope_primary_blocks_1_5"] & table["performed_condition_label"].eq("overt") & table["objective_error_finite"]
    derived_error_columns = [
        "log_error_raw", "log_error_person_mean", "log_error_within",
        "log_error_between", "log_error_between_center",
    ]
    if table.loc[allowed_error, derived_error_columns].isna().any().any():
        raise RuntimeError("valid_error_scope_derived_value_missing")
    if table.loc[~allowed_error, derived_error_columns].notna().any().any():
        raise RuntimeError("derived_error_value_present_outside_valid_scope")
    if table.loc[~allowed_error, "error_derivation_scope"].notna().any():
        raise RuntimeError("error_derivation_scope_present_outside_valid_scope")
    if set(table["condition_c"]) != {-0.5, 0.5} or set(table["familiarity_c"]) != {-0.5, 0.5}:
        raise RuntimeError("general_table_factor_domain_mismatch")
    beta = table["model_role"].isin(["beta_pre", "beta_pmbr"])
    if set(table.loc[beta, "phase_c"]) != {-0.5, 0.5} or table.loc[~beta, "phase_c"].notna().any():
        raise RuntimeError("general_table_phase_domain_mismatch")
    phase_counts = views["h3_beta_primary"].groupby("canonical_event_key")["phase"].nunique()
    if len(phase_counts) != EXPECTED_PRIMARY_TRIALS or not phase_counts.eq(2).all():
        raise RuntimeError("h3_two_phase_per_trial_mismatch")
    if table["canonical_event_key"].nunique() != EXPECTED_TRIALS:
        raise RuntimeError("accepted_trial_representation_mismatch")
    participant_counts = table.drop_duplicates("canonical_event_key").groupby("behavioural_id").size()
    if (
        int(participant_counts.get(4, 0)) != EXPECTED_P4_TRIALS
        or int(participant_counts.get(60, 0)) != EXPECTED_P60_TRIALS
        or int(participant_counts.get(70, 0)) != EXPECTED_P70_TRIALS
    ):
        raise RuntimeError("accepted_low_count_or_override_participant_missing")
    p70 = table.loc[table["behavioural_id"].eq(70)]
    if not (
        set(p70["original_handedness"]) == {"a"}
        and set(p70["analysis_hand"]) == {"right"}
        and set(p70["mapping_source"]) == {"owner_recollection"}
        and set(p70["mapping_override"]) == {True}
    ):
        raise RuntimeError("participant_70_model_table_provenance_mismatch")
    trial_surface = table.drop_duplicates("canonical_event_key")
    if int(trial_surface["eeg_source_id"].eq(49).sum()) != EXPECTED_FILE49_TRIALS:
        raise RuntimeError("file49_trial_count_mismatch")
    filenames = trial_surface["source_recording_filename"].astype(str)
    if trial_surface["eeg_source_id"].eq(86).any() or filenames.str.contains("54_1", regex=False).any():
        raise RuntimeError("file54_1_or_id86_unexpectedly_present")
    value_identity = source_value_identity(table, selected_source_rows)
    view_strict_rows = {
        name: int(frame["strict_clean_eligibility"].sum())
        for name, frame in views.items()
    }
    if view_strict_rows != dict(EXPECTED_VIEW_STRICT_ROWS):
        raise RuntimeError(f"deterministic_view_strict_count_mismatch:{view_strict_rows}")
    primary_estimands = estimands.loc[estimands["primary_family_member"]]
    supplementary = estimands.loc[~estimands["primary_family_member"]]
    if (
        tuple(primary_estimands["estimand_id"]) != PRIMARY_ESTIMAND_IDS
        or tuple(supplementary["estimand_id"]) != SUPPLEMENTARY_ESTIMAND_IDS
        or estimands["rope_authorized"].any()
        or estimands["equivalence_test_authorized"].any()
        or not primary_estimands["future_familywise_support_interval_percent"].eq(FAMILYWISE_SUPPORT_PERCENT).all()
    ):
        raise RuntimeError("estimand_registry_validation_failed")
    validation = {
        "status": "pass",
        "general_rows": len(table),
        "accepted_trials": table["canonical_event_key"].nunique(),
        "participants": table["behavioural_id"].nunique(),
        "recordings": table["recording_stem"].nunique(),
        "strict_clean_trials": int(trial_surface["strict_clean_eligibility"].sum()),
        "duration_warning_trials": int(trial_surface["duration_warning_flag"].sum()),
        "view_rows": {name: len(frame) for name, frame in views.items()},
        "view_strict_clean_rows": view_strict_rows,
        "predictor_derivations": dict(predictor_facts),
        "source_value_identity": value_identity,
        "model_matrix_rank": rank_validation.to_dict("records"),
        "nine_primary_estimands": primary_estimands["estimand_id"].tolist(),
        "supplementary_estimands": supplementary["estimand_id"].tolist(),
        "familywise_support_interval_percent": FAMILYWISE_SUPPORT_PERCENT,
        "no_rope_or_equivalence_test": bool(
            (~estimands["rope_authorized"] & ~estimands["equivalence_test_authorized"]).all()
        ),
        "participant_70_provenance": {
            "original_handedness": "a",
            "analysis_hand": "right",
            "mapping_source": "owner_recollection",
            "mapping_override": True,
        },
        "participants_4_60_70_retained": True,
        "file49_retained": True,
        "file54_1_absent": True,
        "id86_absent": True,
        "no_eligibility_change": True,
        "no_model_or_inferential_output": True,
    }
    expected_top = {
        "general_rows": EXPECTED_GENERAL_ROWS,
        "accepted_trials": EXPECTED_TRIALS,
        "participants": EXPECTED_PARTICIPANTS,
        "recordings": EXPECTED_RECORDINGS,
        "strict_clean_trials": EXPECTED_STRICT_CLEAN,
        "duration_warning_trials": EXPECTED_WARNINGS,
    }
    if any(int(validation[key]) != value for key, value in expected_top.items()):
        raise RuntimeError(f"general_table_validation_count_mismatch:{validation}")
    return validation


def count_summary(table: pd.DataFrame, views: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Return compact exact row/trial/participant counts for every output view."""

    frames: list[tuple[str, pd.DataFrame]] = [("general_table", table), *views.items()]
    rows = []
    for name, frame in frames:
        rows.append(
            {
                "surface": name,
                "rows": len(frame),
                "unique_trials": frame["canonical_event_key"].nunique(),
                "participants": frame["behavioural_id"].nunique(),
                "recordings": frame["recording_stem"].nunique(),
                "strict_clean_rows": int(frame["strict_clean_eligibility"].sum()),
                "duration_warning_rows": int(frame["duration_warning_flag"].sum()),
            }
        )
    return pd.DataFrame(rows)


def make_staging_directory(output_root: Path) -> Path:
    """Create a sibling staging directory for atomic whole-stage publication."""

    output_root.parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=f".{output_root.name}.tmp-", dir=output_root.parent))


def archive_existing_output(output_root: Path) -> Path | None:
    """Move a non-reusable existing output aside without deleting it."""

    if not output_root.exists():
        return None
    history = output_root.parent / f"{output_root.name}_history"
    history.mkdir(parents=True, exist_ok=True)
    stamp = utc_now().replace(":", "").replace("-", "").replace(".", "")
    destination = history / f"stale_{stamp}"
    os.replace(output_root, destination)
    return destination


def publish_staging_directory(staging: Path, output_root: Path) -> None:
    """Atomically publish a validated staging directory at the final path."""

    if output_root.exists():
        raise RuntimeError("atomic_model_table_target_already_exists")
    if staging.parent.resolve() != output_root.parent.resolve():
        raise RuntimeError("atomic_model_table_staging_not_sibling")
    os.replace(staging, output_root)


__all__ = [
    "ACCURACY_BETWEEN_CENTER", "ALL_MODEL_ROLES", "EXPECTED_GENERAL_ROWS",
    "EXPECTED_VIEW_ROWS", "FAMILYWISE_SUPPORT_PERCENT", "LOG_ERROR_BETWEEN_CENTER",
    "MODEL_TABLE_NAMESPACE", "MODEL_TABLE_SCHEMA_VERSION", "PRIMARY_ESTIMAND_IDS",
    "SUPPLEMENTARY_ESTIMAND_IDS", "archive_existing_output", "atomic_write_json",
    "atomic_write_parquet", "atomic_write_text", "build_general_table", "count_summary",
    "controlled_reason_registry", "deterministic_views", "estimand_registry",
    "factor_coding_registry", "file_descriptor", "frame_descriptor",
    "load_and_validate_model_table_config", "make_staging_directory",
    "model_matrix_rank_validation", "model_role_registry", "publish_staging_directory",
    "require_ignored_output_root", "sha256_file", "source_value_identity",
    "validate_general_and_views",
]
