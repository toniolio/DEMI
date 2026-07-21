"""Contracts and helpers for DEMI trial-level band/window EEG features.

This module implements the accepted stage-17 feature contract. It validates
the declared 4--30 Hz bands and fixed response-onset/response-end windows,
computes arithmetic means of already transformed TFR cells, derives the
required unnormalized log-power sensitivity cell-by-cell, retains all physical
channel features, and derives equal-weight predeclared ROIs from those channel
rows.

Inputs:
    The tracked ``feature_config_v1.yaml`` contract, accepted ``tfr_v1`` axes,
    read-only response-onset/response-end raw and dB arrays, linked frozen
    behavioral metadata, and standard_1005 coordinates.

Outputs:
    In-memory channel and ROI feature frames plus deterministic definitions,
    hashes, validation facts, and atomic Parquet/JSON helpers used by the
    numbered production driver.

This module explicitly does not convolve data, create a baseline, use 31--40
Hz for features, split beta, apply CSD, average participants, change trial
eligibility, center predictors, build model matrices, fit models, compute
contrasts, or interpret effects.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping, Sequence
import uuid

import mne
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from epoch_construction import sha256_file
from event_epoch_eligibility import stable_frame_hash


FEATURE_NAMESPACE = "features_v1"
FEATURE_SCHEMA_VERSION = 1
EXPECTED_TRIAL_COUNT = 8_798
EXPECTED_PARTICIPANT_COUNT = 81
EXPECTED_RECORDING_COUNT = 81
EXPECTED_STRICT_CLEAN_COUNT = 8_789
EXPECTED_DURATION_WARNING_COUNT = 9
EXPECTED_FILE49_COUNT = 117
EXPECTED_SHORT_0P5_COUNT = 47
EXPECTED_SHORT_1P0_COUNT = 681
EXPECTED_CHANNEL_COUNT = 30
EXPECTED_CHANNEL_DEFINITIONS = 11
EXPECTED_ROI_DEFINITIONS = 13
EXPECTED_CHANNEL_ROWS = 2_903_340
EXPECTED_ROI_ROWS = 114_374
EXPECTED_P4_TRIALS = 48
EXPECTED_P60_TRIALS = 34

PRIMARY_CHANNELS: tuple[str, ...] = (
    "FP1", "FP2", "F7", "F3", "FZ", "F4", "F8", "FT7", "FC3", "FCZ",
    "FC4", "FT8", "T7", "C3", "CZ", "C4", "T8", "TP7", "CP3", "CPZ",
    "CP4", "TP8", "P7", "P3", "PZ", "P4", "P8", "O1", "OZ", "O2",
)
BAND_LIMITS: Mapping[str, tuple[float, float, int]] = {
    "theta": (4.0, 8.0, 5),
    "alpha": (9.0, 12.0, 4),
    "beta": (13.0, 30.0, 18),
}
WINDOW_LIMITS: Mapping[str, tuple[str, float, float, int, tuple[str, ...]]] = {
    "response_onset_0_1": ("response_onset", 0.0, 1.0, 101, ("theta", "alpha", "beta")),
    "response_onset_0_0p5": ("response_onset", 0.0, 0.5, 51, ("theta", "alpha", "beta")),
    "response_end_0_1": ("response_end", 0.0, 1.0, 101, ("theta", "alpha", "beta")),
    "response_end_pre_0p5_0": ("response_end", -0.5, 0.0, 51, ("beta",)),
    "response_end_pmbr_0p5_1p5": ("response_end", 0.5, 1.5, 101, ("beta",)),
}
FIXED_ROI_CHANNELS: Mapping[str, tuple[str, ...]] = {
    "frontal_theta": ("F3", "FZ", "F4", "FCZ"),
    "posterior_alpha": ("P3", "PZ", "P4", "O1", "OZ", "O2"),
}
RIGHT_TASK_CONTRALATERAL = ("FC3", "C3", "CP3")
RIGHT_TASK_IPSILATERAL = ("FC4", "C4", "CP4")
LEFT_TASK_CONTRALATERAL = RIGHT_TASK_IPSILATERAL
LEFT_TASK_IPSILATERAL = RIGHT_TASK_CONTRALATERAL

TRIAL_REPEAT_COLUMNS: tuple[str, ...] = (
    "canonical_order_index",
    "canonical_event_key",
    "behavioural_row_key",
    "behavioural_id",
    "eeg_source_id",
    "source_recording_filename",
    "recording_stem",
    "source_shard",
    "tfr_row_index",
    "offset_session",
    "offset_block",
    "offset_trial",
    "task_filename",
    "source_file_role",
    "group",
    "performed_condition",
    "condition_semantics",
    "physical",
    "familiarity_repetition_condition",
    "complexity",
    "stimulus_mt",
    "avg_velocity",
    "mt",
    "mt_clip",
    "vresp",
    "accuracy_rating",
    "error",
    "vividness_rating",
    "original_handedness",
    "analysis_hand",
    "motor_laterality_mapping",
    "mapping_source",
    "mapping_override",
    "response_onset_seconds",
    "response_end_seconds",
    "response_duration_seconds",
    "response_duration_lt_0p5",
    "response_duration_lt_1p0",
    "scope_all_accepted",
    "scope_primary_blocks_1_5",
    "scope_imagery_final_overt_bridge",
    "strict_clean_eligibility",
    "duration_warning_flag",
    "continuous_qc_warning",
    "interpolation_count",
    "event_policy_version",
    "ledger_version",
    "continuous_v2_run_id",
    "continuous_v2_manifest_id",
    "continuous_source_sha256",
    "behavioural_inclusion_verified",
    "overt_tracing_filter_provenance",
    "imagery_cleanup_provenance",
    "artifact_red_on_any_flag",
    "artifact_source_any_flag",
    "artifact_source_robust_logp2p_flag",
    "artifact_source_jump_flag",
    "artifact_source_flat_flag",
    "feature_config_sha256",
    "source_tfr_run_manifest_sha256",
    "source_tfr_config_sha256",
    "source_tfr_recording_fingerprint",
)


@dataclass(frozen=True)
class FeatureDefinition:
    """One declared anchor/window/band channel feature."""

    name: str
    source_family: str
    minimum_seconds: float
    maximum_seconds: float
    expected_time_sample_count: int
    band: str
    minimum_hz: float
    maximum_hz: float
    expected_frequency_count: int
    role: str


def inclusive_indices(
    axis: np.ndarray,
    minimum: float,
    maximum: float,
    *,
    expected_count: int,
    label: str,
) -> np.ndarray:
    """Return exact inclusive indices and fail if endpoints/count differ."""

    values = np.asarray(axis, dtype=np.float64)
    indices = np.flatnonzero((values >= minimum) & (values <= maximum))
    if len(indices) != expected_count:
        raise RuntimeError(f"{label}_member_count_mismatch:{len(indices)}")
    selected = values[indices]
    if not (
        math.isclose(float(selected[0]), minimum, abs_tol=1e-12, rel_tol=0.0)
        and math.isclose(float(selected[-1]), maximum, abs_tol=1e-12, rel_tol=0.0)
    ):
        raise RuntimeError(f"{label}_endpoint_mismatch")
    if not np.array_equal(indices, np.arange(indices[0], indices[-1] + 1)):
        raise RuntimeError(f"{label}_members_not_contiguous")
    return indices


def feature_definitions(config: Mapping[str, Any]) -> list[FeatureDefinition]:
    """Expand the tracked window/band contract in deterministic order."""

    definitions: list[FeatureDefinition] = []
    for name, expected in WINDOW_LIMITS.items():
        family, tmin, tmax, time_count, bands = expected
        observed = config["feature_definitions"][name]
        if (
            observed["source_family"] != family
            or float(observed["minimum_seconds"]) != tmin
            or float(observed["maximum_seconds"]) != tmax
            or int(observed["expected_time_sample_count"]) != time_count
            or tuple(observed["bands"]) != bands
            or observed.get("endpoint_inclusive") is not True
        ):
            raise RuntimeError(f"feature_definition_differs:{name}")
        for band in bands:
            fmin, fmax, frequency_count = BAND_LIMITS[band]
            band_config = config["bands"][band]
            if (
                float(band_config["minimum_hz"]) != fmin
                or float(band_config["maximum_hz"]) != fmax
                or int(band_config["expected_frequency_count"]) != frequency_count
                or band_config.get("endpoint_inclusive") is not True
            ):
                raise RuntimeError(f"band_definition_differs:{band}")
            definitions.append(
                FeatureDefinition(
                    name=name,
                    source_family=family,
                    minimum_seconds=tmin,
                    maximum_seconds=tmax,
                    expected_time_sample_count=time_count,
                    band=band,
                    minimum_hz=fmin,
                    maximum_hz=fmax,
                    expected_frequency_count=frequency_count,
                    role=str(observed["role_by_band"][band]),
                )
            )
    if len(definitions) != EXPECTED_CHANNEL_DEFINITIONS:
        raise RuntimeError("channel_feature_definition_count_mismatch")
    return definitions


def load_and_validate_config(path: Path) -> tuple[dict[str, Any], str]:
    """Load the tracked feature contract and enforce accepted invariants."""

    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if config.get("schema_version") != 1 or config.get("pipeline_version") != FEATURE_NAMESPACE:
        raise RuntimeError("feature_config_identity_mismatch")
    accepted = config["accepted_surface"]
    expected_counts = {
        "trial_count": EXPECTED_TRIAL_COUNT,
        "participant_count": EXPECTED_PARTICIPANT_COUNT,
        "recording_count": EXPECTED_RECORDING_COUNT,
        "strict_clean_count": EXPECTED_STRICT_CLEAN_COUNT,
        "duration_warning_count": EXPECTED_DURATION_WARNING_COUNT,
        "file49_trial_count": EXPECTED_FILE49_COUNT,
        "participant_4_trial_count": EXPECTED_P4_TRIALS,
        "participant_60_trial_count": EXPECTED_P60_TRIALS,
        "response_duration_lt_0p5_count": EXPECTED_SHORT_0P5_COUNT,
        "response_duration_lt_1p0_count": EXPECTED_SHORT_1P0_COUNT,
        "channel_count": EXPECTED_CHANNEL_COUNT,
        "channel_feature_definitions_per_trial": EXPECTED_CHANNEL_DEFINITIONS,
        "roi_feature_definitions_per_trial": EXPECTED_ROI_DEFINITIONS,
        "channel_row_count": EXPECTED_CHANNEL_ROWS,
        "roi_row_count": EXPECTED_ROI_ROWS,
    }
    if any(int(accepted[key]) != value for key, value in expected_counts.items()):
        raise RuntimeError("accepted_feature_surface_count_mismatch")
    if any(bool(value) for key, value in config["safety"].items() if key != "process_one_recording_at_a_time"):
        raise RuntimeError("feature_config_enables_unauthorized_operation")
    if config["safety"].get("process_one_recording_at_a_time") is not True:
        raise RuntimeError("feature_config_must_process_one_recording_at_a_time")
    if config["representations"]["sensitivity"]["formula"] != "10 * log10(raw_power)":
        raise RuntimeError("log_power_formula_mismatch")
    if config["representations"]["aggregation"]["cell_transform_before_aggregation"] is not True:
        raise RuntimeError("aggregation_order_mismatch")
    if config["diagnostic_frequency_surface"]["feature_rows_created"] is not False:
        raise RuntimeError("31_40_hz_feature_creation_forbidden")
    if "low_beta" in config["bands"] or "high_beta" in config["bands"]:
        raise RuntimeError("beta_split_forbidden")
    definitions = feature_definitions(config)
    if any(item.maximum_hz > 30.0 for item in definitions):
        raise RuntimeError("feature_frequency_above_30_hz")
    validate_roi_config(config)
    validate_task_hand_config(config)
    return config, sha256_file(path)


def validate_roi_config(config: Mapping[str, Any]) -> None:
    """Require the exact fixed and task-hand-normalized ROI declarations."""

    rois = config["rois"]
    if tuple(rois["frontal_theta"]["physical_channels"]) != FIXED_ROI_CHANNELS["frontal_theta"]:
        raise RuntimeError("frontal_theta_roi_membership_mismatch")
    if tuple(rois["posterior_alpha"]["physical_channels"]) != FIXED_ROI_CHANNELS["posterior_alpha"]:
        raise RuntimeError("posterior_alpha_roi_membership_mismatch")
    contra = rois["contralateral_motor_beta"]
    ipsi = rois["ipsilateral_motor_beta"]
    if (
        tuple(contra["right_task_hand_channels"]) != RIGHT_TASK_CONTRALATERAL
        or tuple(contra["left_task_hand_channels"]) != LEFT_TASK_CONTRALATERAL
        or tuple(ipsi["right_task_hand_channels"]) != RIGHT_TASK_IPSILATERAL
        or tuple(ipsi["left_task_hand_channels"]) != LEFT_TASK_IPSILATERAL
    ):
        raise RuntimeError("motor_roi_membership_mismatch")
    expected_features = (
        "response_onset_0_1",
        "response_onset_0_0p5",
        "response_end_pre_0p5_0",
        "response_end_pmbr_0p5_1p5",
        "response_end_0_1",
    )
    if tuple(contra["features"]) != expected_features or tuple(ipsi["features"]) != expected_features:
        raise RuntimeError("motor_roi_feature_membership_mismatch")


def validate_task_hand_config(config: Mapping[str, Any]) -> None:
    """Require frozen r/l mappings and the explicit participant-70 override."""

    mappings = config["task_hand_mapping"]
    for code, hand in (("r", "right"), ("l", "left")):
        item = mappings["frozen_codes"][code]
        if (
            item["analysis_hand"] != hand
            or item["motor_laterality_mapping"] != f"{hand}_hand_task"
            or item["mapping_source"] != "frozen_handedness"
            or item["mapping_override"] is not False
        ):
            raise RuntimeError(f"frozen_task_hand_mapping_mismatch:{code}")
    override = mappings["owner_overrides"].get(70) or mappings["owner_overrides"].get("70")
    if override != {
        "required_original_handedness": "a",
        "analysis_hand": "right",
        "motor_laterality_mapping": "right_hand_task",
        "mapping_source": "owner_recollection",
        "mapping_override": True,
    }:
        raise RuntimeError("participant_70_task_hand_override_mismatch")


def apply_task_hand_mapping(trials: pd.DataFrame, config: Mapping[str, Any]) -> pd.DataFrame:
    """Add task-hand provenance while retaining the frozen handedness code."""

    result = trials.copy()
    result["original_handedness"] = result["handedness"].astype(str)
    result["analysis_hand"] = ""
    result["motor_laterality_mapping"] = ""
    result["mapping_source"] = ""
    result["mapping_override"] = False
    frozen = config["task_hand_mapping"]["frozen_codes"]
    for code, mapping in frozen.items():
        mask = result["original_handedness"].eq(str(code))
        for column in (
            "analysis_hand", "motor_laterality_mapping", "mapping_source", "mapping_override"
        ):
            result.loc[mask, column] = mapping[column]
    override_config = config["task_hand_mapping"]["owner_overrides"]
    for raw_id, mapping in override_config.items():
        participant = int(raw_id)
        mask = result["behavioural_id"].eq(participant)
        if not mask.any():
            raise RuntimeError(f"task_hand_override_participant_absent:{participant}")
        observed = set(result.loc[mask, "original_handedness"])
        required = str(mapping["required_original_handedness"])
        if observed != {required}:
            raise RuntimeError(f"task_hand_override_original_code_mismatch:{participant}:{observed}")
        for column in (
            "analysis_hand", "motor_laterality_mapping", "mapping_source", "mapping_override"
        ):
            result.loc[mask, column] = mapping[column]
    unresolved = result.loc[~result["analysis_hand"].isin(["right", "left"])]
    if len(unresolved):
        codes = unresolved[["behavioural_id", "original_handedness"]].drop_duplicates().to_dict("records")
        raise RuntimeError(f"unresolved_task_hand_mapping:{codes}")
    result["mapping_override"] = result["mapping_override"].astype(bool)
    return result


def standard_1005_channel_table(channels: Sequence[str] = PRIMARY_CHANNELS) -> pd.DataFrame:
    """Return deterministic standard_1005 Cartesian coordinates for 30 channels."""

    montage = mne.channels.make_standard_montage("standard_1005")
    positions = {name.upper(): value for name, value in montage.get_positions()["ch_pos"].items()}
    rows = []
    for index, channel in enumerate(channels):
        if channel not in positions:
            raise RuntimeError(f"standard_1005_coordinate_missing:{channel}")
        x, y, z = (float(value) for value in positions[channel])
        rows.append(
            {
                "channel_index": index,
                "physical_channel": channel,
                "montage": "standard_1005",
                "coordinate_frame": "head",
                "x_m": x,
                "y_m": y,
                "z_m": z,
            }
        )
    frame = pd.DataFrame(rows)
    if len(frame) != EXPECTED_CHANNEL_COUNT or not np.isfinite(frame[["x_m", "y_m", "z_m"]]).all().all():
        raise RuntimeError("channel_coordinate_contract_mismatch")
    return frame


def definition_table(config: Mapping[str, Any], frequencies: np.ndarray, times: np.ndarray) -> pd.DataFrame:
    """Return exact feature definitions with resolved source-axis membership."""

    rows = []
    for definition in feature_definitions(config):
        frequency_index = inclusive_indices(
            frequencies,
            definition.minimum_hz,
            definition.maximum_hz,
            expected_count=definition.expected_frequency_count,
            label=f"{definition.band}_frequency",
        )
        time_index = inclusive_indices(
            times,
            definition.minimum_seconds,
            definition.maximum_seconds,
            expected_count=definition.expected_time_sample_count,
            label=f"{definition.name}_time",
        )
        rows.append(
            {
                "feature_definition": definition.name,
                "source_family": definition.source_family,
                "window_minimum_seconds": definition.minimum_seconds,
                "window_maximum_seconds": definition.maximum_seconds,
                "window_endpoint_inclusive": True,
                "time_sample_count": len(time_index),
                "first_time_index": int(time_index[0]),
                "last_time_index": int(time_index[-1]),
                "band": definition.band,
                "frequency_minimum_hz": definition.minimum_hz,
                "frequency_maximum_hz": definition.maximum_hz,
                "frequency_endpoint_inclusive": True,
                "frequency_count": len(frequency_index),
                "first_frequency_index": int(frequency_index[0]),
                "last_frequency_index": int(frequency_index[-1]),
                "role_classification": definition.role,
                "aggregation_order": "cell_transform_then_equal_frequency_time_mean",
                "creates_31_40_hz_features": False,
            }
        )
    return pd.DataFrame(rows)


def roi_definition_table(config: Mapping[str, Any]) -> pd.DataFrame:
    """Return one row per predeclared ROI-feature combination."""

    rows = []
    for roi_name, roi in config["rois"].items():
        for feature_name in roi["features"]:
            if roi_name in FIXED_ROI_CHANNELS:
                channel_contract = ";".join(roi["physical_channels"])
            else:
                channel_contract = (
                    "right=" + ";".join(roi["right_task_hand_channels"])
                    + "|left=" + ";".join(roi["left_task_hand_channels"])
                )
            rows.append(
                {
                    "roi_name": roi_name,
                    "feature_definition": feature_name,
                    "band": roi["band"],
                    "normalized_laterality": roi["normalized_laterality"],
                    "channel_contract": channel_contract,
                    "weighting": "equal_weight_arithmetic_mean",
                    "role": roi["role"],
                }
            )
    frame = pd.DataFrame(rows)
    if len(frame) != EXPECTED_ROI_DEFINITIONS:
        raise RuntimeError("roi_definition_count_mismatch")
    return frame


def summarize_transformed_cells(
    db_power: np.ndarray,
    raw_power: np.ndarray,
    frequency_indices: np.ndarray,
    time_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Mean already-transformed dB and cell-level log raw power over cells."""

    frequency_slice = slice(int(frequency_indices[0]), int(frequency_indices[-1]) + 1)
    time_slice = slice(int(time_indices[0]), int(time_indices[-1]) + 1)
    db_cells = np.asarray(db_power[:, :, frequency_slice, time_slice], dtype=np.float64)
    raw_cells = np.asarray(raw_power[:, :, frequency_slice, time_slice], dtype=np.float64)
    if not np.isfinite(db_cells).all() or not np.isfinite(raw_cells).all() or np.any(raw_cells <= 0):
        raise RuntimeError("nonfinite_or_nonpositive_source_feature_cells")
    db_values = db_cells.mean(axis=(2, 3), dtype=np.float64)
    log_values = (10.0 * np.log10(raw_cells)).mean(axis=(2, 3), dtype=np.float64)
    if not np.isfinite(db_values).all() or not np.isfinite(log_values).all():
        raise RuntimeError("nonfinite_channel_feature_value")
    return db_values, log_values


def _repeat_trial_columns(trials: pd.DataFrame, channel_count: int) -> dict[str, np.ndarray]:
    """Repeat trial provenance once per physical channel."""

    missing = [column for column in TRIAL_REPEAT_COLUMNS if column not in trials.columns]
    if missing:
        raise RuntimeError(f"trial_feature_metadata_columns_missing:{missing}")
    return {
        column: np.repeat(trials[column].to_numpy(), channel_count)
        for column in TRIAL_REPEAT_COLUMNS
    }


def build_channel_rows(
    trials: pd.DataFrame,
    channels: pd.DataFrame,
    definition: FeatureDefinition,
    value_db: np.ndarray,
    value_log_power: np.ndarray,
    *,
    source_raw_sha256: str,
    source_db_sha256: str,
) -> pd.DataFrame:
    """Build one definition's trial-by-physical-channel long rows."""

    trial_count = len(trials)
    channel_count = len(channels)
    expected_shape = (trial_count, channel_count)
    if value_db.shape != expected_shape or value_log_power.shape != expected_shape:
        raise RuntimeError("channel_feature_value_shape_mismatch")
    data: dict[str, Any] = _repeat_trial_columns(trials, channel_count)
    for column in ("physical_channel", "channel_index", "montage", "coordinate_frame", "x_m", "y_m", "z_m"):
        data[column] = np.tile(channels[column].to_numpy(), trial_count)
    data.update(
        {
            "feature_definition": definition.name,
            "feature_anchor": definition.source_family,
            "window_minimum_seconds": definition.minimum_seconds,
            "window_maximum_seconds": definition.maximum_seconds,
            "window_endpoint_inclusive": True,
            "time_sample_count": definition.expected_time_sample_count,
            "band": definition.band,
            "frequency_minimum_hz": definition.minimum_hz,
            "frequency_maximum_hz": definition.maximum_hz,
            "frequency_endpoint_inclusive": True,
            "frequency_count": definition.expected_frequency_count,
            "role_classification": definition.role,
            "aggregation_order": "cell_transform_then_equal_frequency_time_mean",
            "value_db": value_db.reshape(-1),
            "value_log_power": value_log_power.reshape(-1),
            "source_tfr_raw_sha256": source_raw_sha256,
            "source_tfr_db_sha256": source_db_sha256,
        }
    )
    frame = pd.DataFrame(data)
    if len(frame) != trial_count * channel_count or frame[["value_db", "value_log_power"]].isna().any().any():
        raise RuntimeError("channel_feature_row_construction_failed")
    return frame


def _roi_from_selection(
    selection: pd.DataFrame,
    *,
    roi_name: str,
    normalized_laterality: str,
    role: str,
    channels: Sequence[str],
) -> pd.DataFrame:
    """Average a declared channel selection and preserve representative provenance."""

    key = "canonical_event_key"
    expected = len(channels)
    counts = selection.groupby(key, sort=False)["physical_channel"].nunique()
    if len(counts) == 0 or not counts.eq(expected).all():
        raise RuntimeError(f"roi_contributing_channel_count_mismatch:{roi_name}")
    if selection.duplicated([key, "physical_channel"]).any():
        raise RuntimeError(f"roi_duplicate_channel_rows:{roi_name}")
    exclude = {"physical_channel", "channel_index", "montage", "coordinate_frame", "x_m", "y_m", "z_m", "value_db", "value_log_power"}
    provenance_columns = [
        column for column in selection.columns if column not in exclude and column != key
    ]
    representative = selection.groupby(key, sort=False, as_index=False)[provenance_columns].first()
    values = selection.groupby(key, sort=False, as_index=False).agg(
        value_db=("value_db", "mean"), value_log_power=("value_log_power", "mean")
    )
    result = representative.merge(values, on=key, how="left", validate="one_to_one")
    result["roi_name"] = roi_name
    result["normalized_laterality"] = normalized_laterality
    result["physical_source_channels"] = ";".join(channels)
    result["contributing_channel_count"] = expected
    result["roi_weighting"] = "equal_weight_arithmetic_mean"
    result["roi_role"] = role
    return result


def derive_roi_rows(channel_rows: pd.DataFrame, config: Mapping[str, Any]) -> pd.DataFrame:
    """Derive all ROI rows only from the completed channel feature rows."""

    outputs: list[pd.DataFrame] = []
    for roi_name in ("frontal_theta", "posterior_alpha"):
        roi = config["rois"][roi_name]
        channels = tuple(roi["physical_channels"])
        for feature_name in roi["features"]:
            selected = channel_rows.loc[
                channel_rows["feature_definition"].eq(feature_name)
                & channel_rows["band"].eq(roi["band"])
                & channel_rows["physical_channel"].isin(channels)
            ]
            outputs.append(
                _roi_from_selection(
                    selected,
                    roi_name=roi_name,
                    normalized_laterality=str(roi["normalized_laterality"]),
                    role=str(roi["role"]),
                    channels=channels,
                )
            )
    for roi_name in ("contralateral_motor_beta", "ipsilateral_motor_beta"):
        roi = config["rois"][roi_name]
        for feature_name in roi["features"]:
            hand_parts = []
            for hand in ("right", "left"):
                channels = tuple(roi[f"{hand}_task_hand_channels"])
                selected = channel_rows.loc[
                    channel_rows["feature_definition"].eq(feature_name)
                    & channel_rows["band"].eq("beta")
                    & channel_rows["analysis_hand"].eq(hand)
                    & channel_rows["physical_channel"].isin(channels)
                ]
                if selected.empty:
                    continue
                hand_parts.append(
                    _roi_from_selection(
                        selected,
                        roi_name=roi_name,
                        normalized_laterality=str(roi["normalized_laterality"]),
                        role=str(roi["role"]),
                        channels=channels,
                    )
                )
            if not hand_parts:
                raise RuntimeError(f"motor_roi_has_no_task_hand_rows:{roi_name}:{feature_name}")
            outputs.append(pd.concat(hand_parts, ignore_index=True))
    result = pd.concat(outputs, ignore_index=True)
    result = result.sort_values(
        ["canonical_order_index", "roi_name", "feature_definition"], kind="stable"
    ).reset_index(drop=True)
    expected_rows = channel_rows["canonical_event_key"].nunique() * EXPECTED_ROI_DEFINITIONS
    if len(result) != expected_rows:
        raise RuntimeError(f"roi_feature_row_count_mismatch:{len(result)}:{expected_rows}")
    if result.duplicated(["canonical_event_key", "roi_name", "feature_definition", "band"]).any():
        raise RuntimeError("duplicate_roi_feature_identity")
    if not np.isfinite(result[["value_db", "value_log_power"]]).all().all():
        raise RuntimeError("nonfinite_roi_feature_value")
    return result


def roi_equality_validation(channel_rows: pd.DataFrame, roi_rows: pd.DataFrame) -> dict[str, float | int]:
    """Recalculate every ROI from channel rows and report maximum differences."""

    identity = ["canonical_event_key", "roi_name", "feature_definition", "band"]
    requested = roi_rows[
        [
            *identity,
            "physical_source_channels",
            "contributing_channel_count",
            "value_db",
            "value_log_power",
        ]
    ].copy()
    requested["physical_channel"] = requested["physical_source_channels"].str.split(";")
    requested = requested.explode("physical_channel", ignore_index=True)
    source = channel_rows[
        [
            "canonical_event_key",
            "feature_definition",
            "band",
            "physical_channel",
            "value_db",
            "value_log_power",
        ]
    ].rename(columns={"value_db": "source_value_db", "value_log_power": "source_value_log_power"})
    linked = requested.merge(
        source,
        on=["canonical_event_key", "feature_definition", "band", "physical_channel"],
        how="left",
        validate="many_to_one",
    )
    recalculated = linked.groupby(identity, sort=False, as_index=False).agg(
        observed_channel_count=("source_value_db", "count"),
        recalculated_db=("source_value_db", "mean"),
        recalculated_log_power=("source_value_log_power", "mean"),
        expected_channel_count=("contributing_channel_count", "first"),
        persisted_db=("value_db", "first"),
        persisted_log_power=("value_log_power", "first"),
    )
    if not recalculated["observed_channel_count"].eq(recalculated["expected_channel_count"]).all():
        raise RuntimeError("roi_recalculation_channel_count_mismatch")
    differences = np.column_stack(
        [
            np.abs(recalculated["persisted_db"] - recalculated["recalculated_db"]),
            np.abs(
                recalculated["persisted_log_power"]
                - recalculated["recalculated_log_power"]
            ),
        ]
    )
    return {
        "rows_checked": len(roi_rows),
        "maximum_absolute_db_difference": float(differences[:, 0].max(initial=0.0)),
        "maximum_absolute_log_power_difference": float(differences[:, 1].max(initial=0.0)),
    }


def require_ignored_output_root(repo_root: Path, output_root: Path) -> None:
    """Require the output beneath this repository and ignored by Git."""

    resolved_repo = repo_root.resolve()
    resolved_output = output_root.resolve()
    if resolved_output == resolved_repo or resolved_repo not in resolved_output.parents:
        raise RuntimeError("unsafe_feature_output_root")
    relative = resolved_output.relative_to(resolved_repo).as_posix()
    result = subprocess.run(
        ["git", "check-ignore", "--quiet", relative], cwd=resolved_repo, check=False
    )
    if result.returncode != 0:
        raise RuntimeError("feature_output_root_not_ignored")


def atomic_write_parquet(path: Path, frame: pd.DataFrame, *, compression: str = "zstd") -> None:
    """Write one Parquet file beside its destination and atomically replace."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        frame.to_parquet(temporary, index=False, compression=compression)
        reopened = pd.read_parquet(temporary)
        if len(reopened) != len(frame) or list(reopened.columns) != list(frame.columns):
            raise RuntimeError(f"parquet_reopen_contract_mismatch:{path}")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Write deterministic JSON atomically."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        json.loads(temporary.read_text(encoding="utf-8"))
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


def atomic_combine_parquet(shards: Sequence[Path], destination: Path) -> dict[str, Any]:
    """Stream compatible Parquet shards into one atomically published file."""

    if not shards:
        raise RuntimeError("no_parquet_shards_to_combine")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.parent / f".{destination.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    writer: pq.ParquetWriter | None = None
    row_count = 0
    try:
        for shard in shards:
            table = pq.read_table(shard)
            if writer is None:
                writer = pq.ParquetWriter(temporary, table.schema, compression="zstd")
            elif table.schema != writer.schema:
                raise RuntimeError(f"aggregate_parquet_schema_mismatch:{shard}")
            writer.write_table(table)
            row_count += table.num_rows
        assert writer is not None
        writer.close()
        writer = None
        metadata = pq.read_metadata(temporary)
        if metadata.num_rows != row_count:
            raise RuntimeError("aggregate_parquet_row_count_mismatch")
        os.replace(temporary, destination)
    finally:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)
    return {
        "row_count": row_count,
        "size_bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "format": "parquet",
    }


def file_descriptor(path: Path, *, root: Path, row_count: int | None = None) -> dict[str, Any]:
    """Return a relative-path, size, hash, and optional row-count descriptor."""

    descriptor: dict[str, Any] = {
        "relative_path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if row_count is not None:
        descriptor["row_count"] = int(row_count)
    return descriptor


def source_snapshot(paths: Iterable[Path], *, root: Path) -> list[dict[str, Any]]:
    """Capture deterministic size/mtime evidence for immutable source files."""

    return [
        {
            "relative_path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
        }
        for path in sorted((Path(item) for item in paths), key=lambda item: item.as_posix())
    ]


def compare_source_snapshots(before: Sequence[Mapping[str, Any]], after: Sequence[Mapping[str, Any]]) -> None:
    """Fail closed if a TFR source path, size, or mtime changes."""

    if list(before) != list(after):
        raise RuntimeError("source_tfr_size_or_mtime_changed")


def frame_descriptor(path: Path, frame: pd.DataFrame, *, root: Path) -> dict[str, Any]:
    """Describe a persisted frame including stable logical content hash."""

    descriptor = file_descriptor(path, root=root, row_count=len(frame))
    descriptor["content_sha256"] = stable_frame_hash(frame)
    descriptor["format"] = "parquet"
    return descriptor
