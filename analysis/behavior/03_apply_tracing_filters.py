"""Apply the historical DEMI TraceLab tracing filters.

This script reads the reconstructed TraceLab tables created by
``analysis/behavior/02_reconstruct_tracing_tables.py`` and applies the
sample-level and trial-level tracing cleanup rules from the historical R
pipeline. The work here is behavioural/tracing infrastructure for the EEG
reanalysis: later event-alignment code needs a deterministic record of where a
physical trace really began, how long the filtered trace lasted, and which
trials were affected by the old TraceLab cleanup rules before any comparison
with EEG triggers is attempted.

The old authority source is ``_Scripts/01_preprocessing.R``, with thresholds
from ``_Scripts/_settings.R`` and predicate helpers from
``_Scripts/_functions/filters.R`` and ``_Scripts/_functions/complexity.R``.
This Python script deliberately ports that solved tracing-filter logic rather
than reinterpreting the raw TraceLab files. The old R column ``id`` maps to
``participant_id`` in the reconstructed Python tables; the script keeps
``participant_id`` in outputs and comments whenever an old-R-to-Python mapping
matters.

Old logic ported now:

- origin extraction from the first TraceLab figure point for each trial;
- re-zeroing tracing time within trial before grouped lag/lead operations;
- distance from prior tracing sample and distance from origin;
- repeated-point removal;
- no-shape trial flagging from figure-to-tracing size ratio;
- failed trial-end trimming with both old ``trial_done`` parameter sets;
- touchscreen glitch point flagging, with the hard-coded ``(239, 1079)``
  coordinate case preserved;
- false-start sample flagging and ``start_shift`` recovery;
- hand-noise sample flagging;
- incomplete-tracing trial flagging from end gap, size ratio, length ratio, and
  lateral-shift difference;
- excessive time/distance gap trial flagging;
- edge-hit trial flagging without dropping those trials;
- final filtered movement time and response path-length summaries needed by
  later event-offset reconstruction.

Work intentionally deferred:

- task-row joins and task-to-tracing event-offset construction;
- tracing/figure interpolation and accuracy/error metrics;
- imagery-trial cleanup from the old behavioural analysis path;
- raw EEG annotation inspection, EEG preprocessing, and EEG quality-control
  decisions;
- participant-level decisions of any kind;
- old plotting side effects from ``plot_filters <- TRUE``.

Expected inputs, all local-only:

- ``_Data/behavior/tracing_tables/points.parquet``;
- ``_Data/behavior/tracing_tables/segments.parquet``;
- ``_Data/behavior/tracing_tables/frames.parquet``;
- ``_Data/behavior/tracing_tables/tracings.parquet``;
- ``_Data/behavior/tracing_tables/zip_parse_manifest.csv``.

Generated local outputs:

- ``_Data/behavior/tracing_filters/tracings_with_filter_flags.parquet``;
- ``_Data/behavior/tracing_filters/filtered_tracings.parquet``;
- ``_Data/behavior/tracing_filters/tracing_trial_filter_summary.parquet``;
- ``_Data/behavior/tracing_filters/trace_filter_info.parquet``;
- one parquet audit table per filter stage;
- ``_Data/behavior/tracing_filters/tracing_filter_summary.md``.

Safety boundaries:

- all source tables are read-only inputs;
- outputs are written only below ``_Data/behavior/tracing_filters/``;
- the script does not inspect raw EEG, preprocess EEG, or join to EEG events;
- the script does not revise the published behavioural analysis;
- filter outputs are audit records, not participant-level decisions.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/03_apply_tracing_filters.py

Paths are resolved relative to this file, so the command also works from
another current working directory.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


TRIAL_KEY = ["participant_id", "session", "block", "trial"]

TRACING_TABLE_DIR = Path("_Data") / "behavior" / "tracing_tables"
FILTER_OUTPUT_DIR = Path("_Data") / "behavior" / "tracing_filters"

INPUT_FILES = {
    "points": "points.parquet",
    "segments": "segments.parquet",
    "frames": "frames.parquet",
    "tracings": "tracings.parquet",
    "zip_manifest": "zip_parse_manifest.csv",
}

REQUIRED_COLUMNS = {
    "points": {*TRIAL_KEY, "point_index", "x", "y"},
    "segments": {*TRIAL_KEY, "segment_index", "start.x", "start.y", "end.x", "end.y", "ctrl.x", "ctrl.y"},
    "frames": {*TRIAL_KEY, "sample_index", "x", "y", "time"},
    "tracings": {*TRIAL_KEY, "sample_index", "x", "y", "time"},
    "zip_manifest": {
        *TRIAL_KEY,
        "source_zip_basename",
        "issue_codes",
        "frames_parse_status",
        "frames_row_count",
        "tracings_parse_status",
        "tracings_row_count",
    },
}

SCREEN_RES_X = 1920
SCREEN_RES_Y = 1080

NO_SHAPE_PARAMS = {
    "min_size_ratio": 5.0,
}

DONE_FILTER_PARAMS = {
    "origin_radius": 50.0,
    "end_radius": 400.0,
    "pause_radius": 200.0,
    "min_prop": 0.62,
    "end_prop": 0.82,
    "min_pause": 0.100,
}

DONE_FILTER_PARAMS_2 = {
    "origin_radius": 50.0,
    "end_radius": 400.0,
    "pause_radius": 200.0,
    "min_prop": 0.62,
    "end_prop": 0.55,
    "min_pause": 0.400,
}

GLITCH_FILTER_PARAMS = {
    "min_angle_diff": 88.0,
    "min_angle_diff_alt": 70.0,
    "min_dist": 173.0,
}

FALSE_START_PARAMS = {
    "start_radius": 80.0,
    "max_prop": 0.20,
    "min_pause": 0.064,
}

HAND_NOISE_PARAMS = {
    "min_sharp_turnsum": 90.0,
    "max_size": 100.0,
    "max_timediff": 0.100,
    "min_dist_a": 120.0,
    "min_dist_b": 300.0,
    "min_angle_diff": 30.0,
    "min_end_dist": 100.0,
}

INCOMPLETE_PARAMS = {
    "min_end_gap": 300.0,
    "min_size_ratio": 2.5,
    "min_length_ratio": 1.45,
    "min_shift_diff": 0.71,
}

GAP_FILTER_PARAMS = {
    "min_pause": 0.170,
    "min_timegap": 0.070,
    "min_dist_b": 500.0,
    "min_dist_c": 250.0,
    "min_dist_d": 180.0,
    "min_turnsum_c": 90.0,
    "min_turnsum_d": 120.0,
}

# This watchlist is copied from prior local inventory reviews at Tony's
# request. It is used only in the local summary to surface rows needing human
# review; it is not a decision rule.
REVIEW_WATCHLIST_IDS = {
    "incomplete_task_marker": [2, 12, 20, 22, 27, 39, 88, 89, 99],
    "no_local_figure_zip_in_prior_manifest": [9, 10, 12, 13, 14, 20, 24, 26, 36, 78, 96],
    "nonstandard_trace_table_coverage": [2, 22, 27, 39, 88, 89, 99, 100, 999],
}

STAGE_OUTPUT_FILES = {
    "repeated_point": "repeated_points.parquet",
    "no_shape": "no_shape_trials.parquet",
    "failed_end": "failed_end_trials.parquet",
    "glitch": "glitch_trials.parquet",
    "false_start": "false_start_trials.parquet",
    "hand_noise": "hand_noise_trials.parquet",
    "incomplete": "incomplete_trials.parquet",
    "large_gap": "large_gap_trials.parquet",
    "hit_edge": "edge_hit_trials.parquet",
    "no_tracing_payload": "no_tracing_payload_trials.parquet",
}

STAGE_ORDER = {
    "repeated_point": 1,
    "no_shape": 2,
    "failed_end": 3,
    "glitch": 4,
    "false_start": 5,
    "hand_noise": 6,
    "incomplete": 7,
    "large_gap": 8,
    "hit_edge": 9,
    "no_tracing_payload": 10,
}

FLAG_COLUMNS = [
    "repeated_point",
    "no_shape_trial",
    "failed_end_done",
    "glitch_xy",
    "glitch",
    "false_start",
    "hnoise",
    "incomplete_trial",
    "gap",
    "large_gap_trial",
    "hit_edge_trial",
    "kept_after_filters",
]


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Args:
        None.

    Returns:
        Absolute repository root path.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def require_pyarrow() -> None:
    """Verify that pyarrow-backed parquet I/O is available.

    Args:
        None.

    Returns:
        None.

    Side effects:
        Imports ``pyarrow``. Raises ``RuntimeError`` with setup context if the
        active environment cannot read and write parquet files.
    """

    try:
        import pyarrow  # noqa: F401
    except ImportError as error:
        raise RuntimeError(
            "pyarrow is required for tracing-filter parquet inputs and outputs. "
            "Run this script in the project environment."
        ) from error


def relative_to_repo(path: Path, repo_root: Path) -> str:
    """Render a path relative to the repository when possible.

    Args:
        path: Path to render.
        repo_root: Repository root used for relative paths.

    Returns:
        POSIX path text.

    Side effects:
        None.
    """

    try:
        return path.relative_to(repo_root).as_posix()
    except ValueError:
        return path.as_posix()


def validate_columns(table_name: str, frame: pd.DataFrame) -> None:
    """Check that a loaded input table has the columns required by this port.

    Args:
        table_name: Logical table name such as ``tracings``.
        frame: Loaded input table.

    Returns:
        None.

    Side effects:
        Raises ``RuntimeError`` when a required column is absent.
    """

    required = REQUIRED_COLUMNS[table_name]
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise RuntimeError(f"{table_name} input is missing required column(s): {', '.join(missing)}")


def read_input_tables(repo_root: Path) -> dict[str, pd.DataFrame]:
    """Read all reconstructed TraceLab tables needed by the filter port.

    Args:
        repo_root: Repository root used to resolve local input paths.

    Returns:
        Dictionary of input ``DataFrame`` objects keyed by logical table name.

    Side effects:
        Reads parquet and CSV files below ``_Data/behavior/tracing_tables/``.
        Raises ``RuntimeError`` with a clear path if any input is missing or
        malformed.
    """

    require_pyarrow()
    input_dir = repo_root / TRACING_TABLE_DIR
    tables: dict[str, pd.DataFrame] = {}

    for table_name, filename in INPUT_FILES.items():
        path = input_dir / filename
        if not path.exists():
            raise RuntimeError(
                f"required input is missing: {relative_to_repo(path, repo_root)}. "
                "Run analysis/behavior/02_reconstruct_tracing_tables.py first."
            )

        try:
            if path.suffix == ".parquet":
                frame = pd.read_parquet(path)
            else:
                frame = pd.read_csv(path)
        except Exception as error:  # noqa: BLE001 - add path context.
            raise RuntimeError(f"failed to read {relative_to_repo(path, repo_root)}: {error}") from error

        validate_columns(table_name, frame)
        tables[table_name] = coerce_trial_key_columns(frame)

    return tables


def coerce_trial_key_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Coerce participant/session/block/trial columns to nullable integers.

    Args:
        frame: Input table containing the trial key columns.

    Returns:
        Copy of ``frame`` with normalized key dtypes.

    Side effects:
        None.
    """

    out = frame.copy()
    for column in TRIAL_KEY:
        out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
    return out


def sort_by_trial_and_order(frame: pd.DataFrame, order_column: str) -> pd.DataFrame:
    """Sort rows by trial key and explicit within-trial order.

    Args:
        frame: Table to sort.
        order_column: Column carrying TraceLab row order, such as
            ``sample_index`` or ``point_index``.

    Returns:
        Stable sorted copy.

    Side effects:
        None.
    """

    sort_columns = [*TRIAL_KEY, order_column]
    return frame.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)


def line_length(start_x: Any, start_y: Any, end_x: Any, end_y: Any) -> Any:
    """Compute Euclidean distance between coordinate pairs.

    Args:
        start_x: Starting x coordinate or vector.
        start_y: Starting y coordinate or vector.
        end_x: Ending x coordinate or vector.
        end_y: Ending y coordinate or vector.

    Returns:
        Scalar or vector of Euclidean distances.

    Side effects:
        None.
    """

    return np.sqrt((end_x - start_x) ** 2 + (end_y - start_y) ** 2)


def safe_divide(numerator: Any, denominator: Any) -> Any:
    """Divide while preserving numpy's ``inf``/``nan`` behaviour.

    Args:
        numerator: Scalar or vector numerator.
        denominator: Scalar or vector denominator.

    Returns:
        Division result with divide-by-zero warnings suppressed.

    Side effects:
        None.
    """

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.asarray(numerator, dtype="float64") / np.asarray(denominator, dtype="float64")


def angle_diffs_degrees(frame: pd.DataFrame, skip: int = 0) -> pd.Series:
    """Port old ``get_angle_diffs`` and return degrees.

    Args:
        frame: Sorted tracing table with ``x`` and ``y`` columns.
        skip: Number of following segments to combine with the current segment,
            matching the old R helper's ``skip`` argument.

    Returns:
        ``Series`` of signed angle differences in degrees.

    Side effects:
        None.
    """

    grouped = frame.groupby(TRIAL_KEY, sort=False, dropna=False)
    dx = frame["x"] - grouped["x"].shift(1)
    dy = frame["y"] - grouped["y"].shift(1)

    if skip > 0:
        dx2 = dx + grouped_apply_shift(dx, frame, -skip)
        dy2 = dy + grouped_apply_shift(dy, frame, -skip)
    else:
        dx2 = dx
        dy2 = dy

    prev_dx = grouped_apply_shift(dx, frame, 1)
    prev_dy = grouped_apply_shift(dy, frame, 1)
    theta = np.arctan2(dy2, dx2) - np.arctan2(prev_dy, prev_dx)
    theta = np.where(np.abs(theta) > np.pi, theta - (2 * np.pi) * np.sign(theta), theta)
    return pd.Series((theta / np.pi) * 180.0, index=frame.index)


def grouped_apply_shift(values: pd.Series, frame: pd.DataFrame, periods: int) -> pd.Series:
    """Shift an already-computed vector within trial groups.

    Args:
        values: Series aligned to ``frame``.
        frame: DataFrame containing trial key columns for grouping.
        periods: Shift amount; negative values are old R ``lead`` operations.

    Returns:
        Grouped shift result aligned to ``frame``.

    Side effects:
        None.
    """

    temp = pd.DataFrame({**{column: frame[column] for column in TRIAL_KEY}, "value": values})
    return temp.groupby(TRIAL_KEY, sort=False, dropna=False)["value"].shift(periods)


def build_origins(points: pd.DataFrame) -> pd.DataFrame:
    """Extract the first TraceLab figure point as the trial origin.

    Args:
        points: Reconstructed ``points`` table.

    Returns:
        One row per trial with ``origin.x`` and ``origin.y``.

    Side effects:
        None.
    """

    sorted_points = sort_by_trial_and_order(points, "point_index")
    origins = sorted_points.groupby(TRIAL_KEY, sort=False, dropna=False).first().reset_index()
    return origins[[*TRIAL_KEY, "x", "y"]].rename(columns={"x": "origin.x", "y": "origin.y"})


def add_basic_trace_metrics(frame: pd.DataFrame, zero_time: bool) -> pd.DataFrame:
    """Add grouped time, distance, and origin-distance metrics.

    Args:
        frame: Sorted tracing samples with origin columns.
        zero_time: Whether to subtract the first within-trial time value.

    Returns:
        Copy with ``time``, ``dist``, ``origin.dist``, and ``timediff`` updated.

    Side effects:
        None.
    """

    out = sort_by_trial_and_order(frame, "sample_index").copy()
    grouped = out.groupby(TRIAL_KEY, sort=False, dropna=False)

    if zero_time:
        out["time"] = out["raw_time"] - grouped["raw_time"].transform("first")

    lag_x = grouped["x"].shift(1)
    lag_y = grouped["y"].shift(1)
    lag_time = grouped["time"].shift(1)
    out["dist"] = line_length(lag_x, lag_y, out["x"], out["y"])
    out["origin.dist"] = line_length(out["origin.x"], out["origin.y"], out["x"], out["y"])
    out["timediff"] = (out["time"] - lag_time).fillna(0.0)
    return out


def prepare_response_samples(tracings: pd.DataFrame, origins: pd.DataFrame) -> pd.DataFrame:
    """Create the old ``responsedat`` starting table.

    Args:
        tracings: Reconstructed physical tracing samples.
        origins: Trial origin table from ``build_origins``.

    Returns:
        Tracing samples with stable row IDs, old-R-style zeroed time, origin
        coordinates, distance from prior sample, and distance from origin.

    Side effects:
        Raises ``RuntimeError`` if any tracing sample lacks a matching origin.
    """

    sorted_tracings = sort_by_trial_and_order(tracings, "sample_index").copy()
    sorted_tracings.insert(0, "tracing_row_id", np.arange(1, len(sorted_tracings) + 1, dtype="int64"))
    sorted_tracings["raw_time"] = sorted_tracings["time"].astype("float64")

    merged = sorted_tracings.merge(origins, on=TRIAL_KEY, how="left", validate="many_to_one")
    missing_origin = merged["origin.x"].isna() | merged["origin.y"].isna()
    if missing_origin.any():
        examples = merged.loc[missing_origin, TRIAL_KEY].drop_duplicates().head(10).to_dict("records")
        raise RuntimeError(f"{int(missing_origin.sum())} tracing rows lack origin points; first keys: {examples}")

    return add_basic_trace_metrics(merged, zero_time=True)


def initialize_flag_table(responses: pd.DataFrame) -> pd.DataFrame:
    """Create a full-row filter flag table before stage filtering begins.

    Args:
        responses: Initial tracing response table from ``prepare_response_samples``.

    Returns:
        One row per reconstructed tracing sample, indexed by ``tracing_row_id``.

    Side effects:
        None.
    """

    keep_columns = [
        "tracing_row_id",
        *TRIAL_KEY,
        "sample_index",
        "x",
        "y",
        "raw_time",
        "time",
        "origin.x",
        "origin.y",
        "dist",
        "origin.dist",
    ]
    flags = responses[keep_columns].rename(
        columns={
            "time": "time_from_original_trace_start",
            "dist": "initial_dist",
            "origin.dist": "initial_origin.dist",
        }
    )
    for column in FLAG_COLUMNS:
        flags[column] = False
    flags["filter_stage_removed"] = ""
    flags["filter_stage_removed_order"] = pd.NA
    return flags.set_index("tracing_row_id", drop=False)


def set_sample_flag(flags: pd.DataFrame, rows: pd.DataFrame, flag_column: str) -> None:
    """Set a sample-level flag for rows selected during one stage.

    Args:
        flags: Full tracing flag table indexed by ``tracing_row_id``.
        rows: Tracing rows whose row IDs should be flagged.
        flag_column: Boolean flag column to set.

    Returns:
        None.

    Side effects:
        Mutates ``flags`` in place.
    """

    if rows.empty:
        return
    flags.loc[rows["tracing_row_id"].to_numpy(), flag_column] = True


def set_trial_flag(flags: pd.DataFrame, trials: pd.DataFrame, flag_column: str) -> None:
    """Set a trial-level flag for all original samples from selected trials.

    Args:
        flags: Full tracing flag table indexed by ``tracing_row_id``.
        trials: Trial-level table containing the trial key columns.
        flag_column: Boolean flag column to set.

    Returns:
        None.

    Side effects:
        Mutates ``flags`` in place.
    """

    if trials.empty:
        return
    mask = trial_membership_mask(flags, trials)
    flags.loc[mask, flag_column] = True


def set_removed_stage(flags: pd.DataFrame, rows: pd.DataFrame, stage_name: str) -> None:
    """Record the first stage that removed each affected sample.

    Args:
        flags: Full tracing flag table indexed by ``tracing_row_id``.
        rows: Tracing rows removed by ``stage_name``.
        stage_name: Stable stage name.

    Returns:
        None.

    Side effects:
        Mutates ``flags`` in place.
    """

    if rows.empty:
        return
    row_ids = rows["tracing_row_id"].to_numpy()
    unset = flags.loc[row_ids, "filter_stage_removed"].eq("")
    unset_ids = flags.loc[row_ids].loc[unset, "tracing_row_id"].to_numpy()
    if len(unset_ids) == 0:
        return
    flags.loc[unset_ids, "filter_stage_removed"] = stage_name
    flags.loc[unset_ids, "filter_stage_removed_order"] = STAGE_ORDER[stage_name]


def trial_membership_mask(rows: pd.DataFrame, trials: pd.DataFrame) -> pd.Series:
    """Return a mask for rows whose trial key appears in ``trials``.

    Args:
        rows: Sample-level or trial-level DataFrame with trial key columns.
        trials: Trial-level DataFrame with trial key columns.

    Returns:
        Boolean Series aligned to ``rows``.

    Side effects:
        None.
    """

    if rows.empty or trials.empty:
        return pd.Series(False, index=rows.index)
    row_keys = pd.MultiIndex.from_frame(rows[TRIAL_KEY])
    trial_keys = pd.MultiIndex.from_frame(trials[TRIAL_KEY].drop_duplicates())
    return pd.Series(row_keys.isin(trial_keys), index=rows.index)


def build_frame_info(frames: pd.DataFrame) -> pd.DataFrame:
    """Summarize figure frame geometry for filter comparisons.

    Args:
        frames: Reconstructed animation-frame table.

    Returns:
        One row per trial with figure size, lateral shift, and path length.

    Side effects:
        None.
    """

    sorted_frames = sort_by_trial_and_order(frames, "sample_index").copy()
    grouped = sorted_frames.groupby(TRIAL_KEY, sort=False, dropna=False)
    sorted_frames["seglen"] = line_length(
        grouped["x"].shift(1),
        grouped["y"].shift(1),
        sorted_frames["x"],
        sorted_frames["y"],
    )

    summary = grouped.agg(
        frame_rows=("sample_index", "size"),
        fig_x_min=("x", "min"),
        fig_x_max=("x", "max"),
        fig_y_min=("y", "min"),
        fig_y_max=("y", "max"),
        fig_len=("seglen", "sum"),
    ).reset_index()
    summary["fig_width"] = summary["fig_x_max"] - summary["fig_x_min"]
    summary["fig_height"] = summary["fig_y_max"] - summary["fig_y_min"]
    summary["fig_max_size"] = summary[["fig_width", "fig_height"]].max(axis=1)
    # Old R measured lateral shift against the screen midpoint, not the
    # figure's own origin. Preserve that convention exactly.
    summary["fig_lat_shift"] = (safe_divide(summary["fig_x_min"] - (SCREEN_RES_X / 2), summary["fig_width"]) + 0.5) * 2
    return summary[
        [
            *TRIAL_KEY,
            "frame_rows",
            "fig_x_min",
            "fig_x_max",
            "fig_y_min",
            "fig_y_max",
            "fig_width",
            "fig_height",
            "fig_max_size",
            "fig_lat_shift",
            "fig_len",
        ]
    ]


def build_raw_tracing_info(tracings: pd.DataFrame) -> pd.DataFrame:
    """Summarize raw reconstructed tracing samples before any filter stage.

    Args:
        tracings: Reconstructed physical tracing samples.

    Returns:
        One row per trial with raw sample counts and raw TraceLab time range.

    Side effects:
        None.
    """

    if tracings.empty:
        return pd.DataFrame(columns=[*TRIAL_KEY, "input_tracing_rows", "raw_trace_onset", "raw_trace_end"])

    sorted_tracings = sort_by_trial_and_order(tracings, "sample_index")
    summary = sorted_tracings.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
        input_tracing_rows=("sample_index", "size"),
        raw_trace_onset=("time", "min"),
        raw_trace_end=("time", "max"),
    ).reset_index()
    summary["raw_trace_duration"] = summary["raw_trace_end"] - summary["raw_trace_onset"]
    return summary


def apply_repeated_point_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop repeated points exactly where old R used ``dist > 0``.

    Args:
        responses: Current tracing samples with ``dist`` from prior sample.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)`` where audit rows are repeated samples.

    Side effects:
        Mutates ``flags`` by marking repeated samples and their removal stage.
    """

    repeated_mask = responses["dist"].notna() & (responses["dist"] <= 0)
    repeated_rows = responses.loc[repeated_mask].copy()
    set_sample_flag(flags, repeated_rows, "repeated_point")
    set_removed_stage(flags, repeated_rows, "repeated_point")

    audit = repeated_rows[
        [
            "tracing_row_id",
            *TRIAL_KEY,
            "sample_index",
            "x",
            "y",
            "raw_time",
            "time",
            "dist",
            "origin.dist",
        ]
    ].copy()
    audit["stage"] = "repeated_point"
    audit["stage_order"] = STAGE_ORDER["repeated_point"]

    kept = responses.loc[~repeated_mask].copy()
    return kept, audit


def summarize_trace_shape(responses: pd.DataFrame) -> pd.DataFrame:
    """Summarize current tracing size and duration for trial-level filters.

    Args:
        responses: Current tracing samples.

    Returns:
        One row per trial with sample count, duration, width, height, and max
        dimension.

    Side effects:
        None.
    """

    if responses.empty:
        return pd.DataFrame(columns=[*TRIAL_KEY, "samples", "duration", "trace_width", "trace_height", "max_size"])

    summary = responses.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
        samples=("tracing_row_id", "size"),
        duration=("time", "max"),
        trace_x_min=("x", "min"),
        trace_x_max=("x", "max"),
        trace_y_min=("y", "min"),
        trace_y_max=("y", "max"),
    ).reset_index()
    summary["trace_width"] = summary["trace_x_max"] - summary["trace_x_min"]
    summary["trace_height"] = summary["trace_y_max"] - summary["trace_y_min"]
    summary["max_size"] = summary[["trace_width", "trace_height"]].max(axis=1)
    return summary


def apply_no_shape_stage(
    responses: pd.DataFrame,
    frame_info: pd.DataFrame,
    flags: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag and drop trials whose tracing has no formed shape.

    Args:
        responses: Current tracing samples after repeated-point removal.
        frame_info: Trial-level figure geometry summaries.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)`` where audit rows are no-shape trials.

    Side effects:
        Mutates ``flags`` with the ``no_shape_trial`` flag and removal stage.
    """

    trace_shape = summarize_trace_shape(responses)
    audit = trace_shape.merge(frame_info[[*TRIAL_KEY, "fig_max_size"]], on=TRIAL_KEY, how="left")
    audit["size_ratio"] = safe_divide(audit["fig_max_size"], audit["max_size"])
    audit["no_shape"] = audit["size_ratio"] > NO_SHAPE_PARAMS["min_size_ratio"]
    audit = audit.loc[audit["no_shape"]].copy()
    audit["stage"] = "no_shape"
    audit["stage_order"] = STAGE_ORDER["no_shape"]

    affected_mask = trial_membership_mask(responses, audit)
    affected_rows = responses.loc[affected_mask].copy()
    set_trial_flag(flags, audit, "no_shape_trial")
    set_removed_stage(flags, affected_rows, "no_shape")
    kept = responses.loc[~affected_mask].copy()
    return kept, audit


def trial_done_for_group(group: pd.DataFrame, params: dict[str, float]) -> pd.Series:
    """Port the old ``trial_done`` predicate for one trial.

    Args:
        group: One trial's current tracing samples.
        params: Done-filter threshold dictionary from ``_settings.R``.

    Returns:
        Boolean Series aligned to ``group`` marking points after a failed
        trial-end attempt.

    Side effects:
        None.
    """

    origin_dist = group["origin.dist"].astype("float64")
    timediff = group["timediff"].astype("float64")
    n_samples = len(group)
    if n_samples == 0:
        return pd.Series(False, index=group.index)

    prop = pd.Series(np.arange(1, n_samples + 1, dtype="float64") / n_samples, index=group.index)
    lag_origin = origin_dist.shift(1)
    lag_on_origin = (lag_origin < params["origin_radius"]).fillna(False)
    not_first = lag_origin.notna()
    on_origin = origin_dist < params["origin_radius"]
    further_from_origin = not_first & (origin_dist > lag_origin)

    entered = not_first & on_origin & (~lag_on_origin)
    on_last_segment = pd.Series(
        np.cumsum((origin_dist > params["end_radius"]).to_numpy()[::-1])[::-1] == 0,
        index=group.index,
    )
    eligible = prop >= params["min_prop"]
    returned = pd.Series(np.cumsum((entered & on_last_segment & eligible).to_numpy()) > 0, index=group.index)
    missed = pd.Series(np.cumsum((returned & further_from_origin).to_numpy()) > 0, index=group.index)

    stopped = timediff.notna() & (timediff > params["min_pause"])
    near_end = (prop >= params["end_prop"]) & (lag_origin < params["pause_radius"])
    stopped_near = pd.Series(np.cumsum((stopped & near_end).fillna(False).to_numpy()) > 0, index=group.index)

    return (missed | stopped_near).fillna(False)


def apply_done_flags(responses: pd.DataFrame, params: dict[str, float]) -> pd.Series:
    """Apply ``trial_done_for_group`` across all trials.

    Args:
        responses: Current tracing samples with ``timediff``.
        params: Done-filter threshold dictionary.

    Returns:
        Boolean Series aligned to ``responses``.

    Side effects:
        None.
    """

    done = pd.Series(False, index=responses.index)
    for _, group in responses.groupby(TRIAL_KEY, sort=False, dropna=False):
        done.loc[group.index] = trial_done_for_group(group, params)
    return done


def summarize_failed_end_trials(responses: pd.DataFrame) -> pd.DataFrame:
    """Create the old-style failed trial-end audit table.

    Args:
        responses: Current tracing samples with a boolean ``done`` column.

    Returns:
        One row per trial with at least one ``done`` sample.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    for key_values, group in responses.groupby(TRIAL_KEY, sort=False, dropna=False):
        done = group["done"].fillna(False)
        num_done = int(done.sum())
        if num_done == 0:
            continue
        not_done_time = group.loc[~done, "time"]
        max_not_done = float(not_done_time.max()) if len(not_done_time) else np.nan
        max_time = float(group["time"].max())
        row = dict(zip(TRIAL_KEY, key_values, strict=True))
        row.update(
            {
                "samples": int(len(group)),
                "num_done": num_done,
                "mt_diff": max_not_done - max_time,
                "prop_remaining": 1.0 - (num_done / len(group)),
                "stage": "failed_end",
                "stage_order": STAGE_ORDER["failed_end"],
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def apply_failed_end_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Trim samples after failed trial-end attempts.

    Args:
        responses: Current tracing samples after no-shape trial removal.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)``.

    Side effects:
        Mutates ``flags`` with ``failed_end_done`` and removal-stage markers.
    """

    current = add_basic_trace_metrics(responses, zero_time=False)
    current["done1"] = apply_done_flags(current, DONE_FILTER_PARAMS)
    current["done2"] = apply_done_flags(current, DONE_FILTER_PARAMS_2)
    current["done"] = current["done1"] | current["done2"]

    done_rows = current.loc[current["done"]].copy()
    set_sample_flag(flags, done_rows, "failed_end_done")
    set_removed_stage(flags, done_rows, "failed_end")

    audit = summarize_failed_end_trials(current)
    kept = current.loc[~current["done"]].drop(columns=["done1", "done2", "done"]).copy()
    return kept, audit


def glitch_flags_for_group(group: pd.DataFrame) -> pd.Series:
    """Port the old ``is_glitch`` predicate for one trial.

    Args:
        group: One trial's current tracing samples with ``angle_diff``.

    Returns:
        Boolean Series aligned to ``group`` marking touchscreen glitch points.

    Side effects:
        None.
    """

    min_dist = GLITCH_FILTER_PARAMS["min_dist"]
    min_angle_diff = GLITCH_FILTER_PARAMS["min_angle_diff"]
    min_angle_diff_alt = GLITCH_FILTER_PARAMS["min_angle_diff_alt"]

    x = group["x"].astype("float64")
    y = group["y"].astype("float64")
    angle_diff = group["angle_diff"].astype("float64")
    origin_dist = group["origin.dist"].astype("float64")
    dx = x - x.shift(1)
    dy = y - y.shift(1)
    dist = line_length(x.shift(1), y.shift(1), x, y)

    eligible = dist.notna() & dist.shift(-1).notna() & angle_diff.shift(-1).notna()
    both_large = (dist > min_dist) & (dist.shift(-1) > min_dist)
    sharp_angle = angle_diff.shift(-1).abs() > min_angle_diff
    angle_glitch = (eligible & both_large & sharp_angle).fillna(False)

    if int(angle_glitch.sum()) > 1:
        temp = pd.DataFrame({**{column: group[column] for column in TRIAL_KEY}, "x": x, "y": y})
        temp_dx = pd.Series(dx, index=group.index)
        temp_dy = pd.Series(dy, index=group.index)
        dx2 = temp_dx + temp_dx.shift(-1)
        dy2 = temp_dy + temp_dy.shift(-1)
        theta = np.arctan2(dy2, dx2) - np.arctan2(temp_dy.shift(1), temp_dx.shift(1))
        theta = np.where(np.abs(theta) > np.pi, theta - (2 * np.pi) * np.sign(theta), theta)
        angle_diff2 = pd.Series((theta / np.pi) * 180.0, index=temp.index)

        pre_glitch = angle_diff2.shift(-1).notna() & angle_glitch.shift(-1).fillna(False)
        sharp_angle_alt = (~pre_glitch) | (angle_diff2.shift(-1).abs() > min_angle_diff_alt)
        angle_glitch = angle_glitch & (angle_glitch.shift(-2).fillna(False) | sharp_angle_alt.fillna(False))

    before = x.shift(-1).notna() & angle_glitch.shift(-1).fillna(False) & ((x == x.shift(-1)) | (y == y.shift(-1)))
    after = x.shift(1).notna() & angle_glitch.shift(1).fillna(False) & ((x == x.shift(1)) | (y == y.shift(1)))
    double_glitch = (before | after).fillna(False)

    start_glitch = dist.isna() & (dist.shift(-1) > min_dist) & (origin_dist > min_dist)
    end_glitch = dist.shift(-1).isna() & ((origin_dist - origin_dist.shift(1)) > min_dist)

    return (angle_glitch | double_glitch | start_glitch.fillna(False) | end_glitch.fillna(False)).fillna(False)


def apply_glitch_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag and remove touchscreen glitch points.

    Args:
        responses: Current tracing samples after failed trial-end trimming.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)``.

    Side effects:
        Mutates ``flags`` with ``glitch_xy``, ``glitch``, and removal-stage
        markers.
    """

    current = add_basic_trace_metrics(responses, zero_time=False)
    current["angle_diff"] = angle_diffs_degrees(current)
    current["glitch"] = False
    for _, group in current.groupby(TRIAL_KEY, sort=False, dropna=False):
        current.loc[group.index, "glitch"] = glitch_flags_for_group(group)

    current["glitch_xy"] = (current["x"] == 239) & (current["y"] == 1079)
    current["glitch"] = current["glitch"] | current["glitch_xy"]

    glitch_rows = current.loc[current["glitch"]].copy()
    set_sample_flag(flags, current.loc[current["glitch_xy"]], "glitch_xy")
    set_sample_flag(flags, glitch_rows, "glitch")
    set_removed_stage(flags, glitch_rows, "glitch")

    if current.empty:
        audit = pd.DataFrame()
    else:
        audit = current.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
            samples=("tracing_row_id", "size"),
            glitches=("glitch", "sum"),
            non_xy_glitches=("glitch", lambda values: int(values.sum())),
        ).reset_index()
        non_xy_counts = current.assign(non_xy_glitch=current["glitch"] & (~current["glitch_xy"])).groupby(
            TRIAL_KEY, sort=False, dropna=False
        )["non_xy_glitch"].sum().reset_index(name="non_xy_glitches")
        audit = audit.drop(columns=["non_xy_glitches"]).merge(non_xy_counts, on=TRIAL_KEY, how="left")
        audit = audit.loc[audit["glitches"] > 0].copy()
    audit["stage"] = "glitch"
    audit["stage_order"] = STAGE_ORDER["glitch"]

    kept = current.loc[~current["glitch"]].drop(columns=["glitch", "glitch_xy"]).copy()
    return kept, audit


def false_start_for_group(group: pd.DataFrame) -> pd.Series:
    """Port the old ``false_start`` predicate for one trial.

    Args:
        group: One trial's current tracing samples with ``timediff``.

    Returns:
        Boolean Series aligned to ``group`` marking samples before the recovered
        trace start.

    Side effects:
        None.
    """

    n_samples = len(group)
    if n_samples == 0:
        return pd.Series(False, index=group.index)

    prop = pd.Series(np.arange(1, n_samples + 1, dtype="float64") / n_samples, index=group.index)
    stopped = group["timediff"].notna() & (group["timediff"] > FALSE_START_PARAMS["min_pause"])
    near_origin = (prop <= FALSE_START_PARAMS["max_prop"]) & (group["origin.dist"] < FALSE_START_PARAMS["start_radius"])
    events = (stopped & near_origin).fillna(False)
    event_count = int(events.sum())
    false_start = pd.Series(np.cumsum(events.to_numpy()) < event_count, index=group.index)
    return false_start.fillna(False)


def apply_false_start_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag false-start samples and recover old ``start_shift`` values.

    Args:
        responses: Current tracing samples after glitch removal.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)``.

    Side effects:
        Mutates ``flags`` with ``false_start`` and removal-stage markers.
    """

    current = add_basic_trace_metrics(responses, zero_time=False)
    current["false_start"] = False
    for _, group in current.groupby(TRIAL_KEY, sort=False, dropna=False):
        current.loc[group.index, "false_start"] = false_start_for_group(group)

    false_rows = current.loc[current["false_start"]].copy()
    set_sample_flag(flags, false_rows, "false_start")
    set_removed_stage(flags, false_rows, "false_start")

    rows: list[dict[str, Any]] = []
    for key_values, group in current.groupby(TRIAL_KEY, sort=False, dropna=False):
        flagged = int(group["false_start"].sum())
        if flagged == 0:
            continue
        first_kept_time = group.loc[~group["false_start"], "time"]
        start_shift = float(first_kept_time.iloc[0]) if len(first_kept_time) else np.nan
        row = dict(zip(TRIAL_KEY, key_values, strict=True))
        row.update(
            {
                "samples": int(len(group)),
                "flagged": flagged,
                "start_shift": start_shift,
                "stage": "false_start",
                "stage_order": STAGE_ORDER["false_start"],
            }
        )
        rows.append(row)
    audit = pd.DataFrame(rows)

    kept = current.loc[~current["false_start"]].drop(columns=["false_start"]).copy()
    return kept, audit


def hand_noise_for_group(group: pd.DataFrame) -> pd.Series:
    """Port the old ``hand_noise`` predicate for one trial.

    Args:
        group: One trial's current tracing samples with distance, timing, angle,
            and origin-distance metrics.

    Returns:
        Boolean Series aligned to ``group`` marking hand-noise samples.

    Side effects:
        None.
    """

    params = HAND_NOISE_PARAMS
    x = group["x"].astype("float64")
    y = group["y"].astype("float64")
    timediff = group["timediff"].astype("float64")
    angle_diff = group["angle_diff"].astype("float64")
    origin_dist = group["origin.dist"].astype("float64")
    dist = line_length(x.shift(1), y.shift(1), x, y)

    sharp_turn = (angle_diff.abs() + angle_diff.shift(-1).abs()) > params["min_sharp_turnsum"]
    any_sharp_jumps = bool(((dist > params["min_dist_a"]) & sharp_turn).fillna(False).any())

    sorted_dist = dist.dropna().sort_values(ascending=False)
    second_largest = float(sorted_dist.iloc[1]) if len(sorted_dist) >= 2 else np.nan

    eligible = any_sharp_jumps & angle_diff.notna() & (timediff < params["max_timediff"])
    thresh1 = (dist > params["min_dist_a"]) & (angle_diff.abs() > params["min_angle_diff"])
    thresh2 = dist > params["min_dist_b"]
    large_enough = thresh1 | thresh2
    large_jump = (eligible & large_enough & (dist >= second_largest)).fillna(False)
    jumps = pd.Series(np.cumsum(large_jump.to_numpy()), index=group.index)
    is_noise = ((jumps % 2) == 1) & ((int(jumps.max()) if len(jumps) else 0) % 2 == 0)

    ends_near_origin = bool(len(origin_dist) and origin_dist.iloc[-1] < params["min_end_dist"])
    eligible_end = x.shift(1).notna() & (not ends_near_origin)
    away_from_origin = (eligible_end & ((origin_dist - origin_dist.shift(1)) > params["min_end_dist"])).fillna(False)
    away_cumsum = pd.Series(np.cumsum(away_from_origin.to_numpy()), index=group.index)
    end_noise = away_cumsum == max(1, int(away_from_origin.sum()))

    noise = (is_noise | end_noise).fillna(False)

    # The old R helper constrained candidate noise regions by bounding-box size
    # so broad missing-data regions would not be relabelled as hand noise.
    if bool(noise.any()):
        noise_width = float(x.loc[noise].max() - x.loc[noise].min())
        noise_height = float(y.loc[noise].max() - y.loc[noise].min())
        too_big = max(noise_width, noise_height) > params["max_size"]
        if too_big:
            noise = pd.Series(False, index=group.index)

    return noise


def apply_hand_noise_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag and remove hand-noise samples.

    Args:
        responses: Current tracing samples after false-start removal.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)``.

    Side effects:
        Mutates ``flags`` with ``hnoise`` and removal-stage markers.
    """

    current = add_basic_trace_metrics(responses, zero_time=False)
    current["angle_diff"] = angle_diffs_degrees(current)
    current["hnoise"] = False
    for _, group in current.groupby(TRIAL_KEY, sort=False, dropna=False):
        current.loc[group.index, "hnoise"] = hand_noise_for_group(group)

    noise_rows = current.loc[current["hnoise"]].copy()
    set_sample_flag(flags, noise_rows, "hnoise")
    set_removed_stage(flags, noise_rows, "hand_noise")

    if current.empty:
        audit = pd.DataFrame()
    else:
        audit = current.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
            samples=("tracing_row_id", "size"),
            flagged=("hnoise", "sum"),
        ).reset_index()
        audit = audit.loc[audit["flagged"] > 0].copy()
    audit["stage"] = "hand_noise"
    audit["stage_order"] = STAGE_ORDER["hand_noise"]

    kept = current.loc[~current["hnoise"]].drop(columns=["hnoise"]).copy()
    return kept, audit


def summarize_trace_incomplete_inputs(responses: pd.DataFrame) -> pd.DataFrame:
    """Summarize current traces for the incomplete-tracing filter.

    Args:
        responses: Current tracing samples after hand-noise removal.

    Returns:
        One row per trial with end gap, size, lateral shift, and path length.

    Side effects:
        None.
    """

    if responses.empty:
        return pd.DataFrame(columns=[*TRIAL_KEY, "end_gap", "max_size", "lat_shift", "trace_len"])

    current = add_basic_trace_metrics(responses, zero_time=False)
    rows: list[dict[str, Any]] = []
    for key_values, group in current.groupby(TRIAL_KEY, sort=False, dropna=False):
        width = float(group["x"].max() - group["x"].min())
        height = float(group["y"].max() - group["y"].min())
        max_size = max(width, height)
        lat_shift = (safe_divide(group["x"].min() - (SCREEN_RES_X / 2), width) + 0.5) * 2
        row = dict(zip(TRIAL_KEY, key_values, strict=True))
        row.update(
            {
                "samples": int(len(group)),
                "end_gap": float(abs(group["origin.dist"].iloc[-1] - group["origin.dist"].iloc[0])),
                "trace_width": width,
                "trace_height": height,
                "max_size": max_size,
                "lat_shift": float(np.asarray(lat_shift).item()),
                "trace_len": float(group["dist"].sum(skipna=True)),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def apply_incomplete_stage(
    responses: pd.DataFrame,
    frame_info: pd.DataFrame,
    flags: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag and drop trials with incomplete tracings.

    Args:
        responses: Current tracing samples after hand-noise removal.
        frame_info: Trial-level figure geometry summaries.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)``.

    Side effects:
        Mutates ``flags`` with ``incomplete_trial`` and removal-stage markers.
    """

    trace_info = summarize_trace_incomplete_inputs(responses)
    audit = trace_info.merge(
        frame_info[[*TRIAL_KEY, "fig_max_size", "fig_lat_shift", "fig_len"]],
        on=TRIAL_KEY,
        how="left",
    )
    audit["size_ratio"] = safe_divide(audit["fig_max_size"], audit["max_size"])
    audit["shift_diff"] = audit["lat_shift"].abs() - audit["fig_lat_shift"].abs()
    audit["len_ratio"] = safe_divide(audit["fig_len"], audit["trace_len"])
    audit["too_small"] = audit["size_ratio"] > INCOMPLETE_PARAMS["min_size_ratio"]
    audit["no_return"] = audit["end_gap"] > INCOMPLETE_PARAMS["min_end_gap"]
    audit["accidental_end"] = (audit["len_ratio"] > INCOMPLETE_PARAMS["min_length_ratio"]) & (
        audit["shift_diff"] > INCOMPLETE_PARAMS["min_shift_diff"]
    )
    audit["incomplete"] = audit["too_small"] | audit["no_return"] | audit["accidental_end"]
    audit = audit.loc[audit["incomplete"]].copy()
    audit["stage"] = "incomplete"
    audit["stage_order"] = STAGE_ORDER["incomplete"]

    affected_mask = trial_membership_mask(responses, audit)
    affected_rows = responses.loc[affected_mask].copy()
    set_trial_flag(flags, audit, "incomplete_trial")
    set_removed_stage(flags, affected_rows, "incomplete")
    kept = responses.loc[~affected_mask].copy()
    return kept, audit


def apply_large_gap_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Flag and drop trials with excessive time or distance gaps.

    Args:
        responses: Current tracing samples after incomplete-tracing removal.
        flags: Full tracing flag table to update.

    Returns:
        ``(kept_rows, audit_rows)``.

    Side effects:
        Mutates ``flags`` with sample ``gap`` and trial ``large_gap_trial``
        markers.
    """

    current = add_basic_trace_metrics(responses, zero_time=False)
    current["angle_diff"] = angle_diffs_degrees(current)
    grouped = current.groupby(TRIAL_KEY, sort=False, dropna=False)
    current["turn_sum"] = current["angle_diff"].abs() + grouped["angle_diff"].shift(-1).abs()

    eligible = current["dist"].notna() & (current["timediff"] > GAP_FILTER_PARAMS["min_timegap"])
    lift_t = current["timediff"] > GAP_FILTER_PARAMS["min_pause"]
    lift_dt = eligible & (current["dist"] > GAP_FILTER_PARAMS["min_dist_b"])
    lift_dat = eligible & (current["dist"] > GAP_FILTER_PARAMS["min_dist_c"]) & (
        current["turn_sum"] > GAP_FILTER_PARAMS["min_turnsum_c"]
    )
    lift_dat2 = eligible & (current["dist"] > GAP_FILTER_PARAMS["min_dist_d"]) & (
        current["turn_sum"] > GAP_FILTER_PARAMS["min_turnsum_d"]
    )
    current["gap"] = (lift_t | lift_dt | lift_dat | lift_dat2).fillna(False)

    gap_rows = current.loc[current["gap"]].copy()
    set_sample_flag(flags, gap_rows, "gap")

    if current.empty:
        audit = pd.DataFrame()
    else:
        audit = current.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
            samples=("tracing_row_id", "size"),
            longest_pause=("timediff", "max"),
            longest_jump=("dist", "max"),
            any_gaps=("gap", "any"),
        ).reset_index()
        audit = audit.loc[audit["any_gaps"]].copy()
    audit["stage"] = "large_gap"
    audit["stage_order"] = STAGE_ORDER["large_gap"]

    affected_mask = trial_membership_mask(current, audit)
    affected_rows = current.loc[affected_mask].copy()
    set_trial_flag(flags, audit, "large_gap_trial")
    set_removed_stage(flags, affected_rows, "large_gap")
    kept = current.loc[~affected_mask].drop(columns=["gap"]).copy()
    return kept, audit


def apply_edge_stage(responses: pd.DataFrame, flags: pd.DataFrame) -> pd.DataFrame:
    """Flag remaining trials where a tracing touches a screen edge.

    Args:
        responses: Current tracing samples after large-gap trial removal.
        flags: Full tracing flag table to update.

    Returns:
        Trial-level edge-hit audit table. Rows are not removed.

    Side effects:
        Mutates ``flags`` with ``hit_edge_trial`` markers.
    """

    if responses.empty:
        audit = pd.DataFrame(
            columns=[*TRIAL_KEY, "hit_top_edge", "hit_bottom_edge", "hit_left_edge", "hit_right_edge"]
        )
    else:
        audit = responses.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
            hit_top_edge=("y", lambda values: bool((values == 0).any())),
            hit_bottom_edge=("y", lambda values: bool((values == (SCREEN_RES_Y - 1)).any())),
            hit_left_edge=("x", lambda values: bool((values == 0).any())),
            hit_right_edge=("x", lambda values: bool((values == (SCREEN_RES_X - 1)).any())),
        ).reset_index()
        audit = audit.loc[
            audit[["hit_top_edge", "hit_bottom_edge", "hit_left_edge", "hit_right_edge"]].any(axis=1)
        ].copy()

    audit["stage"] = "hit_edge"
    audit["stage_order"] = STAGE_ORDER["hit_edge"]
    set_trial_flag(flags, audit, "hit_edge_trial")
    return audit


def finalize_filtered_tracings(responses: pd.DataFrame, flags: pd.DataFrame) -> pd.DataFrame:
    """Recompute final filtered tracing metrics after all removal stages.

    Args:
        responses: Current samples after trial and sample filters.
        flags: Full tracing flag table to update with retained-row markers.

    Returns:
        Final filtered tracing table with final zeroed time and derived metrics.

    Side effects:
        Mutates ``flags`` by marking rows that remain after all removal stages.
    """

    if responses.empty:
        return responses.copy()

    final = sort_by_trial_and_order(responses, "sample_index").copy()
    final["time_from_original_trace_start"] = final["time"]
    grouped = final.groupby(TRIAL_KEY, sort=False, dropna=False)
    final["time"] = final["time"] - grouped["time"].transform("first")
    final = add_basic_trace_metrics(final, zero_time=False)
    final["angle_diff"] = angle_diffs_degrees(final)

    set_sample_flag(flags, final, "kept_after_filters")
    flags.loc[final["tracing_row_id"].to_numpy(), "filter_stage_removed"] = "kept_after_filters"
    flags.loc[final["tracing_row_id"].to_numpy(), "filter_stage_removed_order"] = pd.NA
    return final


def summarize_final_traces(filtered: pd.DataFrame) -> pd.DataFrame:
    """Summarize final filtered tracing duration and path length.

    Args:
        filtered: Final filtered tracing table.

    Returns:
        One row per retained tracing trial with ``mt_clip``, ``PLresp``, and
        ``vresp``.

    Side effects:
        None.
    """

    if filtered.empty:
        return pd.DataFrame(columns=[*TRIAL_KEY, "final_tracing_rows", "mt_clip", "PLresp", "vresp"])

    summary = filtered.groupby(TRIAL_KEY, sort=False, dropna=False).agg(
        final_tracing_rows=("tracing_row_id", "size"),
        mt_clip=("time", "max"),
        PLresp=("dist", "sum"),
        first_kept_time_from_original_trace_start=("time_from_original_trace_start", "min"),
        last_kept_time_from_original_trace_start=("time_from_original_trace_start", "max"),
    ).reset_index()
    summary["vresp"] = safe_divide(summary["PLresp"], summary["mt_clip"])
    return summary


def no_tracing_payload_trials(zip_manifest: pd.DataFrame) -> pd.DataFrame:
    """Create an audit table for zip trials with no non-NA tracing samples.

    Args:
        zip_manifest: One-row-per-zip parse manifest from reconstruction.

    Returns:
        Trial-level rows where the tracing payload had no parsed samples.

    Side effects:
        None.
    """

    audit = zip_manifest.loc[
        (zip_manifest["tracings_parse_status"] != "ok") | (zip_manifest["tracings_row_count"].fillna(0) == 0)
    ].copy()
    columns = [
        *TRIAL_KEY,
        "source_zip_basename",
        "issue_codes",
        "tracings_parse_status",
        "tracings_row_count",
    ]
    audit = audit[[column for column in columns if column in audit.columns]]
    audit["stage"] = "no_tracing_payload"
    audit["stage_order"] = STAGE_ORDER["no_tracing_payload"]
    return audit


def run_filter_pipeline(tables: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Run all tracing-filter stages and collect outputs.

    Args:
        tables: Input tables returned by ``read_input_tables``.

    Returns:
        Dictionary containing filtered tracing outputs, trial summaries, and
        stage audit tables.

    Side effects:
        None. File writing happens later.
    """

    points = tables["points"]
    frames = tables["frames"]
    tracings = tables["tracings"]
    zip_manifest = tables["zip_manifest"]

    origins = build_origins(points)
    frame_info = build_frame_info(frames)
    raw_tracing_info = build_raw_tracing_info(tracings)
    responses = prepare_response_samples(tracings, origins)
    flags = initialize_flag_table(responses)

    stage_audits: dict[str, pd.DataFrame] = {}

    print("Applying repeated-point filter")
    responses, stage_audits["repeated_point"] = apply_repeated_point_stage(responses, flags)

    print("Applying no-shape trial filter")
    responses, stage_audits["no_shape"] = apply_no_shape_stage(responses, frame_info, flags)

    print("Applying failed trial-end trimming")
    responses, stage_audits["failed_end"] = apply_failed_end_stage(responses, flags)

    print("Applying touchscreen glitch filter")
    responses, stage_audits["glitch"] = apply_glitch_stage(responses, flags)

    print("Applying false-start filter")
    responses, stage_audits["false_start"] = apply_false_start_stage(responses, flags)

    print("Applying hand-noise filter")
    responses, stage_audits["hand_noise"] = apply_hand_noise_stage(responses, flags)

    print("Applying incomplete-tracing trial filter")
    responses, stage_audits["incomplete"] = apply_incomplete_stage(responses, frame_info, flags)

    print("Applying excessive-gap trial filter")
    responses, stage_audits["large_gap"] = apply_large_gap_stage(responses, flags)

    print("Flagging edge-hit trials")
    stage_audits["hit_edge"] = apply_edge_stage(responses, flags)
    stage_audits["no_tracing_payload"] = no_tracing_payload_trials(zip_manifest)

    filtered = finalize_filtered_tracings(responses, flags)
    final_summary = summarize_final_traces(filtered)
    trial_summary = build_trial_filter_summary(
        zip_manifest=zip_manifest,
        raw_tracing_info=raw_tracing_info,
        final_summary=final_summary,
        stage_audits=stage_audits,
    )
    trace_filter_info = build_trace_filter_info(stage_audits)

    outputs: dict[str, pd.DataFrame] = {
        "tracings_with_filter_flags": flags.reset_index(drop=True),
        "filtered_tracings": filtered.reset_index(drop=True),
        "tracing_trial_filter_summary": trial_summary,
        "trace_filter_info": trace_filter_info,
    }
    outputs.update({f"audit_{name}": audit for name, audit in stage_audits.items()})
    return outputs


def audit_count_by_trial(audit: pd.DataFrame, value_column: str, output_column: str) -> pd.DataFrame:
    """Return a compact trial-key audit count table.

    Args:
        audit: Stage audit table.
        value_column: Column to copy into the compact output.
        output_column: Desired output column name.

    Returns:
        Trial-key table with one renamed value column, or an empty typed table.

    Side effects:
        None.
    """

    if audit.empty or value_column not in audit.columns:
        return pd.DataFrame(columns=[*TRIAL_KEY, output_column])
    return audit[[*TRIAL_KEY, value_column]].rename(columns={value_column: output_column})


def build_trial_filter_summary(
    zip_manifest: pd.DataFrame,
    raw_tracing_info: pd.DataFrame,
    final_summary: pd.DataFrame,
    stage_audits: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Build one row per zip/trial with filter-stage facts.

    Args:
        zip_manifest: Reconstruction parse manifest.
        raw_tracing_info: Raw tracing sample summary.
        final_summary: Final filtered tracing summary.
        stage_audits: Stage audit tables.

    Returns:
        Trial-level filter summary suitable for later event-offset work.

    Side effects:
        None.
    """

    manifest_columns = [
        *TRIAL_KEY,
        "source_zip_basename",
        "issue_codes",
        "frames_parse_status",
        "frames_row_count",
        "tracings_parse_status",
        "tracings_row_count",
    ]
    summary = zip_manifest[[column for column in manifest_columns if column in zip_manifest.columns]].copy()
    summary = summary.merge(raw_tracing_info, on=TRIAL_KEY, how="left")
    summary = summary.merge(final_summary, on=TRIAL_KEY, how="left")

    compact_tables = [
        audit_count_by_trial(stage_audits["repeated_point"].groupby(TRIAL_KEY, dropna=False).size().reset_index(name="repeated_point_rows"), "repeated_point_rows", "repeated_point_rows"),
        audit_count_by_trial(stage_audits["failed_end"], "num_done", "failed_end_rows"),
        audit_count_by_trial(stage_audits["glitch"], "glitches", "glitch_rows"),
        audit_count_by_trial(stage_audits["false_start"], "flagged", "false_start_rows"),
        audit_count_by_trial(stage_audits["false_start"], "start_shift", "start_shift"),
        audit_count_by_trial(stage_audits["hand_noise"], "flagged", "hand_noise_rows"),
    ]
    for compact in compact_tables:
        summary = summary.merge(compact, on=TRIAL_KEY, how="left")

    for stage_name, flag_column in [
        ("no_shape", "no_shape_trial"),
        ("incomplete", "incomplete_trial"),
        ("large_gap", "large_gap_trial"),
        ("hit_edge", "hit_edge_trial"),
        ("no_tracing_payload", "no_tracing_payload"),
    ]:
        audit = stage_audits[stage_name]
        flag_table = audit[TRIAL_KEY].drop_duplicates().copy() if not audit.empty else pd.DataFrame(columns=TRIAL_KEY)
        flag_table[flag_column] = True
        summary = summary.merge(flag_table, on=TRIAL_KEY, how="left")
        summary[flag_column] = summary[flag_column].fillna(False)

    count_columns = [
        "input_tracing_rows",
        "final_tracing_rows",
        "repeated_point_rows",
        "failed_end_rows",
        "glitch_rows",
        "false_start_rows",
        "hand_noise_rows",
    ]
    for column in count_columns:
        if column in summary.columns:
            summary[column] = summary[column].fillna(0).astype("int64")

    summary["has_reconstructed_tracing"] = summary["input_tracing_rows"] > 0
    summary["has_filtered_tracing"] = summary["final_tracing_rows"] > 0
    return summary.sort_values(TRIAL_KEY, kind="mergesort").reset_index(drop=True)


def build_trace_filter_info(stage_audits: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Stack stage audit tables into a long old-compatible audit table.

    Args:
        stage_audits: Audit tables produced by each filter stage.

    Returns:
        Long table with ``stage`` and ``stage_order`` columns.

    Side effects:
        None.
    """

    pieces: list[pd.DataFrame] = []
    for stage_name, audit in stage_audits.items():
        piece = audit.copy()
        if "stage" not in piece.columns:
            piece["stage"] = stage_name
        if "stage_order" not in piece.columns:
            piece["stage_order"] = STAGE_ORDER[stage_name]
        piece["old_trace_filter_info_name"] = stage_name
        pieces.append(piece)

    if not pieces:
        return pd.DataFrame(columns=[*TRIAL_KEY, "stage", "stage_order", "old_trace_filter_info_name"])

    combined = pd.concat(pieces, ignore_index=True, sort=False)
    sort_columns = ["stage_order", *TRIAL_KEY]
    existing_sort_columns = [column for column in sort_columns if column in combined.columns]
    return combined.sort_values(existing_sort_columns, kind="mergesort").reset_index(drop=True)


def stage_counts(outputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Calculate row and trial counts for each audit stage.

    Args:
        outputs: Output dictionary from ``run_filter_pipeline``.

    Returns:
        Table with stage names, audit row counts, and affected trial counts.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    for stage_name in STAGE_OUTPUT_FILES:
        audit = outputs.get(f"audit_{stage_name}", pd.DataFrame())
        affected_trials = (
            int(audit[TRIAL_KEY].drop_duplicates().shape[0])
            if not audit.empty and set(TRIAL_KEY).issubset(audit.columns)
            else 0
        )
        rows.append(
            {
                "stage": stage_name,
                "stage_order": STAGE_ORDER[stage_name],
                "audit_rows": int(len(audit)),
                "affected_trials": affected_trials,
            }
        )
    return pd.DataFrame(rows).sort_values("stage_order").reset_index(drop=True)


def watchlist_summary(trial_summary: pd.DataFrame) -> list[str]:
    """Format local review-watchlist facts for the Markdown summary.

    Args:
        trial_summary: Trial-level summary table.

    Returns:
        Markdown bullet lines.

    Side effects:
        None.
    """

    lines: list[str] = []
    present_ids = set(int(value) for value in trial_summary["participant_id"].dropna().unique())
    filtered_ids = set(
        int(value)
        for value in trial_summary.loc[trial_summary["has_filtered_tracing"], "participant_id"].dropna().unique()
    )

    for label, ids in REVIEW_WATCHLIST_IDS.items():
        present = [value for value in ids if value in present_ids]
        final = [value for value in ids if value in filtered_ids]
        absent = [value for value in ids if value not in present_ids]
        lines.append(
            f"- {label}: present in reconstructed zip manifest: {present or 'none'}; "
            f"with final filtered tracing rows: {final or 'none'}; "
            f"not present in reconstructed zip manifest: {absent or 'none'}"
        )
    return lines


def write_outputs(outputs: dict[str, pd.DataFrame], repo_root: Path, started_at: datetime) -> None:
    """Write parquet outputs and the local Markdown summary.

    Args:
        outputs: Output dictionary from ``run_filter_pipeline``.
        repo_root: Repository root used to resolve output paths.
        started_at: Start timestamp for run metadata.

    Returns:
        None.

    Side effects:
        Creates ``_Data/behavior/tracing_filters/`` if needed and writes local
        parquet/Markdown outputs there.
    """

    require_pyarrow()
    output_dir = repo_root / FILTER_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    main_outputs = {
        "tracings_with_filter_flags": "tracings_with_filter_flags.parquet",
        "filtered_tracings": "filtered_tracings.parquet",
        "tracing_trial_filter_summary": "tracing_trial_filter_summary.parquet",
        "trace_filter_info": "trace_filter_info.parquet",
    }

    for output_name, filename in main_outputs.items():
        path = output_dir / filename
        outputs[output_name].to_parquet(path, index=False, engine="pyarrow")
        print(f"Wrote {relative_to_repo(path, repo_root)}")

    for stage_name, filename in STAGE_OUTPUT_FILES.items():
        path = output_dir / filename
        outputs[f"audit_{stage_name}"].to_parquet(path, index=False, engine="pyarrow")
        print(f"Wrote {relative_to_repo(path, repo_root)}")

    summary_path = output_dir / "tracing_filter_summary.md"
    summary_path.write_text(build_markdown_summary(outputs, started_at), encoding="utf-8")
    print(f"Wrote {relative_to_repo(summary_path, repo_root)}")


def build_markdown_summary(outputs: dict[str, pd.DataFrame], started_at: datetime) -> str:
    """Build the local Markdown summary text.

    Args:
        outputs: Output dictionary from ``run_filter_pipeline``.
        started_at: Start timestamp for run metadata.

    Returns:
        Markdown summary text.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    flags = outputs["tracings_with_filter_flags"]
    filtered = outputs["filtered_tracings"]
    trial_summary = outputs["tracing_trial_filter_summary"]
    counts = stage_counts(outputs)

    input_rows = int(len(flags))
    input_trials = int(flags[TRIAL_KEY].drop_duplicates().shape[0]) if input_rows else 0
    output_rows = int(len(filtered))
    output_trials = int(filtered[TRIAL_KEY].drop_duplicates().shape[0]) if output_rows else 0

    stage_lines = [
        f"- {row.stage}: audit rows {int(row.audit_rows):,}; affected trials {int(row.affected_trials):,}"
        for row in counts.itertuples(index=False)
    ]
    watch_lines = watchlist_summary(trial_summary)

    removed_stage_counts = Counter(flags["filter_stage_removed"].fillna(""))
    removed_lines = [
        f"- {stage}: {count:,} rows"
        for stage, count in sorted(removed_stage_counts.items(), key=lambda item: (str(item[0]) == "", str(item[0])))
        if stage
    ]
    if not removed_lines:
        removed_lines = ["- none"]

    return f"""# Tracing Filter Summary

Generated: {finished_at.isoformat(timespec="seconds")}

Started: {started_at.isoformat(timespec="seconds")}

## Scope

This local summary reports the Python port of the historical DEMI TraceLab
tracing filters. It reads reconstructed TraceLab tables from
`_Data/behavior/tracing_tables/` and writes audit outputs under
`_Data/behavior/tracing_filters/`. It does not join to task data or EEG data,
inspect raw EEG, preprocess EEG, or make participant-level decisions.

## Row and Trial Counts

- Input tracing rows: {input_rows:,}
- Input tracing trials with non-NA samples: {input_trials:,}
- Output filtered tracing rows: {output_rows:,}
- Output filtered tracing trials: {output_trials:,}
- Zip/trial rows in trial summary: {len(trial_summary):,}

## Filter Stage Audit Counts

{chr(10).join(stage_lines)}

## First Removal Stage by Row

{chr(10).join(removed_lines)}

## Local Review Watchlist

These IDs are surfaced from prior local inventory reviews for human review only.
They are not decision rules.

{chr(10).join(watch_lines)}

## Known Deviations and TODOs for Exact Parity

- The old R plotting side effects are not reproduced here.
- Accuracy/error metrics, interpolation, and task joins are deferred to later scripts.
- This port uses explicit `sample_index` ordering before every lag/lead operation; old R relied on grouped row order.
- The next parity check should compare these audit tables against the historical R `trace_filter_info.rds` where available.
- Later event-offset work should use `start_shift`, raw trace onset/end, and final `mt_clip` from these local outputs.
"""


def main() -> int:
    """Run the tracing-filter port from reconstructed TraceLab tables.

    Args:
        None.

    Returns:
        Process exit code: ``0`` on success, ``1`` on a handled failure.

    Side effects:
        Reads local reconstructed TraceLab tables and writes local filter
        outputs below ``_Data/behavior/tracing_filters/``.
    """

    started_at = datetime.now()
    repo_root = repo_root_from_script()

    try:
        tables = read_input_tables(repo_root)
        outputs = run_filter_pipeline(tables)
        write_outputs(outputs, repo_root, started_at)
    except RuntimeError as error:
        print(f"FAIL tracing filters: {error}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
