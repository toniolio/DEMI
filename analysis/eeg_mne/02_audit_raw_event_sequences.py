"""Audit raw EEG annotation sequences with historical trigger-cleanup logic.

This script is the second read-only raw event audit for the DEMI EEG
reanalysis. It starts from the local CSV outputs produced by
``analysis/eeg_mne/01_compare_raw_annotations_to_offsets.py`` and compares a
proposed, in-memory event-label sequence against the old-compatible
behavioural/tracing offset table at
``_Data/behavior/event_offsets/event_offsets_old_compatible.csv``.

The second audit is needed because the first raw annotation comparison showed
that the EDF annotations are readable and mostly map to the old trigger
vocabulary, but that the raw trial alignment is not stable enough for epoch
construction. Several files have missing, duplicated, mislabeled, split-file,
concatenated-file, raw-only, or offset-only event structures. This script keeps
those facts visible while asking a narrower technical question: how much of the
current mismatch surface is explained by the historical duplicate/mislabeled
trigger cleanup logic if that logic is evaluated in memory only?

Historical logic ported here:

- ``external/DEMI_EEG_Pipeline/eeg_pipeline.py`` searched semantic annotation
  labels for consecutive repeats. If the repeated label was followed within
  2 ms by the next annotation, the repeated row was treated as a duplicated
  trigger and removed from the cleaned annotation sequence. If the next
  annotation was farther away, the repeated row was treated as a mislabeled
  trigger and advanced to the next event name in the old event order:
  ``stim_on``, ``red_on``, ``trace_start``, ``trace_end``,
  ``accuracy_submit``, ``vividness_submit``. This script ports that rule into
  proposed labels only; it never rewrites an EDF.
- ``_Scripts/_functions/eeg.R::update_events`` counted trials from ``stim_on``
  triggers and bad-start-like ``red_on`` triggers, marked practice-like groups
  using the old ``red_on - stim_on > 3000 ms`` rule, removed practice-like
  groups before old trial renumbering, and then removed bad-start-like groups
  from the joinable event stream. This script rebuilds those trial numbers
  from the proposed labels.
- The same R function printed timing gut checks when the observed
  ``trace_end`` timing differed from ``trace_start + end_trigger`` by more
  than 150 ms, or when ``red_on - stim_on`` differed from task
  ``stimulus_mt`` by more than 250 ms. This script preserves those checks and
  adds nearby-trial timing candidates so count-level joins can be compared with
  timing similarity.

Intentionally deferred:

- No raw EDF is opened, edited, cropped, or saved.
- No EEG signal is loaded.
- No preprocessing, filtering, epoching, interpolation, ICA, or time-frequency
  step is run.
- No repaired EEG file, BIDS file, or annotation file is written.
- No split EDFs are concatenated; IDs 54, 56, and 65 are reported separately.
- The concatenated ID 5 EDF is reported as found; the historical crop-at
  ``file start`` step is not applied here.
- No participant-level retention or scientific event-handling decision is made.
  The script surfaces facts and candidate timing evidence for Tony's review.

Expected inputs, all read from local repo paths:

- ``_Data/eeg/event_alignment/raw_annotation_events.csv``;
- ``_Data/eeg/event_alignment/annotation_trial_summary.csv``;
- ``_Data/eeg/event_alignment/offset_annotation_join_audit.csv``;
- ``_Data/eeg/event_alignment/raw_annotation_alignment_summary.md``;
- ``_Data/behavior/event_offsets/event_offsets_old_compatible.csv``.

Generated local-only outputs:

- ``_Data/eeg/event_sequence_audit/proposed_annotation_events.csv``;
- ``_Data/eeg/event_sequence_audit/proposed_trial_summary.csv``;
- ``_Data/eeg/event_sequence_audit/proposed_offset_join_audit.csv``;
- ``_Data/eeg/event_sequence_audit/first_mismatch_by_file.csv``;
- ``_Data/eeg/event_sequence_audit/nearby_trial_candidate_audit.csv``;
- ``_Data/eeg/event_sequence_audit/event_sequence_audit_summary.md``.

Safety boundaries:

- All outputs are written below ``_Data/eeg/event_sequence_audit/``.
- The script uses ``pathlib`` and ``pandas`` and does not import MNE.
- Existing annotation CSVs and old-compatible offsets are read only.
- Proposed labels are audit columns, not edits to source data.
- Failure messages name missing inputs or columns directly so the audit fails
  before producing partial, misleading outputs.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/02_audit_raw_event_sequences.py

The script also works when launched from another directory because paths are
resolved relative to this file's location in the repository.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


EVENT_ALIGNMENT_DIR = Path("_Data") / "eeg" / "event_alignment"
EVENT_SEQUENCE_AUDIT_DIR = Path("_Data") / "eeg" / "event_sequence_audit"
OLD_COMPATIBLE_OFFSETS_PATH = Path("_Data") / "behavior" / "event_offsets" / "event_offsets_old_compatible.csv"

RAW_ANNOTATION_EVENTS_FILENAME = "raw_annotation_events.csv"
ANNOTATION_TRIAL_SUMMARY_FILENAME = "annotation_trial_summary.csv"
OFFSET_ANNOTATION_JOIN_AUDIT_FILENAME = "offset_annotation_join_audit.csv"
ALIGNMENT_SUMMARY_FILENAME = "raw_annotation_alignment_summary.md"

PROPOSED_ANNOTATION_EVENTS_FILENAME = "proposed_annotation_events.csv"
PROPOSED_TRIAL_SUMMARY_FILENAME = "proposed_trial_summary.csv"
PROPOSED_OFFSET_JOIN_AUDIT_FILENAME = "proposed_offset_join_audit.csv"
FIRST_MISMATCH_BY_FILE_FILENAME = "first_mismatch_by_file.csv"
NEARBY_TRIAL_CANDIDATE_AUDIT_FILENAME = "nearby_trial_candidate_audit.csv"
EVENT_SEQUENCE_AUDIT_SUMMARY_FILENAME = "event_sequence_audit_summary.md"

LOW_ANNOTATION_COUNT_THRESHOLD = 600
PRACTICE_FIG_DURATION_SECONDS_THRESHOLD = 3.0
TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD = 0.150
STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD = 0.250
HISTORICAL_DUPLICATE_TRIGGER_SECONDS = 0.002
NEARBY_TRIAL_WINDOW = 5

WATCHLIST_IDS = {5, 11, 14, 54, 56, 65, 86, 89, 94, 100}
SPLIT_EDF_IDS = {54, 56, 65}
CONCATENATED_EDF_IDS = {5}
FACTUAL_SPECIAL_CASE_IDS = {11, 14, 89, 94, 100}

EXPECTED_TRIAL_EVENT_NAMES = ("stim_on", "red_on", "trace_start", "trace_end")
TASK_EVENT_NAMES = (
    "stim_on",
    "red_on",
    "trace_start",
    "trace_end",
    "accuracy_submit",
    "vividness_submit",
)
HISTORICAL_CLEANUP_EVENT_ORDER = (
    "stim_on",
    "red_on",
    "trace_start",
    "trace_end",
    "accuracy_submit",
    "vividness_submit",
)

INTEGER_AUDIT_COLUMNS = (
    "participant_id",
    "split_part",
    "event_code",
    "raw_trial_sequence",
    "old_update_trial_count",
    "join_trial_count",
    "proposed_raw_trial_sequence",
    "proposed_old_update_trial_count",
    "proposed_join_trial_count",
    "audit_participant_id",
    "audit_trial_count",
    "offset_id",
    "offset_trial_count",
    "trial_count",
    "candidate_offset_trial_count",
    "candidate_offset_delta",
    "candidate_timing_rank",
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Args:
        None.

    Returns:
        Absolute path to the repository root.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def text_or_empty(value: Any) -> str:
    """Return a clean text value, preserving missing values as empty strings.

    Args:
        value: Any scalar value read from a CSV cell.

    Returns:
        ``""`` when the value is missing, otherwise stripped string text.

    Side effects:
        None.
    """

    if pd.isna(value):
        return ""
    return str(value).strip()


def is_true(value: Any) -> bool:
    """Interpret booleans that may have been read back from CSV text.

    Args:
        value: Scalar value that may be a bool, number, string, or missing
            value.

    Returns:
        ``True`` for true-like values and ``False`` otherwise.

    Side effects:
        None.
    """

    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def safe_float(value: Any) -> float:
    """Convert a scalar to float, using ``nan`` for missing/non-numeric values.

    Args:
        value: Scalar value to convert.

    Returns:
        Floating-point value or ``numpy.nan``.

    Side effects:
        None.
    """

    try:
        if pd.isna(value):
            return np.nan
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def int_or_zero(value: Any) -> int:
    """Convert a scalar count to int, using zero for missing values.

    Args:
        value: Scalar value expected to represent a non-negative count.

    Returns:
        Integer count, with missing/non-numeric values treated as zero.

    Side effects:
        None.
    """

    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return 0
    return int(numeric)


def first_nonmissing(values: pd.Series) -> Any:
    """Return the first non-missing value from a Series.

    Args:
        values: Values to scan in existing order.

    Returns:
        First non-missing scalar, or ``numpy.nan`` if every value is missing.

    Side effects:
        None.
    """

    for value in values:
        if not pd.isna(value):
            return value
    return np.nan


def unique_joined_text(values: pd.Series) -> str:
    """Join unique non-empty text values in their first-seen order.

    Args:
        values: Series containing text-like values.

    Returns:
        Semicolon-delimited text.

    Side effects:
        None.
    """

    seen: dict[str, None] = {}
    for value in values:
        text = text_or_empty(value)
        if text:
            seen[text] = None
    return ";".join(seen.keys())


def split_issue_codes(value: Any) -> list[str]:
    """Split a semicolon-delimited issue-code cell into individual codes.

    Args:
        value: Issue-code cell value.

    Returns:
        List of non-empty codes in cell order.

    Side effects:
        None.
    """

    return [code for code in text_or_empty(value).split(";") if code]


def count_issue_codes(values: pd.Series) -> Counter[str]:
    """Count semicolon-delimited issue codes.

    Args:
        values: Series of issue-code strings.

    Returns:
        Counter keyed by issue code.

    Side effects:
        None.
    """

    counts: Counter[str] = Counter()
    for value in values.fillna(""):
        for code in split_issue_codes(value):
            counts[code] += 1
    return counts


def coerce_nullable_integer_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert known integer-like audit columns to nullable integer dtype.

    Args:
        frame: DataFrame with columns that may contain integer-like values and
            missing values.

    Returns:
        Copy with known integer columns coerced to pandas ``Int64`` where
        possible.

    Side effects:
        None.
    """

    out = frame.copy()
    for column in INTEGER_AUDIT_COLUMNS:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
    return out


def read_csv_checked(path: Path, required_columns: set[str], label: str) -> pd.DataFrame:
    """Read a CSV file and verify required columns.

    Args:
        path: CSV path to read.
        required_columns: Column names that must be present.
        label: Human-readable input label for error messages.

    Returns:
        DataFrame read from ``path``.

    Side effects:
        Reads one CSV file. Raises ``RuntimeError`` when the file or required
        columns are missing.
    """

    if not path.exists():
        raise RuntimeError(f"required {label} is missing: {path.as_posix()}")

    frame = pd.read_csv(path, low_memory=False)
    missing = sorted(required_columns.difference(frame.columns))
    if missing:
        raise RuntimeError(f"{label} is missing required column(s): {', '.join(missing)}")
    return frame


def read_old_compatible_offsets(path: Path) -> pd.DataFrame:
    """Read and validate the old-compatible behavioural/tracing offsets.

    Args:
        path: Path to ``event_offsets_old_compatible.csv``.

    Returns:
        Offset DataFrame sorted by participant ID and trial count.

    Side effects:
        Reads a CSV file. Raises ``RuntimeError`` for missing columns or
        duplicate ``id``/``trial_count`` keys.
    """

    required_columns = {
        "id",
        "physical",
        "it",
        "mt",
        "stimulus_mt",
        "trial_count",
        "trace_onset",
        "trace_end",
        "start_shift",
        "real_start",
        "real_end",
        "end_trigger",
    }
    offsets = read_csv_checked(path, required_columns, "old-compatible event-offset CSV")
    offsets = offsets.copy()
    offsets["id"] = pd.to_numeric(offsets["id"], errors="coerce").astype("Int64")
    offsets["trial_count"] = pd.to_numeric(offsets["trial_count"], errors="coerce").astype("Int64")

    duplicate_keys = offsets.duplicated(["id", "trial_count"], keep=False)
    if duplicate_keys.any():
        raise RuntimeError(
            "old-compatible offsets contain "
            f"{int(duplicate_keys.sum())} duplicate id/trial_count row(s)"
        )

    return offsets.sort_values(["id", "trial_count"], kind="mergesort").reset_index(drop=True)


def read_audit_inputs(repo_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """Read the first-audit outputs and old-compatible event offsets.

    Args:
        repo_root: Repository root.

    Returns:
        ``(raw_events, first_trial_summary, first_join_audit, offsets,
        first_summary_text)``.

    Side effects:
        Reads local CSV/Markdown inputs below ``_Data/``.
    """

    alignment_dir = repo_root / EVENT_ALIGNMENT_DIR
    raw_events = read_csv_checked(
        alignment_dir / RAW_ANNOTATION_EVENTS_FILENAME,
        {
            "participant_id",
            "source_filename",
            "file_role",
            "split_part",
            "event_row_order",
            "onset_seconds",
            "annotation_description",
            "event_name",
            "mapped_task_event",
            "trial_numbering_candidate",
        },
        "raw annotation event CSV from script 01",
    )
    first_trial_summary = read_csv_checked(
        alignment_dir / ANNOTATION_TRIAL_SUMMARY_FILENAME,
        {
            "summary_row_kind",
            "participant_id",
            "source_filename",
            "file_role",
            "split_part",
            "annotation_count",
            "mapped_task_event_count",
            "trial_issue_codes",
        },
        "annotation trial-summary CSV from script 01",
    )
    first_join_audit = read_csv_checked(
        alignment_dir / OFFSET_ANNOTATION_JOIN_AUDIT_FILENAME,
        {
            "audit_participant_id",
            "audit_trial_count",
            "join_status",
            "old_trace_epoch_duration_mismatch",
            "old_stimulus_duration_mismatch",
            "missing_expected_event_names",
        },
        "offset annotation join-audit CSV from script 01",
    )
    offsets = read_old_compatible_offsets(repo_root / OLD_COMPATIBLE_OFFSETS_PATH)

    summary_path = alignment_dir / ALIGNMENT_SUMMARY_FILENAME
    first_summary_text = summary_path.read_text(encoding="utf-8") if summary_path.exists() else ""

    return (
        coerce_nullable_integer_columns(raw_events),
        coerce_nullable_integer_columns(first_trial_summary),
        coerce_nullable_integer_columns(first_join_audit),
        offsets,
        first_summary_text,
    )


def build_file_summary_from_first_trial_summary(first_trial_summary: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct one-row-per-file facts from script-01 trial summary rows.

    Args:
        first_trial_summary: ``annotation_trial_summary.csv`` from the first
            raw annotation comparison.

    Returns:
        DataFrame with one row per EDF file and file-level annotation facts.

    Side effects:
        None.
    """

    required_columns = [
        "participant_id",
        "participant_id_padded",
        "source_filename",
        "file_path",
        "file_role",
        "split_part",
        "watchlist_id",
        "read_status",
        "sampling_frequency_hz",
        "recording_duration_seconds",
        "annotation_count",
        "mapped_task_event_count",
        "dominant_event_code_pattern",
        "trial_issue_codes",
    ]
    available_columns = [column for column in required_columns if column in first_trial_summary.columns]
    if "source_filename" not in available_columns:
        raise RuntimeError("first trial summary cannot be reduced to file facts without source_filename")

    sort_columns = [column for column in ["participant_id", "source_filename", "raw_trial_sequence"] if column in first_trial_summary]
    sorted_summary = first_trial_summary.sort_values(sort_columns, kind="mergesort", na_position="last")
    file_summary = sorted_summary[available_columns].drop_duplicates("source_filename", keep="first").copy()
    file_summary = file_summary.rename(columns={"trial_issue_codes": "file_issue_codes"})
    if "watchlist_id" not in file_summary.columns:
        file_summary["watchlist_id"] = file_summary["participant_id"].map(
            lambda value: int(value) in WATCHLIST_IDS if pd.notna(value) else False
        )
    return coerce_nullable_integer_columns(
        file_summary.sort_values(["participant_id", "source_filename"], kind="mergesort", na_position="last")
        .reset_index(drop=True)
    )


def next_historical_event_name(event_name: str) -> str | None:
    """Return the next label in the old preprocessing event-name order.

    Args:
        event_name: Current event name.

    Returns:
        Next historical event name, or ``None`` if ``event_name`` is not in the
        old order or is already the terminal old label.

    Side effects:
        None.
    """

    if event_name not in HISTORICAL_CLEANUP_EVENT_ORDER:
        return None
    index = HISTORICAL_CLEANUP_EVENT_ORDER.index(event_name)
    if index + 1 >= len(HISTORICAL_CLEANUP_EVENT_ORDER):
        return None
    return HISTORICAL_CLEANUP_EVENT_ORDER[index + 1]


def propose_annotation_event_sequence(raw_events: pd.DataFrame) -> pd.DataFrame:
    """Apply historical duplicate/mislabeled-trigger rules in memory.

    Args:
        raw_events: One-row-per-annotation table from script 01. The table must
            already contain the derived semantic ``event_name`` field.

    Returns:
        Copy of ``raw_events`` with raw labels preserved and proposed label,
        action, reason, and trial-numbering candidate columns added.

    Side effects:
        None. No EDF or annotation file is written.
    """

    if raw_events.empty:
        return raw_events.copy()

    out = raw_events.copy()
    out["raw_event_name"] = out["event_name"].map(text_or_empty)
    out["raw_annotation_description"] = out["annotation_description"].map(text_or_empty)
    out["proposed_event_name"] = np.where(
        out["raw_event_name"].isin(HISTORICAL_CLEANUP_EVENT_ORDER),
        out["raw_event_name"],
        "",
    )
    out["proposed_sequence_action"] = np.where(out["proposed_event_name"].eq(""), "skipped", "unchanged")
    out["proposed_sequence_reason"] = np.where(
        out["proposed_event_name"].eq(""),
        "not_in_historical_cleanup_event_order",
        "raw_label_preserved",
    )

    file_groups = out.groupby(["participant_id", "source_filename"], dropna=False, sort=False).groups
    for _, index in file_groups.items():
        file_events = out.loc[index].sort_values("event_row_order", kind="mergesort")
        ordered_indices = list(file_events.index)

        # The historical Python preprocessing loop checked interior annotation
        # rows only: range(1, trigger_count - 1). This preserves that exact
        # scope so the audit can distinguish a ported old rule from any future
        # deliberate adaptation Tony might choose.
        for position in range(1, len(ordered_indices) - 1):
            current_index = ordered_indices[position]
            previous_index = ordered_indices[position - 1]
            next_index = ordered_indices[position + 1]

            current_name = text_or_empty(out.at[current_index, "raw_event_name"])
            previous_name = text_or_empty(out.at[previous_index, "raw_event_name"])
            if not current_name or current_name != previous_name:
                continue

            current_onset = safe_float(out.at[current_index, "onset_seconds"])
            next_onset = safe_float(out.at[next_index, "onset_seconds"])
            next_gap_seconds = next_onset - current_onset

            # Old logic treated rapid same-label repeats as doubled triggers.
            # This audit keeps the row but removes it from the proposed trial
            # sequence by clearing proposed_event_name.
            if np.isfinite(next_gap_seconds) and next_gap_seconds < HISTORICAL_DUPLICATE_TRIGGER_SECONDS:
                out.at[current_index, "proposed_event_name"] = ""
                out.at[current_index, "proposed_sequence_action"] = "duplicated"
                out.at[current_index, "proposed_sequence_reason"] = (
                    "historical_same_label_repeat_with_next_gap_under_2ms"
                )
                continue

            next_name = next_historical_event_name(current_name)
            if next_name is None:
                out.at[current_index, "proposed_sequence_action"] = "unresolved"
                out.at[current_index, "proposed_sequence_reason"] = (
                    "same_label_repeat_has_no_next_historical_event_name"
                )
                continue

            # Old logic treated slower same-label repeats as the next expected
            # event mislabeled as the previous event. The proposed label is
            # advanced, while raw_event_name remains unchanged for review.
            out.at[current_index, "proposed_event_name"] = next_name
            out.at[current_index, "proposed_sequence_action"] = "relabeled"
            out.at[current_index, "proposed_sequence_reason"] = (
                f"historical_same_label_repeat_relabel_to_{next_name}"
            )

    out["proposed_mapped_task_event"] = out["proposed_event_name"].isin(TASK_EVENT_NAMES)
    out["proposed_trial_numbering_candidate"] = (
        out["proposed_mapped_task_event"] & out["proposed_sequence_action"].isin(["unchanged", "relabeled"])
    )
    out["proposed_label_differs_from_raw"] = out["proposed_event_name"].ne(out["raw_event_name"])

    for column, default in [
        ("proposed_raw_trial_sequence", pd.NA),
        ("proposed_old_update_trial_count", pd.NA),
        ("proposed_join_trial_count", pd.NA),
        ("proposed_old_update_removed_reason", ""),
        ("proposed_fig_duration_seconds", np.nan),
        ("proposed_bad_start_like_event", False),
        ("proposed_practice_like_group", False),
        ("proposed_bad_start_like_group", False),
    ]:
        out[column] = default

    return coerce_nullable_integer_columns(out)


def assign_proposed_trial_numbers(proposed_events: pd.DataFrame) -> pd.DataFrame:
    """Apply old ``update_events`` trial-numbering order to proposed labels.

    Args:
        proposed_events: Annotation table after in-memory duplicate/relabel
            audit columns have been added.

    Returns:
        Copy with proposed raw trial sequence, old post-practice trial number,
        join trial number, and practice/bad-start flags filled.

    Side effects:
        None.
    """

    if proposed_events.empty:
        return proposed_events.copy()

    out = proposed_events.copy()
    group_columns = ["participant_id", "source_filename"]

    for _, index in out.groupby(group_columns, dropna=False, sort=False).groups.items():
        file_events = out.loc[index].sort_values("event_row_order", kind="mergesort").copy()
        candidate = file_events[file_events["proposed_trial_numbering_candidate"].fillna(False)].copy()
        if candidate.empty:
            continue

        candidate["original_event_index"] = candidate.index
        candidate["previous_proposed_task_event_name"] = candidate["proposed_event_name"].shift()
        candidate["previous_proposed_task_onset_seconds"] = candidate["onset_seconds"].shift()

        # This mirrors update_events: a red_on that is not immediately preceded
        # by stim_on is treated as a bad-start-like event and starts a group.
        # The group is retained through practice renumbering but will not get a
        # join trial count.
        candidate["proposed_bad_start_like_event"] = candidate["proposed_event_name"].eq("red_on") & ~candidate[
            "previous_proposed_task_event_name"
        ].eq("stim_on")
        starts_new_trial = candidate["proposed_event_name"].eq("stim_on") | candidate[
            "proposed_bad_start_like_event"
        ]
        candidate["proposed_raw_trial_sequence"] = starts_new_trial.cumsum()

        candidate["duration_since_previous_proposed_task_event_seconds"] = (
            candidate["onset_seconds"] - candidate["previous_proposed_task_onset_seconds"]
        )
        candidate.loc[
            candidate["proposed_bad_start_like_event"],
            "duration_since_previous_proposed_task_event_seconds",
        ] = np.nan

        trial_rows: list[dict[str, Any]] = []
        for proposed_raw_trial_sequence, trial_events in candidate.groupby("proposed_raw_trial_sequence", sort=True):
            red_durations = trial_events.loc[
                trial_events["proposed_event_name"].eq("red_on"),
                "duration_since_previous_proposed_task_event_seconds",
            ].dropna()
            fig_duration_seconds = float(red_durations.iloc[0]) if not red_durations.empty else np.nan
            practice_like_group = (
                not pd.isna(fig_duration_seconds)
                and fig_duration_seconds > PRACTICE_FIG_DURATION_SECONDS_THRESHOLD
            )
            bad_start_like_group = bool(trial_events["proposed_bad_start_like_event"].any())
            trial_rows.append(
                {
                    "proposed_raw_trial_sequence": int(proposed_raw_trial_sequence),
                    "proposed_fig_duration_seconds": fig_duration_seconds,
                    "proposed_practice_like_group": practice_like_group,
                    "proposed_bad_start_like_group": bad_start_like_group,
                }
            )

        trial_info = pd.DataFrame(trial_rows)

        # Old ordering matters: practice-like groups are removed first, then
        # remaining groups are renumbered, then bad-start-like groups are
        # omitted from the joinable stream. Later groups keep the old
        # post-practice numbering even when a bad-start group occurred.
        nonpractice = trial_info[~trial_info["proposed_practice_like_group"]].copy()
        nonpractice["proposed_old_update_trial_count"] = np.arange(1, len(nonpractice) + 1, dtype="int64")
        trial_info = trial_info.merge(
            nonpractice[["proposed_raw_trial_sequence", "proposed_old_update_trial_count"]],
            on="proposed_raw_trial_sequence",
            how="left",
        )
        trial_info["proposed_join_trial_count"] = trial_info["proposed_old_update_trial_count"]
        trial_info.loc[trial_info["proposed_bad_start_like_group"], "proposed_join_trial_count"] = np.nan

        trial_info["proposed_old_update_removed_reason"] = ""
        trial_info.loc[
            trial_info["proposed_practice_like_group"], "proposed_old_update_removed_reason"
        ] = "practice_like_group"
        trial_info.loc[
            trial_info["proposed_bad_start_like_group"], "proposed_old_update_removed_reason"
        ] = "bad_start_like_group"
        trial_info.loc[
            trial_info["proposed_practice_like_group"] & trial_info["proposed_bad_start_like_group"],
            "proposed_old_update_removed_reason",
        ] = "practice_like_group;bad_start_like_group"

        original_index = candidate["original_event_index"]
        trial_info_by_sequence = trial_info.set_index("proposed_raw_trial_sequence")
        for column in [
            "proposed_raw_trial_sequence",
            "proposed_old_update_trial_count",
            "proposed_join_trial_count",
            "proposed_old_update_removed_reason",
            "proposed_fig_duration_seconds",
            "proposed_practice_like_group",
            "proposed_bad_start_like_group",
            "proposed_bad_start_like_event",
        ]:
            if column in {"proposed_raw_trial_sequence", "proposed_bad_start_like_event"}:
                values = candidate[column]
            else:
                values = candidate["proposed_raw_trial_sequence"].map(trial_info_by_sequence[column])
            out.loc[original_index, column] = values.to_numpy()

    return coerce_nullable_integer_columns(out)


def proposed_trial_issue_codes(row: dict[str, Any]) -> list[str]:
    """Create per-trial issue/status codes for the proposed sequence surface.

    Args:
        row: Proposed trial-summary row dictionary.

    Returns:
        List of short issue/status codes.

    Side effects:
        None.
    """

    codes: list[str] = []
    if is_true(row.get("practice_like_group", False)):
        codes.append("practice_like_group")
    if is_true(row.get("bad_start_like_group", False)):
        codes.append("bad_start_like_group")
    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        if row.get(f"missing_{event_name}", False):
            codes.append(f"missing_{event_name}")
        if int_or_zero(row.get(f"{event_name}_count", 0)) > 1:
            codes.append(f"extra_{event_name}")
    if int_or_zero(row.get("relabeled_event_count", 0)) > 0:
        codes.append("proposed_relabel_present")
    if row.get("file_role") == "split_part":
        codes.append("split_file_not_concatenated")
    if row.get("file_role") == "concatenated":
        codes.append("concatenated_file_not_repaired")
    if is_true(row.get("watchlist_id", False)):
        codes.append("watchlist_id")
    return codes


def summarize_proposed_trial_events(trial_events: pd.DataFrame) -> dict[str, Any]:
    """Summarize one proposed-label trial group.

    Args:
        trial_events: Proposed task-event rows from one EDF and proposed raw
            trial group.

    Returns:
        Dictionary containing proposed event counts, key event onsets, missing
        event flags, timing fields, and proposed cleanup action counts.

    Side effects:
        None.
    """

    counts = Counter(trial_events["proposed_event_name"])
    row: dict[str, Any] = {
        "sequence_surface": "proposed",
        "summary_row_kind": "trial_group",
        "participant_id": first_nonmissing(trial_events["participant_id"]),
        "participant_id_padded": first_nonmissing(trial_events["participant_id_padded"]),
        "source_filename": first_nonmissing(trial_events["source_filename"]),
        "file_path": first_nonmissing(trial_events["file_path"]),
        "file_role": first_nonmissing(trial_events["file_role"]),
        "split_part": first_nonmissing(trial_events["split_part"]),
        "watchlist_id": is_true(first_nonmissing(trial_events["watchlist_id"])),
        "read_status": first_nonmissing(trial_events["read_status"]),
        "sampling_frequency_hz": first_nonmissing(trial_events["sampling_frequency_hz"]),
        "recording_duration_seconds": first_nonmissing(trial_events["recording_duration_seconds"]),
        "annotation_count": np.nan,
        "mapped_task_event_count": int(len(trial_events)),
        "dominant_event_code_pattern": "",
        "raw_trial_sequence": first_nonmissing(trial_events["proposed_raw_trial_sequence"]),
        "old_update_trial_count": first_nonmissing(trial_events["proposed_old_update_trial_count"]),
        "join_trial_count": first_nonmissing(trial_events["proposed_join_trial_count"]),
        "proposed_raw_trial_sequence": first_nonmissing(trial_events["proposed_raw_trial_sequence"]),
        "proposed_old_update_trial_count": first_nonmissing(trial_events["proposed_old_update_trial_count"]),
        "proposed_join_trial_count": first_nonmissing(trial_events["proposed_join_trial_count"]),
        "old_update_removed_reason": unique_joined_text(trial_events["proposed_old_update_removed_reason"]),
        "practice_like_group": is_true(first_nonmissing(trial_events["proposed_practice_like_group"])),
        "bad_start_like_group": is_true(first_nonmissing(trial_events["proposed_bad_start_like_group"])),
        "fig_duration_seconds": first_nonmissing(trial_events["proposed_fig_duration_seconds"]),
        "trial_event_count": int(len(trial_events)),
        "relabeled_event_count": int(trial_events["proposed_sequence_action"].eq("relabeled").sum()),
        "duplicated_event_count": int(trial_events["proposed_sequence_action"].eq("duplicated").sum()),
        "unresolved_event_count": int(trial_events["proposed_sequence_action"].eq("unresolved").sum()),
        "proposed_sequence_actions": unique_joined_text(trial_events["proposed_sequence_action"]),
        "proposed_sequence_reasons": unique_joined_text(trial_events["proposed_sequence_reason"]),
    }

    for event_name in TASK_EVENT_NAMES:
        event_rows = trial_events[trial_events["proposed_event_name"].eq(event_name)]
        row[f"{event_name}_count"] = int(counts.get(event_name, 0))
        row[f"{event_name}_onset_seconds"] = first_nonmissing(event_rows["onset_seconds"]) if not event_rows.empty else np.nan

    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        row[f"missing_{event_name}"] = row[f"{event_name}_count"] == 0

    missing_expected = [event_name for event_name in EXPECTED_TRIAL_EVENT_NAMES if row[f"missing_{event_name}"]]
    row["missing_expected_event_names"] = ";".join(missing_expected)
    extra_expected = [event_name for event_name in EXPECTED_TRIAL_EVENT_NAMES if row[f"{event_name}_count"] > 1]
    row["extra_expected_event_names"] = ";".join(extra_expected)

    trace_start = row["trace_start_onset_seconds"]
    trace_end = row["trace_end_onset_seconds"]
    if pd.notna(trace_start) and pd.notna(trace_end):
        row["trace_epoch_duration_seconds"] = float(trace_end) - float(trace_start)
    else:
        row["trace_epoch_duration_seconds"] = np.nan

    row["trial_issue_codes"] = ";".join(proposed_trial_issue_codes(row))
    return row


def file_without_proposed_trials_summary_row(file_row: pd.Series) -> dict[str, Any]:
    """Create a proposed trial-summary row for a file with no trial groups.

    Args:
        file_row: One row from the reconstructed file summary.

    Returns:
        Dictionary matching the proposed trial-summary schema as closely as
        possible, with trial fields left missing.

    Side effects:
        None.
    """

    row: dict[str, Any] = {
        "sequence_surface": "proposed",
        "summary_row_kind": "file_without_proposed_trial_groups",
        "participant_id": file_row["participant_id"],
        "participant_id_padded": file_row.get("participant_id_padded", ""),
        "source_filename": file_row["source_filename"],
        "file_path": file_row.get("file_path", ""),
        "file_role": file_row["file_role"],
        "split_part": file_row["split_part"],
        "watchlist_id": is_true(file_row.get("watchlist_id", False)),
        "read_status": file_row.get("read_status", ""),
        "sampling_frequency_hz": file_row.get("sampling_frequency_hz", np.nan),
        "recording_duration_seconds": file_row.get("recording_duration_seconds", np.nan),
        "annotation_count": file_row.get("annotation_count", np.nan),
        "mapped_task_event_count": 0,
        "dominant_event_code_pattern": file_row.get("dominant_event_code_pattern", ""),
        "raw_trial_sequence": pd.NA,
        "old_update_trial_count": pd.NA,
        "join_trial_count": pd.NA,
        "proposed_raw_trial_sequence": pd.NA,
        "proposed_old_update_trial_count": pd.NA,
        "proposed_join_trial_count": pd.NA,
        "old_update_removed_reason": "",
        "practice_like_group": False,
        "bad_start_like_group": False,
        "fig_duration_seconds": np.nan,
        "trial_event_count": 0,
        "relabeled_event_count": 0,
        "duplicated_event_count": 0,
        "unresolved_event_count": 0,
        "proposed_sequence_actions": "",
        "proposed_sequence_reasons": "",
        "trace_epoch_duration_seconds": np.nan,
        "missing_expected_event_names": ";".join(EXPECTED_TRIAL_EVENT_NAMES),
        "extra_expected_event_names": "",
    }
    for event_name in TASK_EVENT_NAMES:
        row[f"{event_name}_count"] = 0
        row[f"{event_name}_onset_seconds"] = np.nan
    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        row[f"missing_{event_name}"] = True

    inherited_codes = split_issue_codes(file_row.get("file_issue_codes", ""))
    inherited_codes.append("no_proposed_trial_groups_for_file")
    row["trial_issue_codes"] = ";".join(dict.fromkeys(code for code in inherited_codes if code))
    return row


def build_proposed_trial_summary(proposed_events: pd.DataFrame, file_summary: pd.DataFrame) -> pd.DataFrame:
    """Build one proposed-label trial-summary table.

    Args:
        proposed_events: Annotation event table after proposed labels and
            proposed old trial numbers have been assigned.
        file_summary: One-row-per-EDF file facts recovered from script 01.

    Returns:
        DataFrame with one row per proposed trial group plus factual file rows
        for EDFs without proposed trial groups.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    if not proposed_events.empty:
        candidates = proposed_events[proposed_events["proposed_trial_numbering_candidate"].fillna(False)].copy()
        candidates = candidates[candidates["proposed_raw_trial_sequence"].notna()]
        for _, trial_events in candidates.groupby(
            ["participant_id", "source_filename", "proposed_raw_trial_sequence"],
            dropna=False,
            sort=True,
        ):
            rows.append(summarize_proposed_trial_events(trial_events.sort_values("event_row_order", kind="mergesort")))

    files_with_trial_rows = {row["source_filename"] for row in rows}
    for _, file_row in file_summary.iterrows():
        if file_row["source_filename"] not in files_with_trial_rows:
            rows.append(file_without_proposed_trials_summary_row(file_row))

    trial_summary = pd.DataFrame(rows)
    if trial_summary.empty:
        return trial_summary

    for _, file_row in file_summary.iterrows():
        mask = trial_summary["source_filename"].eq(file_row["source_filename"])
        for column in [
            "annotation_count",
            "dominant_event_code_pattern",
            "recording_duration_seconds",
            "sampling_frequency_hz",
        ]:
            if column in file_row.index:
                trial_summary.loc[mask, column] = file_row[column]

    sort_columns = ["participant_id", "source_filename", "raw_trial_sequence"]
    trial_summary = trial_summary.sort_values(sort_columns, kind="mergesort", na_position="last")
    return coerce_nullable_integer_columns(trial_summary.reset_index(drop=True))


def participant_file_fact_maps(file_summary: pd.DataFrame) -> dict[str, set[int]]:
    """Create participant-level file fact maps for join audit issue codes.

    Args:
        file_summary: One-row-per-EDF file facts.

    Returns:
        Dictionary mapping fact names to participant ID sets.

    Side effects:
        None.
    """

    maps = {
        "has_raw_edf_file": set(),
        "has_zero_annotation_file": set(),
        "has_low_annotation_file": set(),
        "has_split_file": set(),
        "has_concatenated_file": set(),
    }
    for _, row in file_summary.iterrows():
        participant_id = row.get("participant_id")
        if pd.isna(participant_id):
            continue
        pid = int(participant_id)
        annotation_count = int_or_zero(row.get("annotation_count", 0))
        maps["has_raw_edf_file"].add(pid)
        if annotation_count == 0:
            maps["has_zero_annotation_file"].add(pid)
        if 0 < annotation_count < LOW_ANNOTATION_COUNT_THRESHOLD:
            maps["has_low_annotation_file"].add(pid)
        if row.get("file_role") == "split_part":
            maps["has_split_file"].add(pid)
        if row.get("file_role") == "concatenated":
            maps["has_concatenated_file"].add(pid)
    return maps


def join_issue_codes(row: pd.Series) -> str:
    """Build semicolon-delimited issue codes for one proposed join row.

    Args:
        row: Row from the full outer proposed raw/offset join audit.

    Returns:
        Semicolon-delimited issue/status codes.

    Side effects:
        None.
    """

    codes: list[str] = []
    if row.get("join_status") == "raw_annotation_without_offset":
        codes.append("raw_annotation_without_offset")
    if row.get("join_status") == "offset_without_raw_annotation_trial":
        codes.append("offset_without_raw_annotation_trial")
    if is_true(row.get("raw_join_key_duplicate", False)):
        codes.append("duplicate_raw_annotation_join_key")
    if is_true(row.get("old_trace_epoch_duration_mismatch", False)):
        codes.append("old_trace_epoch_duration_mismatch")
    if is_true(row.get("old_stimulus_duration_mismatch", False)):
        codes.append("old_stimulus_duration_mismatch")
    if int_or_zero(row.get("relabeled_event_count", 0)) > 0:
        codes.append("proposed_relabel_present")

    for event_name in split_issue_codes(row.get("missing_expected_event_names", "")):
        codes.append(f"missing_{event_name}")
    for event_name in split_issue_codes(row.get("extra_expected_event_names", "")):
        codes.append(f"extra_{event_name}")

    if is_true(row.get("participant_has_zero_annotation_file", False)):
        codes.append("participant_has_zero_annotation_file")
    if is_true(row.get("participant_has_low_annotation_file", False)):
        codes.append("participant_has_low_annotation_file")
    if is_true(row.get("participant_has_split_file", False)):
        codes.append("split_file_not_concatenated")
    if is_true(row.get("participant_has_concatenated_file", False)):
        codes.append("concatenated_file_not_repaired")
    if not is_true(row.get("participant_has_raw_edf_file", True)):
        codes.append("participant_has_no_raw_edf_file")
    if is_true(row.get("watchlist_id", False)):
        codes.append("watchlist_id")

    return ";".join(dict.fromkeys(code for code in codes if code))


def alignment_problem_codes(row: pd.Series) -> str:
    """Return issue codes that represent direct alignment problems.

    Args:
        row: Proposed join-audit row.

    Returns:
        Semicolon-delimited codes focused on event sequence, join, and timing
        problems. General watchlist/split/concatenated facts are not counted
        here unless they also create a direct join or timing problem.

    Side effects:
        None.
    """

    codes: list[str] = []
    if row.get("join_status") == "raw_annotation_without_offset":
        codes.append("raw_annotation_without_offset")
    if row.get("join_status") == "offset_without_raw_annotation_trial":
        codes.append("offset_without_raw_annotation_trial")
    if is_true(row.get("raw_join_key_duplicate", False)):
        codes.append("duplicate_raw_annotation_join_key")
    if is_true(row.get("old_trace_epoch_duration_mismatch", False)):
        codes.append("old_trace_epoch_duration_mismatch")
    if is_true(row.get("old_stimulus_duration_mismatch", False)):
        codes.append("old_stimulus_duration_mismatch")
    for event_name in split_issue_codes(row.get("missing_expected_event_names", "")):
        codes.append(f"missing_{event_name}")
    for event_name in split_issue_codes(row.get("extra_expected_event_names", "")):
        codes.append(f"extra_{event_name}")
    return ";".join(dict.fromkeys(code for code in codes if code))


def build_proposed_offset_join_audit(
    proposed_trial_summary: pd.DataFrame,
    offsets: pd.DataFrame,
    file_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Compare proposed annotation-derived trials with old offsets.

    Args:
        proposed_trial_summary: Per-trial summary rebuilt from proposed labels.
        offsets: Old-compatible behavioural/tracing offsets.
        file_summary: One-row-per-EDF facts for participant-level file context.

    Returns:
        Full outer join audit keyed by participant ID and trial count.

    Side effects:
        None.
    """

    raw_joinable = proposed_trial_summary[
        proposed_trial_summary["summary_row_kind"].eq("trial_group")
        & proposed_trial_summary["join_trial_count"].notna()
    ].copy()

    raw_joinable["participant_id"] = pd.to_numeric(raw_joinable["participant_id"], errors="coerce").astype("Int64")
    raw_joinable["join_trial_count"] = pd.to_numeric(raw_joinable["join_trial_count"], errors="coerce").astype("Int64")
    raw_joinable["raw_join_key_duplicate"] = raw_joinable.duplicated(
        ["participant_id", "join_trial_count"], keep=False
    )

    raw_columns = [
        "participant_id",
        "participant_id_padded",
        "source_filename",
        "file_path",
        "file_role",
        "split_part",
        "watchlist_id",
        "summary_row_kind",
        "raw_trial_sequence",
        "old_update_trial_count",
        "join_trial_count",
        "proposed_raw_trial_sequence",
        "proposed_old_update_trial_count",
        "proposed_join_trial_count",
        "raw_join_key_duplicate",
        "practice_like_group",
        "bad_start_like_group",
        "fig_duration_seconds",
        "trial_event_count",
        "relabeled_event_count",
        "duplicated_event_count",
        "unresolved_event_count",
        "proposed_sequence_actions",
        "proposed_sequence_reasons",
        "stim_on_count",
        "red_on_count",
        "trace_start_count",
        "trace_end_count",
        "accuracy_submit_count",
        "vividness_submit_count",
        "stim_on_onset_seconds",
        "red_on_onset_seconds",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
        "accuracy_submit_onset_seconds",
        "vividness_submit_onset_seconds",
        "trace_epoch_duration_seconds",
        "missing_expected_event_names",
        "extra_expected_event_names",
        "trial_issue_codes",
    ]
    raw_joinable = raw_joinable[[column for column in raw_columns if column in raw_joinable.columns]].copy()

    offsets_for_join = offsets.rename(columns={"id": "offset_id", "trial_count": "offset_trial_count"})

    # As in the first audit, the full outer merge keeps both directions visible:
    # annotation-derived trials with no old-compatible offset and offsets with
    # no annotation-derived trial.
    joined = raw_joinable.merge(
        offsets_for_join,
        left_on=["participant_id", "join_trial_count"],
        right_on=["offset_id", "offset_trial_count"],
        how="outer",
        validate="many_to_one",
        indicator=True,
    )

    joined["audit_participant_id"] = joined["participant_id"].combine_first(joined["offset_id"]).astype("Int64")
    joined["audit_trial_count"] = joined["join_trial_count"].combine_first(joined["offset_trial_count"]).astype("Int64")
    joined["join_status"] = joined["_merge"].map(
        {
            "both": "raw_annotation_and_offset",
            "left_only": "raw_annotation_without_offset",
            "right_only": "offset_without_raw_annotation_trial",
        }
    )

    # Timing gut check 1 from update_events: observed trace_end should be close
    # to trace_start + end_trigger. This is a timing-similarity check, not a
    # count check.
    joined["trace_end_expected_from_offset_seconds"] = joined["trace_start_onset_seconds"] + joined["end_trigger"]
    joined["trace_epoch_end_difference_seconds"] = (
        joined["trace_end_onset_seconds"] - joined["trace_end_expected_from_offset_seconds"]
    ).abs()
    joined["old_trace_epoch_duration_mismatch"] = (
        joined["trace_epoch_end_difference_seconds"] > TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD
    ).fillna(False)

    # Timing gut check 2 from update_events: red_on - stim_on should be close
    # to task stimulus_mt. Large differences often mark event-sequence drift.
    joined["stimulus_duration_difference_seconds"] = joined["fig_duration_seconds"] - joined["stimulus_mt"]
    joined["old_stimulus_duration_mismatch"] = (
        joined["stimulus_duration_difference_seconds"].abs() > STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD
    ).fillna(False)

    fact_maps = participant_file_fact_maps(file_summary)
    joined["participant_has_raw_edf_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_raw_edf_file"] if pd.notna(value) else False
    )
    joined["participant_has_zero_annotation_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_zero_annotation_file"] if pd.notna(value) else False
    )
    joined["participant_has_low_annotation_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_low_annotation_file"] if pd.notna(value) else False
    )
    joined["participant_has_split_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_split_file"] if pd.notna(value) else False
    )
    joined["participant_has_concatenated_file"] = joined["audit_participant_id"].map(
        lambda value: int(value) in fact_maps["has_concatenated_file"] if pd.notna(value) else False
    )
    joined["watchlist_id"] = joined["audit_participant_id"].map(
        lambda value: int(value) in WATCHLIST_IDS if pd.notna(value) else False
    )
    joined["join_audit_issue_codes"] = joined.apply(join_issue_codes, axis=1)
    joined["alignment_problem_codes"] = joined.apply(alignment_problem_codes, axis=1)
    joined["possible_alignment_problem"] = joined["alignment_problem_codes"].ne("")
    joined["clean_timing_row"] = (
        joined["join_status"].eq("raw_annotation_and_offset")
        & ~joined["old_trace_epoch_duration_mismatch"].fillna(False)
        & ~joined["old_stimulus_duration_mismatch"].fillna(False)
        & joined["missing_expected_event_names"].fillna("").eq("")
        & joined["extra_expected_event_names"].fillna("").eq("")
        & ~joined["raw_join_key_duplicate"].fillna(False)
    )

    preferred = [
        "audit_participant_id",
        "audit_trial_count",
        "join_status",
        "join_audit_issue_codes",
        "alignment_problem_codes",
        "possible_alignment_problem",
        "clean_timing_row",
        "source_filename",
        "file_role",
        "split_part",
        "watchlist_id",
        "participant_has_raw_edf_file",
        "participant_has_zero_annotation_file",
        "participant_has_low_annotation_file",
        "participant_has_split_file",
        "participant_has_concatenated_file",
        "raw_trial_sequence",
        "old_update_trial_count",
        "join_trial_count",
        "proposed_raw_trial_sequence",
        "proposed_old_update_trial_count",
        "proposed_join_trial_count",
        "offset_id",
        "offset_trial_count",
        "raw_join_key_duplicate",
        "fig_duration_seconds",
        "stimulus_mt",
        "stimulus_duration_difference_seconds",
        "old_stimulus_duration_mismatch",
        "trace_start_onset_seconds",
        "trace_end_onset_seconds",
        "trace_epoch_duration_seconds",
        "end_trigger",
        "trace_end_expected_from_offset_seconds",
        "trace_epoch_end_difference_seconds",
        "old_trace_epoch_duration_mismatch",
        "physical",
        "it",
        "mt",
        "trace_onset",
        "trace_end",
        "start_shift",
        "real_start",
        "real_end",
        "trial_event_count",
        "relabeled_event_count",
        "duplicated_event_count",
        "unresolved_event_count",
        "proposed_sequence_actions",
        "proposed_sequence_reasons",
        "stim_on_count",
        "red_on_count",
        "trace_start_count",
        "trace_end_count",
        "accuracy_submit_count",
        "vividness_submit_count",
        "missing_expected_event_names",
        "extra_expected_event_names",
        "trial_issue_codes",
    ]
    ordered = [column for column in preferred if column in joined.columns] + [
        column for column in joined.columns if column not in preferred and column != "_merge"
    ]
    return coerce_nullable_integer_columns(
        joined[ordered].sort_values(
            ["audit_participant_id", "audit_trial_count", "source_filename"],
            kind="mergesort",
            na_position="last",
        )
    ).reset_index(drop=True)


def candidate_timing_score(trace_difference: float, stimulus_difference: float) -> float:
    """Calculate a simple nearby-candidate timing score.

    Args:
        trace_difference: Absolute trace-end timing difference in seconds.
        stimulus_difference: Absolute stimulus-duration difference in seconds.

    Returns:
        Sum of available absolute differences, or ``nan`` if neither timing
        difference is available.

    Side effects:
        None.
    """

    parts = [value for value in [trace_difference, stimulus_difference] if np.isfinite(value)]
    return float(sum(parts)) if parts else np.nan


def build_nearby_trial_candidate_audit(join_audit: pd.DataFrame, offsets: pd.DataFrame) -> pd.DataFrame:
    """Build timing-similarity candidates around problematic nominal joins.

    Args:
        join_audit: Proposed full outer join audit.
        offsets: Old-compatible behavioural/tracing offsets.

    Returns:
        DataFrame with nearby offset trial candidates for proposed raw rows that
        have direct alignment problem codes.

    Side effects:
        None.
    """

    columns = [
        "problem_row_number",
        "audit_participant_id",
        "source_filename",
        "file_role",
        "split_part",
        "raw_trial_sequence",
        "join_trial_count",
        "join_status",
        "alignment_problem_codes",
        "candidate_offset_trial_count",
        "candidate_offset_delta",
        "candidate_trace_end_difference_seconds",
        "candidate_stimulus_difference_seconds",
        "candidate_timing_score_seconds",
        "candidate_trace_within_old_threshold",
        "candidate_stimulus_within_old_threshold",
        "candidate_clean_by_timing",
        "nominal_candidate_timing_score_seconds",
        "candidate_better_than_nominal",
        "candidate_timing_rank",
    ]

    problem_rows = join_audit[
        join_audit["possible_alignment_problem"].fillna(False)
        & join_audit["source_filename"].notna()
        & join_audit["join_trial_count"].notna()
    ].copy()
    if problem_rows.empty:
        return pd.DataFrame(columns=columns)

    offsets_by_id = {int(pid): group.copy() for pid, group in offsets.groupby("id", dropna=True)}
    rows: list[dict[str, Any]] = []

    for problem_row_number, (_, problem) in enumerate(problem_rows.iterrows(), start=1):
        participant_id = int(problem["audit_participant_id"])
        nominal_trial = int(problem["join_trial_count"])
        participant_offsets = offsets_by_id.get(participant_id)
        if participant_offsets is None:
            continue

        nearby_offsets = participant_offsets[
            participant_offsets["trial_count"].between(
                nominal_trial - NEARBY_TRIAL_WINDOW,
                nominal_trial + NEARBY_TRIAL_WINDOW,
            )
        ].copy()
        if nearby_offsets.empty:
            continue

        trace_start = safe_float(problem.get("trace_start_onset_seconds"))
        trace_end = safe_float(problem.get("trace_end_onset_seconds"))
        fig_duration = safe_float(problem.get("fig_duration_seconds"))

        for _, candidate in nearby_offsets.iterrows():
            candidate_trial = int(candidate["trial_count"])
            expected_trace_end = trace_start + safe_float(candidate.get("end_trigger"))
            trace_difference = abs(trace_end - expected_trace_end) if np.isfinite(trace_end) and np.isfinite(expected_trace_end) else np.nan
            stimulus_difference = (
                abs(fig_duration - safe_float(candidate.get("stimulus_mt")))
                if np.isfinite(fig_duration) and np.isfinite(safe_float(candidate.get("stimulus_mt")))
                else np.nan
            )
            score = candidate_timing_score(trace_difference, stimulus_difference)
            trace_within = bool(np.isfinite(trace_difference) and trace_difference <= TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD)
            stimulus_within = bool(
                np.isfinite(stimulus_difference)
                and stimulus_difference <= STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD
            )
            rows.append(
                {
                    "problem_row_number": problem_row_number,
                    "audit_participant_id": participant_id,
                    "source_filename": problem.get("source_filename", ""),
                    "file_role": problem.get("file_role", ""),
                    "split_part": problem.get("split_part", pd.NA),
                    "raw_trial_sequence": problem.get("raw_trial_sequence", pd.NA),
                    "join_trial_count": nominal_trial,
                    "join_status": problem.get("join_status", ""),
                    "alignment_problem_codes": problem.get("alignment_problem_codes", ""),
                    "candidate_offset_trial_count": candidate_trial,
                    "candidate_offset_delta": candidate_trial - nominal_trial,
                    "candidate_trace_end_difference_seconds": trace_difference,
                    "candidate_stimulus_difference_seconds": stimulus_difference,
                    "candidate_timing_score_seconds": score,
                    "candidate_trace_within_old_threshold": trace_within,
                    "candidate_stimulus_within_old_threshold": stimulus_within,
                    "candidate_clean_by_timing": trace_within and stimulus_within,
                }
            )

    if not rows:
        return pd.DataFrame(columns=columns)

    candidates = pd.DataFrame(rows)
    nominal_scores = candidates.loc[
        candidates["candidate_offset_delta"].eq(0),
        ["problem_row_number", "candidate_timing_score_seconds"],
    ].rename(columns={"candidate_timing_score_seconds": "nominal_candidate_timing_score_seconds"})
    candidates = candidates.merge(nominal_scores, on="problem_row_number", how="left")
    candidates["candidate_better_than_nominal"] = (
        candidates["candidate_timing_score_seconds"].notna()
        & candidates["nominal_candidate_timing_score_seconds"].notna()
        & candidates["candidate_timing_score_seconds"].lt(candidates["nominal_candidate_timing_score_seconds"])
    )
    candidates["candidate_timing_rank"] = candidates.groupby("problem_row_number")[
        "candidate_timing_score_seconds"
    ].rank(method="first", na_option="bottom")
    candidates = candidates.sort_values(
        ["audit_participant_id", "source_filename", "problem_row_number", "candidate_timing_rank"],
        kind="mergesort",
    )
    return coerce_nullable_integer_columns(candidates[columns].reset_index(drop=True))


def problem_has_missing_or_extra_trigger(problem_codes: str) -> bool:
    """Return whether problem codes contain missing or extra trigger structure.

    Args:
        problem_codes: Semicolon-delimited alignment problem codes.

    Returns:
        ``True`` when any code starts with ``missing_`` or ``extra_``.

    Side effects:
        None.
    """

    return any(code.startswith("missing_") or code.startswith("extra_") for code in split_issue_codes(problem_codes))


def build_first_mismatch_by_file(join_audit: pd.DataFrame, file_summary: pd.DataFrame) -> pd.DataFrame:
    """Find the first direct alignment problem for each EDF file.

    Args:
        join_audit: Proposed full outer join audit.
        file_summary: One-row-per-EDF file facts.

    Returns:
        DataFrame with one row per EDF file describing the first direct
        alignment problem, if any, and simple cascade-after-first counts.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    sort_file_summary = file_summary.sort_values(["participant_id", "source_filename"], kind="mergesort")

    for _, file_row in sort_file_summary.iterrows():
        source_filename = file_row["source_filename"]
        file_rows = join_audit[join_audit["source_filename"].eq(source_filename)].copy()
        base = {
            "participant_id": file_row.get("participant_id", pd.NA),
            "source_filename": source_filename,
            "file_role": file_row.get("file_role", ""),
            "split_part": file_row.get("split_part", pd.NA),
            "annotation_count": file_row.get("annotation_count", np.nan),
            "mapped_task_event_count": file_row.get("mapped_task_event_count", np.nan),
        }

        if file_rows.empty:
            rows.append(
                {
                    **base,
                    "has_alignment_problem": True,
                    "first_problem_codes": "no_proposed_joinable_trial_for_file",
                    "first_raw_trial_sequence": pd.NA,
                    "first_join_trial_count": pd.NA,
                    "first_audit_trial_count": pd.NA,
                    "first_join_status": "",
                    "first_trace_difference_seconds": np.nan,
                    "first_stimulus_difference_seconds": np.nan,
                    "file_join_rows": 0,
                    "file_problem_rows": 0,
                    "rows_after_first_problem": 0,
                    "problem_rows_after_first_problem": 0,
                    "problem_fraction_after_first_problem": np.nan,
                    "first_problem_has_missing_or_extra_trigger": False,
                    "possible_cascade_after_first_problem": False,
                    "possible_cascade_after_missing_or_extra_trigger": False,
                }
            )
            continue

        sorted_rows = file_rows.sort_values(
            ["raw_trial_sequence", "audit_trial_count"],
            kind="mergesort",
            na_position="last",
        ).reset_index(drop=True)
        problem_mask = sorted_rows["possible_alignment_problem"].fillna(False)
        file_problem_rows = int(problem_mask.sum())

        if file_problem_rows == 0:
            rows.append(
                {
                    **base,
                    "has_alignment_problem": False,
                    "first_problem_codes": "",
                    "first_raw_trial_sequence": pd.NA,
                    "first_join_trial_count": pd.NA,
                    "first_audit_trial_count": pd.NA,
                    "first_join_status": "",
                    "first_trace_difference_seconds": np.nan,
                    "first_stimulus_difference_seconds": np.nan,
                    "file_join_rows": len(sorted_rows),
                    "file_problem_rows": 0,
                    "rows_after_first_problem": 0,
                    "problem_rows_after_first_problem": 0,
                    "problem_fraction_after_first_problem": 0.0,
                    "first_problem_has_missing_or_extra_trigger": False,
                    "possible_cascade_after_first_problem": False,
                    "possible_cascade_after_missing_or_extra_trigger": False,
                }
            )
            continue

        first_position = int(np.flatnonzero(problem_mask.to_numpy())[0])
        first_problem = sorted_rows.iloc[first_position]
        after_first = sorted_rows.iloc[first_position + 1 :].copy()
        after_problem_count = int(after_first["possible_alignment_problem"].fillna(False).sum())
        rows_after_first = len(after_first)
        fraction_after = after_problem_count / rows_after_first if rows_after_first else 0.0
        first_codes = text_or_empty(first_problem.get("alignment_problem_codes", ""))
        first_has_missing_or_extra = problem_has_missing_or_extra_trigger(first_codes)
        possible_cascade = after_problem_count >= 3 and fraction_after >= 0.5

        rows.append(
            {
                **base,
                "has_alignment_problem": True,
                "first_problem_codes": first_codes,
                "first_raw_trial_sequence": first_problem.get("raw_trial_sequence", pd.NA),
                "first_join_trial_count": first_problem.get("join_trial_count", pd.NA),
                "first_audit_trial_count": first_problem.get("audit_trial_count", pd.NA),
                "first_join_status": first_problem.get("join_status", ""),
                "first_trace_difference_seconds": first_problem.get("trace_epoch_end_difference_seconds", np.nan),
                "first_stimulus_difference_seconds": first_problem.get("stimulus_duration_difference_seconds", np.nan),
                "file_join_rows": len(sorted_rows),
                "file_problem_rows": file_problem_rows,
                "rows_after_first_problem": rows_after_first,
                "problem_rows_after_first_problem": after_problem_count,
                "problem_fraction_after_first_problem": fraction_after,
                "first_problem_has_missing_or_extra_trigger": first_has_missing_or_extra,
                "possible_cascade_after_first_problem": possible_cascade,
                "possible_cascade_after_missing_or_extra_trigger": first_has_missing_or_extra and possible_cascade,
            }
        )

    out = pd.DataFrame(rows)
    return coerce_nullable_integer_columns(out)


def compute_join_metrics(label: str, join_audit: pd.DataFrame, trial_summary: pd.DataFrame) -> dict[str, Any]:
    """Compute before/after summary metrics for a join-audit surface.

    Args:
        label: Metric label such as ``"before"`` or ``"after"``.
        join_audit: Full outer join audit.
        trial_summary: Trial-summary table for the same event surface.

    Returns:
        Dictionary with row counts, timing counts, and missing-event counts.

    Side effects:
        None.
    """

    join_status = join_audit["join_status"].fillna("")
    clean_timing_rows = (
        join_status.eq("raw_annotation_and_offset")
        & ~join_audit["old_trace_epoch_duration_mismatch"].fillna(False).map(is_true)
        & ~join_audit["old_stimulus_duration_mismatch"].fillna(False).map(is_true)
        & join_audit["missing_expected_event_names"].fillna("").eq("")
        & ~join_audit.get("raw_join_key_duplicate", pd.Series(False, index=join_audit.index)).fillna(False).map(is_true)
    )

    row: dict[str, Any] = {
        "surface": label,
        "matched_rows": int(join_status.eq("raw_annotation_and_offset").sum()),
        "offset_only_rows": int(join_status.eq("offset_without_raw_annotation_trial").sum()),
        "raw_only_rows": int(join_status.eq("raw_annotation_without_offset").sum()),
        "clean_timing_rows": int(clean_timing_rows.sum()),
        "trace_duration_mismatch_rows": int(join_audit["old_trace_epoch_duration_mismatch"].fillna(False).map(is_true).sum()),
        "stimulus_duration_mismatch_rows": int(
            join_audit["old_stimulus_duration_mismatch"].fillna(False).map(is_true).sum()
        ),
        "duplicate_raw_join_key_rows": int(
            join_audit.get("raw_join_key_duplicate", pd.Series(False, index=join_audit.index)).fillna(False).map(is_true).sum()
        ),
    }
    for event_name in EXPECTED_TRIAL_EVENT_NAMES:
        column = f"missing_{event_name}"
        if column in trial_summary.columns:
            row[f"{column}_rows"] = int(trial_summary[column].fillna(False).map(is_true).sum())
        else:
            row[f"{column}_rows"] = 0
    return row


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> list[str]:
    """Format dictionaries as a simple Markdown table.

    Args:
        rows: Row dictionaries.
        columns: ``(key, header)`` pairs controlling output order.

    Returns:
        Markdown table lines.

    Side effects:
        None.
    """

    header = "| " + " | ".join(title for _, title in columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows:
        values = []
        for key, _ in columns:
            value = row.get(key, "")
            if isinstance(value, float) and np.isfinite(value):
                values.append(f"{value:.3f}")
            elif pd.isna(value):
                values.append("")
            else:
                values.append(str(value))
        body.append("| " + " | ".join(values) + " |")
    return [header, divider, *body]


def format_id_list(values: pd.Series | list[Any] | set[Any]) -> str:
    """Format participant IDs as compact comma-delimited text.

    Args:
        values: Values that may contain missing data.

    Returns:
        Sorted comma-delimited integer ID text, or ``none``.

    Side effects:
        None.
    """

    ids = sorted({int(value) for value in values if pd.notna(value)})
    return ", ".join(str(value) for value in ids) if ids else "none"


def build_event_sequence_audit_summary(
    proposed_events: pd.DataFrame,
    proposed_trial_summary: pd.DataFrame,
    proposed_join_audit: pd.DataFrame,
    first_mismatch_by_file: pd.DataFrame,
    nearby_trial_candidates: pd.DataFrame,
    before_metrics: dict[str, Any],
    after_metrics: dict[str, Any],
    offsets: pd.DataFrame,
    started_at: datetime,
) -> str:
    """Build the Markdown summary for the event-sequence audit.

    Args:
        proposed_events: One-row-per-annotation table with proposed labels.
        proposed_trial_summary: Per-trial table rebuilt from proposed labels.
        proposed_join_audit: Proposed full outer raw/offset join audit.
        first_mismatch_by_file: First direct alignment problem by EDF file.
        nearby_trial_candidates: Nearby offset candidate timing audit.
        before_metrics: Metrics from script-01 raw join output.
        after_metrics: Metrics from this proposed-label join output.
        offsets: Old-compatible behavioural/tracing offsets.
        started_at: Timestamp captured at run start.

    Returns:
        Markdown summary text.

    Side effects:
        None.
    """

    finished_at = datetime.now()
    action_counts = proposed_events["proposed_sequence_action"].value_counts(dropna=False).sort_index()
    join_counts = proposed_join_audit["join_status"].value_counts(dropna=False).sort_index()
    trial_issue_counts = count_issue_codes(proposed_trial_summary["trial_issue_codes"])
    join_issue_counts = count_issue_codes(proposed_join_audit["join_audit_issue_codes"])

    files_with_problem = first_mismatch_by_file[first_mismatch_by_file["has_alignment_problem"].fillna(False).map(is_true)]
    possible_cascade_files = first_mismatch_by_file[
        first_mismatch_by_file["possible_cascade_after_missing_or_extra_trigger"].fillna(False).map(is_true)
    ]
    unresolved_events = proposed_events[proposed_events["proposed_sequence_action"].eq("unresolved")]
    split_rows = first_mismatch_by_file[first_mismatch_by_file["participant_id"].isin(SPLIT_EDF_IDS)]
    concatenated_rows = first_mismatch_by_file[first_mismatch_by_file["participant_id"].isin(CONCATENATED_EDF_IDS)]
    factual_special_rows = first_mismatch_by_file[first_mismatch_by_file["participant_id"].isin(FACTUAL_SPECIAL_CASE_IDS)]
    better_candidates = nearby_trial_candidates[
        nearby_trial_candidates["candidate_better_than_nominal"].fillna(False).map(is_true)
    ]

    metric_columns = [
        ("surface", "surface"),
        ("matched_rows", "matched"),
        ("offset_only_rows", "offset-only"),
        ("raw_only_rows", "raw-only"),
        ("clean_timing_rows", "clean timing"),
        ("trace_duration_mismatch_rows", "trace mismatch"),
        ("stimulus_duration_mismatch_rows", "stimulus mismatch"),
        ("duplicate_raw_join_key_rows", "duplicate raw keys"),
        ("missing_stim_on_rows", "missing stim_on"),
        ("missing_red_on_rows", "missing red_on"),
        ("missing_trace_start_rows", "missing trace_start"),
        ("missing_trace_end_rows", "missing trace_end"),
    ]

    lines = [
        "# Raw Event Sequence Audit Summary",
        "",
        f"Generated: {finished_at.isoformat(timespec='seconds')}",
        f"Started: {started_at.isoformat(timespec='seconds')}",
        f"Duration seconds: {(finished_at - started_at).total_seconds():.1f}",
        "",
        "## Scope",
        "",
        "This local audit ports historical duplicate/mislabeled-trigger detection into an in-memory proposed annotation table. It reads script-01 annotation CSVs and old-compatible offsets; it does not open EDFs, write repaired annotations, preprocess EEG, filter EEG, epoch EEG, or make participant-level decisions.",
        "",
        "Historical cleanup rules evaluated:",
        "",
        f"- consecutive same-label trigger with next gap < {HISTORICAL_DUPLICATE_TRIGGER_SECONDS * 1000:.0f} ms: mark as duplicated and remove from the proposed sequence;",
        "- consecutive same-label trigger with a larger next gap: relabel to the next name in the old event order;",
        f"- practice-like groups: `red_on - stim_on > {PRACTICE_FIG_DURATION_SECONDS_THRESHOLD * 1000:.0f} ms`; bad-start-like groups: `red_on` without immediately preceding `stim_on`;",
        f"- timing checks: `trace_end` versus `trace_start + end_trigger` > {TRACE_EPOCH_MISMATCH_SECONDS_THRESHOLD * 1000:.0f} ms, and `red_on - stim_on` versus `stimulus_mt` > {STIMULUS_DURATION_MISMATCH_SECONDS_THRESHOLD * 1000:.0f} ms.",
        "",
        "## Inputs",
        "",
        f"- First-audit annotation events: `{(EVENT_ALIGNMENT_DIR / RAW_ANNOTATION_EVENTS_FILENAME).as_posix()}`",
        f"- First-audit trial summary: `{(EVENT_ALIGNMENT_DIR / ANNOTATION_TRIAL_SUMMARY_FILENAME).as_posix()}`",
        f"- First-audit join audit: `{(EVENT_ALIGNMENT_DIR / OFFSET_ANNOTATION_JOIN_AUDIT_FILENAME).as_posix()}`",
        f"- Old-compatible offsets: `{OLD_COMPATIBLE_OFFSETS_PATH.as_posix()}` ({len(offsets):,} rows)",
        "",
        "## Outputs",
        "",
        f"- `{(EVENT_SEQUENCE_AUDIT_DIR / PROPOSED_ANNOTATION_EVENTS_FILENAME).as_posix()}`",
        f"- `{(EVENT_SEQUENCE_AUDIT_DIR / PROPOSED_TRIAL_SUMMARY_FILENAME).as_posix()}`",
        f"- `{(EVENT_SEQUENCE_AUDIT_DIR / PROPOSED_OFFSET_JOIN_AUDIT_FILENAME).as_posix()}`",
        f"- `{(EVENT_SEQUENCE_AUDIT_DIR / FIRST_MISMATCH_BY_FILE_FILENAME).as_posix()}`",
        f"- `{(EVENT_SEQUENCE_AUDIT_DIR / NEARBY_TRIAL_CANDIDATE_AUDIT_FILENAME).as_posix()}`",
        "",
        "## Before vs After",
        "",
        *markdown_table([before_metrics, after_metrics], metric_columns),
        "",
        "`before` is the first raw annotation join from script 01. `after` is the proposed-label sequence produced by this audit.",
        "",
        "## Proposed Event Actions",
        "",
    ]

    for action, count in action_counts.items():
        lines.append(f"- {action}: {int(count):,}")

    lines.extend(
        [
            "",
            f"- Unresolved proposed-label rows: {len(unresolved_events):,} "
            f"({format_id_list(unresolved_events['participant_id']) if not unresolved_events.empty else 'none'})",
            "",
            "## Proposed Join Status",
            "",
        ]
    )

    for status, count in join_counts.items():
        lines.append(f"- {status}: {int(count):,}")

    lines.extend(["", "Trial issue counts:", ""])
    if trial_issue_counts:
        for issue, count in sorted(trial_issue_counts.items()):
            lines.append(f"- {issue}: {count:,}")
    else:
        lines.append("- none")

    lines.extend(["", "Join issue counts:", ""])
    if join_issue_counts:
        for issue, count in sorted(join_issue_counts.items()):
            lines.append(f"- {issue}: {count:,}")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Timing Candidate Checks",
            "",
            f"- Nearby candidate rows written: {len(nearby_trial_candidates):,}",
            f"- Candidate rows with lower timing score than the nominal join: {len(better_candidates):,}",
            f"- Files with a direct alignment problem row: {len(files_with_problem):,} "
            f"({format_id_list(files_with_problem['participant_id']) if not files_with_problem.empty else 'none'})",
            f"- Files with a possible cascade after a missing/extra trigger: {len(possible_cascade_files):,} "
            f"({format_id_list(possible_cascade_files['participant_id']) if not possible_cascade_files.empty else 'none'})",
            "",
            "The nearby-candidate table reports offset trials within "
            f"+/- {NEARBY_TRIAL_WINDOW} trial counts of each problematic nominal join and ranks candidates by available trace/stimulus timing differences.",
            "",
            "## Split, Concatenated, and Factual Special Cases",
            "",
            f"- Split EDF IDs reported separately: {format_id_list(split_rows['participant_id']) if not split_rows.empty else 'none'}",
            f"- Concatenated EDF ID reported separately: {format_id_list(concatenated_rows['participant_id']) if not concatenated_rows.empty else 'none'}",
            f"- Zero/low/no-offset factual special cases present in file report: {format_id_list(factual_special_rows['participant_id']) if not factual_special_rows.empty else 'none'}",
            "",
            "No split file is concatenated, and ID 5 is not cropped or split by this audit.",
            "",
            "## Safety Boundary",
            "",
            "- Proposed labels are audit fields only.",
            "- No raw EDFs are read or modified.",
            "- No signal data is loaded.",
            "- No preprocessing, filtering, or epoch construction is performed.",
            "- Tony makes scientific/event-handling decisions after reviewing the audit evidence.",
            "",
        ]
    )

    return "\n".join(lines)


def write_outputs(
    proposed_events: pd.DataFrame,
    proposed_trial_summary: pd.DataFrame,
    proposed_join_audit: pd.DataFrame,
    first_mismatch_by_file: pd.DataFrame,
    nearby_trial_candidates: pd.DataFrame,
    summary_text: str,
    output_dir: Path,
) -> None:
    """Write local event-sequence audit outputs.

    Args:
        proposed_events: One-row-per-annotation proposed-label table.
        proposed_trial_summary: Proposed-label trial summary.
        proposed_join_audit: Proposed-label full outer join audit.
        first_mismatch_by_file: First direct alignment problem by file.
        nearby_trial_candidates: Nearby offset timing candidates.
        summary_text: Markdown summary text.
        output_dir: Local output directory below ``_Data/eeg``.

    Returns:
        None.

    Side effects:
        Creates ``output_dir`` if needed and writes CSV/Markdown audit outputs.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    proposed_events.to_csv(output_dir / PROPOSED_ANNOTATION_EVENTS_FILENAME, index=False)
    proposed_trial_summary.to_csv(output_dir / PROPOSED_TRIAL_SUMMARY_FILENAME, index=False)
    proposed_join_audit.to_csv(output_dir / PROPOSED_OFFSET_JOIN_AUDIT_FILENAME, index=False)
    first_mismatch_by_file.to_csv(output_dir / FIRST_MISMATCH_BY_FILE_FILENAME, index=False)
    nearby_trial_candidates.to_csv(output_dir / NEARBY_TRIAL_CANDIDATE_AUDIT_FILENAME, index=False)
    (output_dir / EVENT_SEQUENCE_AUDIT_SUMMARY_FILENAME).write_text(summary_text, encoding="utf-8")


def run_event_sequence_audit(
    repo_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """Run the event-sequence audit in memory.

    Args:
        repo_root: Repository root.

    Returns:
        ``(proposed_events, proposed_trial_summary, proposed_join_audit,
        first_mismatch_by_file, nearby_trial_candidates, summary_text)``.

    Side effects:
        Reads local CSV/Markdown inputs. Does not write files.
    """

    started_at = datetime.now()
    raw_events, first_trial_summary, first_join_audit, offsets, _ = read_audit_inputs(repo_root)
    file_summary = build_file_summary_from_first_trial_summary(first_trial_summary)

    proposed_events = propose_annotation_event_sequence(raw_events)
    proposed_events = assign_proposed_trial_numbers(proposed_events)
    proposed_trial_summary = build_proposed_trial_summary(proposed_events, file_summary)
    proposed_join_audit = build_proposed_offset_join_audit(proposed_trial_summary, offsets, file_summary)
    first_mismatch = build_first_mismatch_by_file(proposed_join_audit, file_summary)
    nearby_candidates = build_nearby_trial_candidate_audit(proposed_join_audit, offsets)

    before_metrics = compute_join_metrics("before", first_join_audit, first_trial_summary)
    after_metrics = compute_join_metrics("after", proposed_join_audit, proposed_trial_summary)
    summary_text = build_event_sequence_audit_summary(
        proposed_events=proposed_events,
        proposed_trial_summary=proposed_trial_summary,
        proposed_join_audit=proposed_join_audit,
        first_mismatch_by_file=first_mismatch,
        nearby_trial_candidates=nearby_candidates,
        before_metrics=before_metrics,
        after_metrics=after_metrics,
        offsets=offsets,
        started_at=started_at,
    )

    return (
        proposed_events,
        proposed_trial_summary,
        proposed_join_audit,
        first_mismatch,
        nearby_candidates,
        summary_text,
    )


def main() -> None:
    """Run the raw event-sequence audit and write local outputs.

    Args:
        None. Paths are resolved from this script's repository location.

    Returns:
        None.

    Side effects:
        Reads script-01 CSV outputs and old-compatible offsets, creates
        ``_Data/eeg/event_sequence_audit/`` if needed, and writes local-only
        audit outputs.
    """

    repo_root = repo_root_from_script()
    output_dir = repo_root / EVENT_SEQUENCE_AUDIT_DIR
    try:
        (
            proposed_events,
            proposed_trial_summary,
            proposed_join_audit,
            first_mismatch,
            nearby_candidates,
            summary_text,
        ) = run_event_sequence_audit(repo_root)
        write_outputs(
            proposed_events=proposed_events,
            proposed_trial_summary=proposed_trial_summary,
            proposed_join_audit=proposed_join_audit,
            first_mismatch_by_file=first_mismatch,
            nearby_trial_candidates=nearby_candidates,
            summary_text=summary_text,
            output_dir=output_dir,
        )
    except Exception as error:
        raise SystemExit(f"FAIL raw event sequence audit: {error}") from error

    print(f"Wrote raw event-sequence audit outputs to {output_dir}")
    print(f"Proposed annotation rows: {len(proposed_events):,}")
    print(f"Proposed join rows: {len(proposed_join_audit):,}")
    print(f"Nearby timing candidate rows: {len(nearby_candidates):,}")


if __name__ == "__main__":
    main()
