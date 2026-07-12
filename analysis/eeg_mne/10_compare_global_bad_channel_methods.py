"""Compare automated continuous global bad-channel methods for DEMI EEG.

This audit script evaluates MNE-native LOF, PyPREP ``NoisyChannels``, and the
integrated PyPREP/PREP pipeline on the same explicitly prepared full-session
EEG input. The prepared input uses EEG channels only, the approved
``standard_1005`` montage, a fixed audit-only 250-Hz resampling rate, and no
notch, reference, or analysis filter. Every method receives an independent
in-memory copy.

The public YAML contract owns all method parameters. A private YAML file owns
the participant/file selection, current descriptive-candidate context,
historical-note context, and selected repeated full-PREP runs. Generated
channel-level and participant-level conclusions remain local under ``_Data``.

Outputs
-------
* ``global_bad_method_runs.csv``: one row per file/method/variant/repeat;
* ``global_bad_channel_details.csv``: channel flags and available scores;
* ``global_bad_method_agreement.csv``: pairwise primary-method agreement;
* ``global_bad_method_determinism.csv``: repeated-run equality checks;
* ``global_bad_method_exceptions.csv``: objective algorithm/review exceptions;
* ``global_bad_method_manifest.json`` and generated Markdown summary.

This is not production preprocessing. It does not modify EDFs, save filtered or
referenced signals, accept a production method, construct epochs, run ICA,
compute CSD, or make participant-level inclusion decisions. Full PREP signal
changes are discarded after their metadata are summarized.

Run from the repository root with::

    ./.venv/bin/python analysis/eeg_mne/10_compare_global_bad_channel_methods.py
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import autoreject
import mne
import numpy as np
import pandas as pd
import pyprep
import sklearn
import yaml

from preprocessing_parameter_audit import (
    CRITERION_KEYS,
    channel_list_json,
    load_audit_contract,
    prepare_continuous_eeg,
    run_full_prep,
    run_lof,
    run_noisy_channels,
    validate_audit_output_path,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = REPO_ROOT / "analysis/eeg_mne/preprocessing_parameter_audit_contract_v1.yaml"
DEFAULT_PRIVATE_CONFIG = REPO_ROOT / "_Private/config/preprocessing_parameter_evidence_pass_v1.yaml"


def load_private_config(path: Path) -> dict[str, Any]:
    """Load and validate the local participant/file selection.

    The file is private because it contains participant-level evidence routing.
    This function reads YAML only and writes nothing.
    """

    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or not config.get("recordings"):
        raise ValueError("Private comparison config must list recordings.")
    names = [row["source_filename"] for row in config["recordings"]]
    if len(names) != len(set(names)):
        raise ValueError("Private comparison config contains duplicate filenames.")
    return config


@contextlib.contextmanager
def quiet_third_party_output() -> Any:
    """Suppress verbose PyPREP/MNE progress output while preserving exceptions."""

    with open(os.devnull, "w", encoding="utf-8") as null:
        with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
            yield


def execute_method(call: Callable[[], dict[str, Any]]) -> tuple[dict[str, Any] | None, str | None, float]:
    """Run one in-memory method, returning result/error and elapsed seconds."""

    started = time.perf_counter()
    try:
        with quiet_third_party_output():
            result = call()
        return result, None, time.perf_counter() - started
    except Exception as exc:  # Each file/method failure must remain auditable.
        return None, f"{type(exc).__name__}: {exc}", time.perf_counter() - started


def candidate_contexts(config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read local descriptive candidates and historical note context.

    Returns raw private/local evidence tables. No candidate is converted into a
    bad-channel label.
    """

    candidates = pd.read_csv(REPO_ROOT / config["channel_qc_candidates"])
    historical = pd.read_csv(REPO_ROOT / config["historical_qc_queue"])
    return candidates, historical


def channel_context(
    source_filename: str,
    channel: str,
    candidates: pd.DataFrame,
    historical: pd.DataFrame,
) -> dict[str, Any]:
    """Return factual candidate and historical-note mention fields for one channel."""

    rows = candidates[
        (candidates["source_filename"] == source_filename)
        & (candidates["normalized_channel_name"] == channel)
    ]
    eeg_rows = rows[rows["configured_channel_type"] == "eeg"]
    historical_rows = historical[historical["raw_file"] == source_filename]
    note = " ".join(historical_rows["historical_masterlist_notes"].fillna("").astype(str))
    return {
        "current_eeg_descriptive_candidate": not eeg_rows.empty,
        "current_candidate_reasons": ";".join(eeg_rows["candidate_reasons"].astype(str)),
        "historical_masterlist_channel_mention": channel.upper() in note.upper(),
        "historical_prep_log_available": False,
    }


def base_run_row(filename: str, groups: list[str], raw: mne.io.BaseRaw) -> dict[str, Any]:
    """Return immutable prepared-input metadata shared by all method rows."""

    return {
        "source_filename": filename,
        "selection_groups_json": json.dumps(sorted(groups), separators=(",", ":")),
        "prepared_n_channels": len(raw.ch_names),
        "prepared_n_times": raw.n_times,
        "prepared_duration_seconds": float(raw.times[-1]),
        "prepared_sfreq_hz": float(raw.info["sfreq"]),
        "prepared_montage": "standard_1005",
        "prepared_eeg_only": True,
        "prepared_reference_applied": False,
        "prepared_notch_applied": False,
        "prepared_analysis_filter_applied": False,
    }


def add_lof_rows(
    filename: str,
    groups: list[str],
    raw: mne.io.BaseRaw,
    contract: dict[str, Any],
    candidates: pd.DataFrame,
    historical: pd.DataFrame,
    runs: list[dict[str, Any]],
    details: list[dict[str, Any]],
) -> None:
    """Run primary/repeated and documented LOF sensitivity settings."""

    settings = contract["global_bad_channel_methods"]["lof"]
    variants = [("primary", settings["primary"], int(settings["repeats"]))]
    variants.extend((row["label"], row, 1) for row in settings["sensitivity"])
    for label, params, repeats in variants:
        for repeat in range(1, repeats + 1):
            result, error, elapsed = execute_method(
                lambda p=params: run_lof(
                    raw,
                    n_neighbors=int(p["n_neighbors"]),
                    metric=str(p["metric"]),
                    threshold=float(p["threshold"]),
                )
            )
            row = {
                **base_run_row(filename, groups, raw),
                "method": "mne_lof",
                "variant": label,
                "repeat": repeat,
                "success": error is None,
                "error": error,
                "elapsed_seconds": elapsed,
                "parameters_json": json.dumps(params, sort_keys=True),
                "bad_all_json": channel_list_json(result["bad_all"]) if result else "[]",
                "bad_count": len(result["bad_all"]) if result else None,
            }
            runs.append(row)
            if result:
                for channel in raw.ch_names:
                    details.append(
                        {
                            "source_filename": filename,
                            "method": "mne_lof",
                            "variant": label,
                            "repeat": repeat,
                            "channel": channel,
                            "flagged": channel in result["bad_all"],
                            "lof_score": result["scores"][channel],
                            **channel_context(filename, channel, candidates, historical),
                        }
                    )


def add_noisy_rows(
    filename: str,
    groups: list[str],
    raw: mne.io.BaseRaw,
    contract: dict[str, Any],
    candidates: pd.DataFrame,
    historical: pd.DataFrame,
    runs: list[dict[str, Any]],
    details: list[dict[str, Any]],
) -> None:
    """Run repeated full-session PyPREP NoisyChannels and retain each criterion."""

    method_root = contract["global_bad_channel_methods"]
    settings = method_root["noisy_channels"]
    seed = int(method_root["random_seed"])
    for repeat in range(1, int(settings["repeats"]) + 1):
        result, error, elapsed = execute_method(lambda: run_noisy_channels(raw, settings, seed))
        row = {
            **base_run_row(filename, groups, raw),
            "method": "pyprep_noisy_channels",
            "variant": "full_criteria",
            "repeat": repeat,
            "success": error is None,
            "error": error,
            "elapsed_seconds": elapsed,
            "parameters_json": json.dumps(settings, sort_keys=True),
            "bad_all_json": channel_list_json(result["bad_all"]) if result else "[]",
            "bad_count": len(result["bad_all"]) if result else None,
        }
        if result:
            for key in CRITERION_KEYS:
                row[f"{key}_json"] = channel_list_json(result["criteria"][key])
        runs.append(row)
        if result:
            for channel_row in result["channel_details"]:
                details.append(
                    {
                        "source_filename": filename,
                        "method": "pyprep_noisy_channels",
                        "variant": "full_criteria",
                        "repeat": repeat,
                        "flagged": channel_row["channel"] in result["bad_all"],
                        **channel_row,
                        **channel_context(
                            filename, channel_row["channel"], candidates, historical
                        ),
                    }
                )


def add_full_prep_rows(
    filename: str,
    groups: list[str],
    raw: mne.io.BaseRaw,
    contract: dict[str, Any],
    repeat_files: set[str],
    candidates: pd.DataFrame,
    historical: pd.DataFrame,
    runs: list[dict[str, Any]],
    details: list[dict[str, Any]],
) -> None:
    """Run integrated PREP once, or twice for the private determinism subset."""

    method_root = contract["global_bad_channel_methods"]
    settings = method_root["full_prep"]
    seed = int(method_root["random_seed"])
    repeats = 2 if filename in repeat_files else int(settings["repeats_default"])
    for repeat in range(1, repeats + 1):
        result, error, elapsed = execute_method(lambda: run_full_prep(raw, settings, seed))
        row = {
            **base_run_row(filename, groups, raw),
            "method": "pyprep_full_prep",
            "variant": "integrated",
            "repeat": repeat,
            "success": error is None,
            "error": error,
            "elapsed_seconds": elapsed,
            "parameters_json": json.dumps(settings, sort_keys=True),
            "bad_all_json": (
                channel_list_json(result["bad_before_interpolation"]) if result else "[]"
            ),
            "bad_count": len(result["bad_before_interpolation"]) if result else None,
        }
        if result:
            row.update(
                {
                    "original_bad_json": channel_list_json(result["original_bad_all"]),
                    "unusable_channels_json": channel_list_json(result["unusable_channels"]),
                    "bad_before_interpolation_json": channel_list_json(
                        result["bad_before_interpolation"]
                    ),
                    "interpolated_channels_json": channel_list_json(
                        result["interpolated_channels"]
                    ),
                    "still_noisy_channels_json": channel_list_json(
                        result["still_noisy_channels"]
                    ),
                    "interpolated_count": result["interpolated_count"],
                    "interpolated_proportion": result["interpolated_proportion"],
                    "reference_kind": result["reference_kind"],
                    "reference_channels_json": channel_list_json(result["reference_channels"]),
                    "rereference_channels_json": channel_list_json(
                        result["rereference_channels"]
                    ),
                    "reference_before_hash": result["reference_before_hash"],
                    "reference_after_hash": result["reference_after_hash"],
                    "reference_before_rms_v": result["reference_before_rms_v"],
                    "reference_after_rms_v": result["reference_after_rms_v"],
                    "original_criteria_json": json.dumps(
                        result["original_criteria"], sort_keys=True
                    ),
                    "before_interpolation_criteria_json": json.dumps(
                        result["before_interpolation_criteria"], sort_keys=True
                    ),
                    "after_interpolation_criteria_json": json.dumps(
                        result["after_interpolation_criteria"], sort_keys=True
                    ),
                }
            )
        runs.append(row)
        if result:
            for channel in raw.ch_names:
                details.append(
                    {
                        "source_filename": filename,
                        "method": "pyprep_full_prep",
                        "variant": "integrated",
                        "repeat": repeat,
                        "channel": channel,
                        "flagged": channel in result["bad_before_interpolation"],
                        "original_bad": channel in result["original_bad_all"],
                        "unusable": channel in result["unusable_channels"],
                        "bad_before_interpolation": channel
                        in result["bad_before_interpolation"],
                        "interpolated": channel in result["interpolated_channels"],
                        "still_noisy": channel in result["still_noisy_channels"],
                        **channel_context(filename, channel, candidates, historical),
                    }
                )


def jaccard(left: set[str], right: set[str]) -> float:
    """Return set Jaccard agreement, defining two empty sets as perfect agreement."""

    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def build_agreement(runs: pd.DataFrame) -> pd.DataFrame:
    """Compare repeat-one primary outputs without treating any method as truth."""

    primary = runs[
        (runs["repeat"] == 1)
        & (runs["success"] == True)  # noqa: E712 - explicit pandas scalar comparison.
        & (
            ((runs["method"] == "mne_lof") & (runs["variant"] == "primary"))
            | (runs["method"] == "pyprep_noisy_channels")
            | (runs["method"] == "pyprep_full_prep")
        )
    ]
    rows: list[dict[str, Any]] = []
    for filename, group in primary.groupby("source_filename"):
        sets = {
            row.method: set(json.loads(row.bad_all_json))
            for row in group.itertuples(index=False)
        }
        methods = sorted(sets)
        for index, left_name in enumerate(methods):
            for right_name in methods[index + 1 :]:
                left, right = sets[left_name], sets[right_name]
                rows.append(
                    {
                        "source_filename": filename,
                        "method_left": left_name,
                        "method_right": right_name,
                        "intersection_json": channel_list_json(left & right),
                        "left_only_json": channel_list_json(left - right),
                        "right_only_json": channel_list_json(right - left),
                        "union_json": channel_list_json(left | right),
                        "jaccard": jaccard(left, right),
                        "historical_prep_not_treated_as_ground_truth": True,
                    }
                )
    return pd.DataFrame(rows)


def build_determinism(runs: pd.DataFrame) -> pd.DataFrame:
    """Summarize whether repeated channel lists and PREP reference hashes match."""

    rows: list[dict[str, Any]] = []
    repeated = runs.groupby(["source_filename", "method", "variant"], dropna=False)
    for keys, group in repeated:
        if len(group) < 2:
            continue
        successful = group[group["success"] == True]  # noqa: E712
        bad_identical = successful["bad_all_json"].nunique(dropna=False) == 1
        hashes = successful.get("reference_after_hash", pd.Series(dtype=object)).dropna()
        reference_identical = None if hashes.empty else hashes.nunique() == 1
        rows.append(
            {
                "source_filename": keys[0],
                "method": keys[1],
                "variant": keys[2],
                "requested_repeats": len(group),
                "successful_repeats": len(successful),
                "bad_channel_lists_identical": bad_identical,
                "reference_after_hashes_identical": reference_identical,
                "deterministic": len(successful) == len(group)
                and bad_identical
                and reference_identical is not False,
            }
        )
    return pd.DataFrame(rows)


def build_exceptions(
    runs: pd.DataFrame,
    determinism: pd.DataFrame,
    warning_proportion: float,
) -> pd.DataFrame:
    """Apply objective exceptional-review criteria, not visual adjudication."""

    rows: list[dict[str, Any]] = []
    for row in runs.itertuples(index=False):
        if not row.success:
            rows.append(
                {
                    "source_filename": row.source_filename,
                    "criterion": "algorithm_error",
                    "method": row.method,
                    "detail": row.error,
                    "requires_exceptional_human_review": True,
                }
            )
        if row.method == "pyprep_full_prep" and row.success:
            still = json.loads(getattr(row, "still_noisy_channels_json", "[]") or "[]")
            proportion = getattr(row, "interpolated_proportion", np.nan)
            if still:
                rows.append(
                    {
                        "source_filename": row.source_filename,
                        "criterion": "still_noisy_after_interpolation",
                        "method": row.method,
                        "detail": channel_list_json(still),
                        "requires_exceptional_human_review": True,
                    }
                )
            if pd.notna(proportion) and float(proportion) > warning_proportion:
                rows.append(
                    {
                        "source_filename": row.source_filename,
                        "criterion": "interpolated_proportion_above_historical_warning",
                        "method": row.method,
                        "detail": f"{float(proportion):.6f} > {warning_proportion:.6f}",
                        "requires_exceptional_human_review": True,
                    }
                )
    if not determinism.empty:
        for row in determinism[determinism["deterministic"] == False].itertuples(index=False):  # noqa: E712
            rows.append(
                {
                    "source_filename": row.source_filename,
                    "criterion": "repeated_run_nondeterminism",
                    "method": row.method,
                    "detail": row.variant,
                    "requires_exceptional_human_review": True,
                }
            )
    columns = (
        "source_filename",
        "criterion",
        "method",
        "detail",
        "requires_exceptional_human_review",
    )
    return pd.DataFrame(rows, columns=columns).drop_duplicates()


def write_outputs(
    output_dir: Path,
    contract: dict[str, Any],
    config: dict[str, Any],
    runs: list[dict[str, Any]],
    details: list[dict[str, Any]],
) -> None:
    """Write method evidence and a concise generated local summary."""

    suffixes = contract["safety"]["allowed_output_suffixes"]
    paths = {
        "runs": output_dir / "global_bad_method_runs.csv",
        "details": output_dir / "global_bad_channel_details.csv",
        "agreement": output_dir / "global_bad_method_agreement.csv",
        "determinism": output_dir / "global_bad_method_determinism.csv",
        "exceptions": output_dir / "global_bad_method_exceptions.csv",
        "manifest": output_dir / "global_bad_method_manifest.json",
        "summary": output_dir / "global_bad_method_summary.md",
    }
    for path in paths.values():
        validate_audit_output_path(path, output_dir, suffixes)
    run_table = pd.DataFrame(runs)
    detail_table = pd.DataFrame(details)
    agreement = build_agreement(run_table)
    determinism = build_determinism(run_table)
    warning = float(
        contract["global_bad_channel_methods"]["full_prep"][
            "historical_interpolation_warning_proportion"
        ]
    )
    exceptions = build_exceptions(run_table, determinism, warning)
    run_table.to_csv(paths["runs"], index=False)
    detail_table.to_csv(paths["details"], index=False)
    agreement.to_csv(paths["agreement"], index=False)
    determinism.to_csv(paths["determinism"], index=False)
    exceptions.to_csv(paths["exceptions"], index=False)

    successful = run_table[run_table["success"] == True]  # noqa: E712
    primary = successful[
        ((successful["method"] == "mne_lof") & (successful["variant"] == "primary"))
        | (successful["method"] != "mne_lof")
    ]
    repeat_one = primary[primary["repeat"] == 1]
    counts = repeat_one.groupby("method")["bad_count"].agg(["count", "sum", "median"])
    genuine_files = sorted(exceptions["source_filename"].unique()) if not exceptions.empty else []
    if counts.empty:
        count_text = "No successful method runs."
    else:
        count_lines = ["| method | files | total flagged | median flagged |", "|---|---:|---:|---:|"]
        for method, row in counts.iterrows():
            count_lines.append(
                f"| {method} | {int(row['count'])} | {int(row['sum'])} | {float(row['median']):.3f} |"
            )
        count_text = "\n".join(count_lines)
    summary = f"""# Automated global bad-channel comparison

Generated: {datetime.now(timezone.utc).isoformat()}

Prepared surface: full continuous EEG, EEG-only, standard_1005, {contract['prepared_continuous_input']['resample_hz']} Hz, no reference, no notch, no analysis filter.

Files evaluated: {run_table['source_filename'].nunique()}. Current EEG descriptive candidate rows covered: {int(detail_table['current_eeg_descriptive_candidate'].sum() > 0) if detail_table.empty else detail_table[detail_table['current_eeg_descriptive_candidate']].drop_duplicates(['source_filename','channel']).shape[0]}.

Primary repeat-one automated counts by method:

{count_text}

Repeated-run checks: {len(determinism)}; deterministic: {int(determinism['deterministic'].sum()) if not determinism.empty else 0}.

Objective exceptional-review files ({len(genuine_files)}): {', '.join(genuine_files) if genuine_files else 'none'}.

Method disagreement is reported descriptively and is not itself treated as ground truth or an automatic review requirement. Auxiliary/stim hard exceptions were protected by EEG-only detector input. Historical masterlist channel mentions are context only; no historical prep_info.csv ground truth was available.

AutoReject RANSAC was not run because its public interface requires Epochs; this task constructs no epochs and does not duplicate PREP RANSAC without a concrete validation purpose.
"""
    paths["summary"].write_text(summary, encoding="utf-8")
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": contract["contract_version"],
        "versions": {
            "python": platform.python_version(),
            "mne": mne.__version__,
            "pyprep": pyprep.__version__,
            "autoreject": autoreject.__version__,
            "scikit_learn": sklearn.__version__,
            "numpy": np.__version__,
        },
        "recordings": config["recordings"],
        "signal_derivatives_written": 0,
        "epochs_constructed": 0,
        "authoritative_signal_changed": False,
        "historical_prep_ground_truth_assumed": False,
        "outputs": [str(path.relative_to(REPO_ROOT)) for path in paths.values()],
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def expected_run_keys(filename: str, contract: dict[str, Any], repeat_files: set[str]) -> set[tuple[str, str, int]]:
    """Return the exact method/variant/repeat keys required for one recording."""

    lof = contract["global_bad_channel_methods"]["lof"]
    keys = {
        ("mne_lof", "primary", repeat)
        for repeat in range(1, int(lof["repeats"]) + 1)
    }
    keys.update(("mne_lof", row["label"], 1) for row in lof["sensitivity"])
    noisy_repeats = int(
        contract["global_bad_channel_methods"]["noisy_channels"]["repeats"]
    )
    keys.update(
        ("pyprep_noisy_channels", "full_criteria", repeat)
        for repeat in range(1, noisy_repeats + 1)
    )
    full_repeats = 2 if filename in repeat_files else int(
        contract["global_bad_channel_methods"]["full_prep"]["repeats_default"]
    )
    keys.update(
        ("pyprep_full_prep", "integrated", repeat)
        for repeat in range(1, full_repeats + 1)
    )
    return keys


def load_checkpoint(output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Load prior local CSV checkpoints when a reporting/process interruption occurred."""

    run_path = output_dir / "global_bad_method_runs.csv"
    detail_path = output_dir / "global_bad_channel_details.csv"
    if not run_path.is_file() or not detail_path.is_file():
        return [], []
    return (
        pd.read_csv(run_path).replace({np.nan: None}).to_dict("records"),
        pd.read_csv(detail_path).replace({np.nan: None}).to_dict("records"),
    )


def main() -> None:
    """Run all three automated methods over the private evidence surface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--private-config", type=Path, default=DEFAULT_PRIVATE_CONFIG)
    parser.add_argument("--only-file", action="append", default=[])
    args = parser.parse_args()
    contract = load_audit_contract(args.contract)
    config = load_private_config(args.private_config)
    output_dir = REPO_ROOT / config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    candidates, historical = candidate_contexts(config)
    repeat_files = set(config["full_prep_repeat_recordings"])
    recordings = config["recordings"]
    if args.only_file:
        recordings = [row for row in recordings if row["source_filename"] in args.only_file]
        if not recordings:
            raise ValueError("--only-file did not match the private selection.")

    runs, details = load_checkpoint(output_dir)
    if runs:
        print(f"Loaded checkpoint with {len(runs)} completed method runs.", flush=True)
    prepared = contract["prepared_continuous_input"]
    for index, recording in enumerate(recordings, start=1):
        filename = recording["source_filename"]
        groups = list(recording["selection_groups"])
        completed = {
            (str(row["method"]), str(row["variant"]), int(row["repeat"]))
            for row in runs
            if row["source_filename"] == filename
        }
        required = expected_run_keys(filename, contract, repeat_files)
        if required.issubset(completed):
            print(f"[{index}/{len(recordings)}] Reusing checkpoint for {filename}", flush=True)
            continue
        print(f"[{index}/{len(recordings)}] Preparing {filename}", flush=True)
        raw = prepare_continuous_eeg(
            REPO_ROOT / config["raw_root"] / filename,
            channel_types=prepared["channel_types"],
            resample_hz=float(prepared["resample_hz"]),
        )
        print(f"[{index}/{len(recordings)}] LOF {filename}", flush=True)
        add_lof_rows(filename, groups, raw, contract, candidates, historical, runs, details)
        print(f"[{index}/{len(recordings)}] NoisyChannels {filename}", flush=True)
        add_noisy_rows(filename, groups, raw, contract, candidates, historical, runs, details)
        print(f"[{index}/{len(recordings)}] Full PREP {filename}", flush=True)
        add_full_prep_rows(
            filename,
            groups,
            raw,
            contract,
            repeat_files,
            candidates,
            historical,
            runs,
            details,
        )
        write_outputs(output_dir, contract, {**config, "recordings": recordings}, runs, details)
        print(f"[{index}/{len(recordings)}] Complete {filename}", flush=True)
    # Also refresh summaries when every requested file was satisfied by a
    # checkpoint and the loop therefore performed no new method call.
    write_outputs(output_dir, contract, {**config, "recordings": recordings}, runs, details)
    print(f"Global bad-channel comparison complete: {output_dir}")


if __name__ == "__main__":
    main()
