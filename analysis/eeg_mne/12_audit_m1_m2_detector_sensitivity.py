"""Audit the systematic M1/M2 effect on continuous global detectors.

The primary 32-EEG-channel comparison found that PyPREP marked M1 and/or M2 in
nearly every selected recording. Acquisition evidence proves that both are
recorded signals but does not establish their online-reference role. This
narrow audit therefore repeats only MNE LOF and PyPREP ``NoisyChannels`` after
removing M1/M2 from the detector pool, then compares the remaining scalp calls
with the completed 32-channel primary runs.

This is a technical-cause sensitivity supporting the later approved 30-channel
detector/reference-source policy, not a production preprocessing implementation.
The public contract fixes the excluded labels and methods; the private config
owns file selection. Full PREP is not repeated because doing so would introduce
a second integrated reference/interpolation pipeline rather than isolate the
detector-pool effect.

Outputs are local CSV/JSON/Markdown evidence. The script writes no signal,
changes no EDF, constructs no epochs, and makes no reference-policy decision.

Run from the repository root with::

    ./.venv/bin/python analysis/eeg_mne/12_audit_m1_m2_detector_sensitivity.py
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import mne
import numpy as np
import pandas as pd
import pyprep
import yaml

from preprocessing_parameter_audit import (
    channel_list_json,
    load_audit_contract,
    prepare_continuous_eeg,
    run_lof,
    run_noisy_channels,
    validate_audit_output_path,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = REPO_ROOT / "analysis/eeg_mne/preprocessing_parameter_audit_contract_v1.yaml"
DEFAULT_PRIVATE_CONFIG = REPO_ROOT / "_Private/config/preprocessing_parameter_evidence_pass_v1.yaml"


def load_private_config(path: Path) -> dict[str, Any]:
    """Load the private file surface and local output paths."""

    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not config.get("recordings"):
        raise ValueError("Private sensitivity config has no recordings.")
    return config


@contextlib.contextmanager
def quiet_third_party_output() -> Any:
    """Suppress progress bars while leaving exceptions visible to the caller."""

    with open(os.devnull, "w", encoding="utf-8") as null:
        with contextlib.redirect_stdout(null), contextlib.redirect_stderr(null):
            yield


def main() -> None:
    """Run fixed 30-channel sensitivities and compare with primary 32-channel calls."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--private-config", type=Path, default=DEFAULT_PRIVATE_CONFIG)
    args = parser.parse_args()
    contract = load_audit_contract(args.contract)
    config = load_private_config(args.private_config)
    output_dir = REPO_ROOT / config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    sensitivity = contract["reference_provenance_sensitivity"]
    excluded = list(sensitivity["excluded_channels"])
    if excluded != ["M1", "M2"] or sensitivity["full_prep_enabled"] is not False:
        raise ValueError("M1/M2 sensitivity contract changed unexpectedly.")

    primary_path = output_dir / "global_bad_method_runs.csv"
    primary = pd.read_csv(primary_path)
    primary = primary[(primary["repeat"] == 1) & (primary["success"] == True)]  # noqa: E712
    primary = primary[
        ((primary["method"] == "mne_lof") & (primary["variant"] == "primary"))
        | (primary["method"] == "pyprep_noisy_channels")
    ]
    baseline = {
        (row.source_filename, row.method): set(json.loads(row.bad_all_json)) - set(excluded)
        for row in primary.itertuples(index=False)
    }

    prepared = contract["prepared_continuous_input"]
    method_root = contract["global_bad_channel_methods"]
    lof_settings = method_root["lof"]["primary"]
    noisy_settings = method_root["noisy_channels"]
    seed = int(method_root["random_seed"])
    rows: list[dict[str, Any]] = []
    for index, recording in enumerate(config["recordings"], start=1):
        filename = recording["source_filename"]
        print(f"[{index}/{len(config['recordings'])}] {filename}", flush=True)
        raw = prepare_continuous_eeg(
            REPO_ROOT / config["raw_root"] / filename,
            channel_types=prepared["channel_types"],
            resample_hz=float(prepared["resample_hz"]),
        )
        raw.drop_channels(excluded)
        calls: dict[str, set[str]] = {}
        with quiet_third_party_output():
            lof = run_lof(
                raw,
                n_neighbors=int(lof_settings["n_neighbors"]),
                metric=str(lof_settings["metric"]),
                threshold=float(lof_settings["threshold"]),
            )
            noisy = run_noisy_channels(raw, noisy_settings, seed)
        calls["mne_lof"] = set(lof["bad_all"])
        calls["pyprep_noisy_channels"] = set(noisy["bad_all"])
        for method, current in calls.items():
            original = baseline[(filename, method)]
            rows.append(
                {
                    "source_filename": filename,
                    "method": method,
                    "primary_32ch_scalp_calls_json": channel_list_json(original),
                    "sensitivity_30ch_calls_json": channel_list_json(current),
                    "primary_32ch_scalp_count": len(original),
                    "sensitivity_30ch_count": len(current),
                    "retained_calls_json": channel_list_json(original & current),
                    "lost_calls_json": channel_list_json(original - current),
                    "gained_calls_json": channel_list_json(current - original),
                    "exact_scalp_call_agreement": original == current,
                }
            )

    table = pd.DataFrame(rows)
    paths = {
        "comparison": output_dir / "m1_m2_detector_sensitivity.csv",
        "manifest": output_dir / "m1_m2_detector_sensitivity_manifest.json",
        "summary": output_dir / "m1_m2_detector_sensitivity_summary.md",
    }
    suffixes = contract["safety"]["allowed_output_suffixes"]
    for path in paths.values():
        validate_audit_output_path(path, output_dir, suffixes)
    table.to_csv(paths["comparison"], index=False)
    aggregates = table.groupby("method").agg(
        files=("source_filename", "count"),
        exact_files=("exact_scalp_call_agreement", "sum"),
        primary_scalp_calls=("primary_32ch_scalp_count", "sum"),
        sensitivity_calls=("sensitivity_30ch_count", "sum"),
    )
    lines = ["| method | files | exact scalp sets | 32-ch scalp calls | 30-ch calls |", "|---|---:|---:|---:|---:|"]
    for method, row in aggregates.iterrows():
        lines.append(
            f"| {method} | {int(row.files)} | {int(row.exact_files)} | {int(row.primary_scalp_calls)} | {int(row.sensitivity_calls)} |"
        )
    paths["summary"].write_text(
        "# M1/M2 detector-pool sensitivity\n\n"
        + "\n".join(lines)
    + "\n\nThis sensitivity isolates detector-pool behavior only. The public contract separately records the accepted 30-channel source-pool decision; acquisition-reference roles remain unavailable.\n",
        encoding="utf-8",
    )
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": contract["contract_version"],
        "versions": {
            "python": platform.python_version(),
            "mne": mne.__version__,
            "pyprep": pyprep.__version__,
            "numpy": np.__version__,
        },
        "excluded_detector_channels": excluded,
        "full_prep_rerun": False,
        "signal_derivatives_written": 0,
        "epochs_constructed": 0,
        "outputs": [str(path.relative_to(REPO_ROOT)) for path in paths.values()],
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"M1/M2 detector sensitivity complete: {output_dir}")


if __name__ == "__main__":
    main()
