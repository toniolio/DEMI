"""Compare explicit 60-Hz removal before a proposed 0.5--45 Hz EEG filter.

This compact audit evaluates two in-memory branches on a private deterministic
set of ordinary and high-60-Hz continuous recordings:

``A``
    No explicit line removal, followed by the same proposed 0.5--45 Hz MNE FIR.
``B``
    A 60-Hz ``spectrum_fit`` removal, followed by that identical FIR.

The full continuous EEG-only input uses the approved ``standard_1005`` montage
and a 250-Hz audit-only resampling rate. Channel-level outputs report residual
60-Hz PSD, theta/alpha/beta PSD, final-branch relative changes, waveform RMS
difference/correlation, and the exact candidate FIR response. A fixed edge
margin is excluded from measurements.

Outputs are CSV/JSON/Markdown evidence under the local audit directory. The
script does not choose a production default, save filtered data, alter EDFs,
fit ICA, construct epochs, or implement the production preprocessing pipeline.

Run from the repository root with::

    ./.venv/bin/python analysis/eeg_mne/11_compare_line_noise_filters.py
"""

from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import mne
import numpy as np
import pandas as pd
import scipy
import yaml

from preprocessing_parameter_audit import (
    compare_line_noise_branches,
    load_audit_contract,
    prepare_continuous_eeg,
    validate_audit_output_path,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = REPO_ROOT / "analysis/eeg_mne/preprocessing_parameter_audit_contract_v1.yaml"
DEFAULT_PRIVATE_CONFIG = REPO_ROOT / "_Private/config/preprocessing_parameter_evidence_pass_v1.yaml"


def load_private_config(path: Path) -> dict[str, Any]:
    """Load the local file selection without exposing it to public config."""

    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    selected = config.get("line_noise_recordings", [])
    if not selected or len(selected) != len(set(selected)):
        raise ValueError("Private config must list unique line-noise recordings.")
    return config


def summarize_file(filename: str, rows: pd.DataFrame) -> dict[str, Any]:
    """Aggregate channel-level branch metrics without hiding extrema."""

    summary: dict[str, Any] = {
        "source_filename": filename,
        "channel_count": len(rows),
        "median_analysis_b_minus_a_rms_proportion": rows[
            "analysis_b_minus_a_rms_proportion"
        ].median(),
        "maximum_analysis_b_minus_a_rms_proportion": rows[
            "analysis_b_minus_a_rms_proportion"
        ].max(),
        "median_analysis_a_b_correlation": rows["analysis_a_b_correlation"].median(),
        "minimum_analysis_a_b_correlation": rows["analysis_a_b_correlation"].min(),
    }
    for band in ("theta", "alpha", "beta", "line_60"):
        for branch in ("raw", "notch_only", "analysis_a", "analysis_b"):
            column = f"{branch}_{band}_mean_psd"
            summary[f"median_{column}"] = rows[column].median()
            summary[f"maximum_{column}"] = rows[column].max()
        change = rows[f"analysis_b_vs_a_{band}_relative_change"].replace(
            [np.inf, -np.inf], np.nan
        )
        summary[f"median_analysis_b_vs_a_{band}_relative_change"] = change.median()
        summary[f"maximum_absolute_analysis_b_vs_a_{band}_relative_change"] = (
            change.abs().max()
        )
    return summary


def main() -> None:
    """Run both line/filter branches and write local comparison evidence."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--private-config", type=Path, default=DEFAULT_PRIVATE_CONFIG)
    args = parser.parse_args()
    contract = load_audit_contract(args.contract)
    config = load_private_config(args.private_config)
    output_dir = REPO_ROOT / config["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    suffixes = contract["safety"]["allowed_output_suffixes"]
    paths = {
        "channels": output_dir / "line_noise_filter_channel_comparison.csv",
        "files": output_dir / "line_noise_filter_file_summary.csv",
        "response": output_dir / "line_noise_filter_response.csv",
        "manifest": output_dir / "line_noise_filter_manifest.json",
        "summary": output_dir / "line_noise_filter_summary.md",
    }
    for path in paths.values():
        validate_audit_output_path(path, output_dir, suffixes)

    prepared = contract["prepared_continuous_input"]
    settings = contract["line_noise_filter_comparison"]
    all_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    response_rows: list[dict[str, Any]] = []
    for index, filename in enumerate(config["line_noise_recordings"], start=1):
        print(f"[{index}/{len(config['line_noise_recordings'])}] {filename}", flush=True)
        raw = prepare_continuous_eeg(
            REPO_ROOT / config["raw_root"] / filename,
            channel_types=prepared["channel_types"],
            resample_hz=float(settings["resample_hz"]),
        )
        channel_rows, response = compare_line_noise_branches(raw, settings)
        for row in channel_rows:
            row["source_filename"] = filename
            row["prepared_sfreq_hz"] = float(raw.info["sfreq"])
            row["recording_duration_seconds"] = float(raw.times[-1])
        table = pd.DataFrame(channel_rows)
        all_rows.extend(channel_rows)
        file_rows.append(summarize_file(filename, table))
        response_rows.append({"source_filename": filename, **response})

    channel_table = pd.DataFrame(all_rows)
    file_table = pd.DataFrame(file_rows)
    response_table = pd.DataFrame(response_rows)
    channel_table.to_csv(paths["channels"], index=False)
    file_table.to_csv(paths["files"], index=False)
    response_table.to_csv(paths["response"], index=False)

    theta_change = file_table["median_analysis_b_vs_a_theta_relative_change"].abs().max()
    alpha_change = file_table["median_analysis_b_vs_a_alpha_relative_change"].abs().max()
    beta_change = file_table["median_analysis_b_vs_a_beta_relative_change"].abs().max()
    line_a = file_table["median_analysis_a_line_60_mean_psd"].median()
    line_b = file_table["median_analysis_b_line_60_mean_psd"].median()
    line_relative = (line_b - line_a) / line_a if line_a else np.nan
    summary = f"""# Line-noise/filter comparison

Generated: {datetime.now(timezone.utc).isoformat()}

Files: {len(file_table)}. Channels per file: {int(file_table['channel_count'].median())}. Both branches used the same proposed {settings['analysis_high_pass_hz']}--{settings['analysis_low_pass_hz']} Hz MNE FIR; branch B first used explicit {settings['line_frequency_hz']} Hz {settings['notch_method']} removal.

Largest absolute file-median B-vs-A relative changes: theta {theta_change:.6g}, alpha {alpha_change:.6g}, beta {beta_change:.6g}.

Across-file median final 60-Hz mean PSD: branch A {line_a:.6g}; branch B {line_b:.6g}; relative change {line_relative:.6g}.

Candidate FIR length: {response_table['filter_length_seconds'].iloc[0]:.3f} s. Response at 60 Hz: {response_table['response_db_60_hz'].iloc[0]:.3f} dB. Measurements exclude {settings['measurement_edge_exclusion_seconds']} s at each recording edge.

These are audit-only comparisons. No filtered signal was written, and the 0.5--45 Hz band remains a proposal rather than an implemented production parameter.
"""
    paths["summary"].write_text(summary, encoding="utf-8")
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": contract["contract_version"],
        "versions": {
            "python": platform.python_version(),
            "mne": mne.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "recordings": config["line_noise_recordings"],
        "parameters": settings,
        "signal_derivatives_written": 0,
        "epochs_constructed": 0,
        "authoritative_signal_changed": False,
        "outputs": [str(path.relative_to(REPO_ROOT)) for path in paths.values()],
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Line-noise/filter comparison complete: {output_dir}")


if __name__ == "__main__":
    main()
