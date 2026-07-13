#!/usr/bin/env python3
"""Run the saved DEMI production continuous-preprocessing validation cohort.

Why this script exists
----------------------
Scripts 00--12 are evidence and audit programs. This script is the first
production continuous-preprocessing driver. It exercises saved derivatives,
atomic publication, objective stop routing, resumability, cache invalidation,
runtime/storage measurement, and derivative reopening on a deterministic
six-recording validation cohort before any full 95-recording run is authorized.

Inputs
------
The script reads raw EDFs from ``_Data/eeg/raw/``, the tracked
``continuous_preprocessing_config_v1.yaml``, and existing local QC/audit CSVs
used only for technical cohort selection. Raw EDFs are opened read-only and
re-hashed after processing.

Outputs
-------
New versioned per-recording FIF/JSON results and run-level JSON/Markdown
evidence are written below the ignored directory
``_Data/eeg/mne_preprocessing/continuous_validation_v1/``. Ordinary terminal
files retain one canonical post-ICA continuous FIF plus the ICA object; review
stops retain a pre-ICA FIF and ICA evidence without an automatic post-ICA file.

Explicit non-goals and safety boundaries
----------------------------------------
This script has no all-recording option. It does not modify or write EDF, repair
events, construct epochs, run AutoReject, compute CSD, change event-policy
eligibility, or make participant inclusion/exclusion decisions. Processing is
strictly sequential. A valid completed or stopped recording is skipped only
after source/config/code/environment and artifact hashes are verified.

Run from the repository root
----------------------------
List the deterministic cohort without writing derivatives::

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py --list-cohort

Run the complete validation cohort::

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py

Use ``--max-files 2`` for a bounded interruption/resume check. Use
``--recording 'demi_01 Data.edf' --force`` only for an explicit focused
recomputation; the prior terminal result is preserved in local history.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
EEG_ANALYSIS_DIR = Path(__file__).resolve().parent
if str(EEG_ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(EEG_ANALYSIS_DIR))

from continuous_preprocessing.contracts import load_config  # noqa: E402
from continuous_preprocessing.runner import (  # noqa: E402
    build_validation_cohort,
    run_validation_cohort,
)


def parse_args() -> argparse.Namespace:
    """Parse the bounded validation-only command-line interface."""

    parser = argparse.ArgumentParser(
        description="Run the saved sequential DEMI continuous-preprocessing validation cohort."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=EEG_ANALYSIS_DIR / "continuous_preprocessing_config_v1.yaml",
        help="Tracked production configuration (default: continuous_preprocessing_config_v1.yaml).",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=None,
        help="Process only the first N selected files for a bounded resume test.",
    )
    parser.add_argument(
        "--recording",
        default=None,
        help="Process one exact filename from the deterministic validation cohort.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Explicitly preserve and recompute existing selected terminal result(s).",
    )
    parser.add_argument(
        "--list-cohort",
        action="store_true",
        help="Print deterministic cohort evidence and exit without writing derivatives.",
    )
    return parser.parse_args()


def main() -> int:
    """Run cohort listing or sequential saved validation processing."""

    args = parse_args()
    config_path = args.config if args.config.is_absolute() else (REPO_ROOT / args.config)
    config = load_config(config_path.resolve())
    if args.list_cohort:
        cohort = build_validation_cohort(REPO_ROOT, config)
        print(json.dumps(cohort, indent=2, sort_keys=True))
        return 0
    aggregate = run_validation_cohort(
        repo_root=REPO_ROOT,
        config_path=config_path.resolve(),
        max_files=args.max_files,
        only_recording=args.recording,
        force=args.force,
    )
    print(f"Run directory: {aggregate['run_directory']}")
    print(f"Invocation counts: {aggregate['invocation_counts']}")
    print(f"Current cohort terminal counts: {aggregate['current_cohort_terminal_counts']}")
    return 0 if aggregate["invocation_counts"].get("failed", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
