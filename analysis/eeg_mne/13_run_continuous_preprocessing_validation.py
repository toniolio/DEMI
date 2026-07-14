#!/usr/bin/env python3
"""Run saved DEMI production continuous preprocessing.

Why this script exists
----------------------
Scripts 00--12 are evidence and audit programs. This script is the first
production continuous-preprocessing driver. Its default remains the separate
six-recording validation cohort. The explicit ``--all-recordings`` mode uses
the same validated package and configuration to process all 95 inventoried
readable EDF files under a separate production root.

Inputs
------
The script reads raw EDFs from ``_Data/eeg/raw/`` and the tracked
``continuous_preprocessing_config_v1.yaml``. Validation selection uses existing
technical audit evidence. Production selection uses the script-00 raw EDF
inventory only; event or analytic eligibility never gates continuous work.

Outputs
-------
New versioned per-recording FIF/JSON results and run-level JSON/Markdown
evidence are written below the ignored directory
``_Data/eeg/mne_preprocessing/continuous_validation_v1/`` for validation or
``_Data/eeg/mne_preprocessing/continuous_v1/`` for production. Ordinary
terminal files retain one canonical post-ICA continuous FIF plus the ICA
object; review stops retain a pre-ICA FIF and ICA evidence without an automatic
post-ICA file.

Explicit non-goals and safety boundaries
----------------------------------------
It does not modify or write EDF, repair events, construct epochs, run
AutoReject, compute CSD, change event-policy eligibility, or make participant
inclusion/exclusion decisions. Processing is strictly sequential. A valid
completed or stopped recording is skipped only after source/config/code/
environment/execution and artifact hashes are verified. Force recomputation
requires one exact recording name.

Run from the repository root
----------------------------
List the deterministic cohort without writing derivatives::

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py --list-cohort

Run the complete validation cohort::

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py

Run the explicitly authorized 95-file production surface::

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/eeg_mne/13_run_continuous_preprocessing_validation.py --all-recordings

Use ``--max-files 2`` for a bounded interruption/resume check. Use
``--all-recordings --recording 'demi_01 Data.edf'`` for a production smoke.
Add ``--force`` only for an explicit focused recomputation; the prior terminal
result is preserved in local history. Use ``--all-recordings --verify-current``
to reopen and rescan all current production FIF/ICA artifacts without
processing signals again. The one-time bounded
``--all-recordings --repair-historical-ica-routing`` mode reuses only the saved
pre-ICA and ICA artifacts from the 25 superseded over-two guardrail stops.
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
    build_production_surface,
    build_validation_cohort,
    run_historical_ica_routing_repair,
    run_continuous_surface,
    verify_current_surface,
)


def parse_args() -> argparse.Namespace:
    """Parse the explicit validation/production command-line interface."""

    parser = argparse.ArgumentParser(
        description="Run saved sequential DEMI continuous preprocessing."
    )
    parser.add_argument(
        "--all-recordings",
        action="store_true",
        help="Use the authorized 95-file production surface and continuous_v1 root.",
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
        help="Process one exact filename from the selected surface.",
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
    parser.add_argument(
        "--list-production-surface",
        action="store_true",
        help="Print the deterministic 95-file production surface and exit.",
    )
    parser.add_argument(
        "--verify-current",
        action="store_true",
        help="Reopen and rescan every current production FIF/ICA artifact.",
    )
    parser.add_argument(
        "--repair-historical-ica-routing",
        action="store_true",
        help="Repair only the 25 superseded over-two ICA stops from saved artifacts.",
    )
    return parser.parse_args()


def main() -> int:
    """List, verify, or sequentially process the selected saved surface."""

    args = parse_args()
    config_path = args.config if args.config.is_absolute() else (REPO_ROOT / args.config)
    config = load_config(config_path.resolve())
    if args.repair_historical_ica_routing:
        if not args.all_recordings:
            raise ValueError("--repair-historical-ica-routing requires --all-recordings.")
        if (
            args.recording
            or args.max_files
            or args.force
            or args.verify_current
            or args.list_cohort
            or args.list_production_surface
        ):
            raise ValueError("ICA routing repair cannot be combined with other processing selectors.")
        aggregate = run_historical_ica_routing_repair(
            repo_root=REPO_ROOT, config_path=config_path.resolve()
        )
        print(f"Repair run directory: {aggregate['run_directory']}")
        print(f"Invocation counts: {aggregate['invocation_counts']}")
        print(f"Current surface terminal counts: {aggregate['current_surface_terminal_counts']}")
        return 0 if aggregate["invocation_counts"].get("failed", 0) == 0 else 1
    if args.list_cohort:
        if args.all_recordings or args.list_production_surface or args.verify_current:
            raise ValueError("--list-cohort cannot be combined with another surface action.")
        cohort = build_validation_cohort(REPO_ROOT, config)
        print(json.dumps(cohort, indent=2, sort_keys=True))
        return 0
    if args.list_production_surface:
        if args.all_recordings or args.verify_current:
            raise ValueError("--list-production-surface is a standalone read-only action.")
        surface = build_production_surface(REPO_ROOT, config)
        print(json.dumps(surface, indent=2, sort_keys=True))
        return 0
    if args.verify_current:
        if not args.all_recordings:
            raise ValueError("--verify-current requires --all-recordings.")
        if args.recording or args.max_files or args.force:
            raise ValueError("--verify-current cannot be combined with processing selectors.")
        verification = verify_current_surface(
            repo_root=REPO_ROOT,
            config_path=config_path.resolve(),
            surface_mode="production",
        )
        print(f"Verification path: {verification['verification_path']}")
        print(
            "All current FIF/ICA artifacts valid: "
            f"{verification['all_current_fif_and_ica_artifacts_reopened_and_valid']}"
        )
        return (
            0
            if verification["all_current_fif_and_ica_artifacts_reopened_and_valid"]
            else 1
        )
    aggregate = run_continuous_surface(
        repo_root=REPO_ROOT,
        config_path=config_path.resolve(),
        surface_mode="production" if args.all_recordings else "validation",
        max_files=args.max_files,
        only_recording=args.recording,
        force=args.force,
    )
    print(f"Run directory: {aggregate['run_directory']}")
    print(f"Invocation counts: {aggregate['invocation_counts']}")
    print(f"Current surface terminal counts: {aggregate['current_surface_terminal_counts']}")
    return 0 if aggregate["invocation_counts"].get("failed", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
