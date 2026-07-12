"""Audit historical DEMI PREP, interpolation, and reference provenance.

This script exists to recover the exact historical preprocessing behavior from
tracked source and settings before the active MNE preprocessing policy is
finalized. It reads the pinned external PyPREP implementation, downstream R
settings/import code, the raw EDF manifest, and any locally recoverable old
``prep_info.csv`` files. It writes source-cited CSV/JSON/Markdown evidence under
the configured local ``_Data`` audit directory.

Inputs
------
* ``analysis/eeg_mne/preprocessing_parameter_audit_contract_v1.yaml``;
* ``external/DEMI_EEG_Pipeline/eeg_pipeline.py`` and ``edf2bids.py``;
* ``_Scripts/_settings.R``, ``04_import_eeg.R``, and ``07_gam_prep.R``;
* ``_Data/eeg/manifest/raw_eeg_file_manifest.csv``;
* optional historical ``prep_info.csv`` files if they remain local.

Outputs
-------
* ``historical_prep_code_evidence.csv``;
* ``historical_prep_output_recovery.csv``;
* ``acquisition_reference_evidence.csv``;
* ``historical_prep_audit_manifest.json``;
* ``historical_prep_audit_summary.md``.

This script does not read EEG samples, preprocess signals, write EDF/FIF files,
construct epochs, or treat historical output as ground truth. Direct code
evidence, bounded inference, and unavailable historical artifacts are labelled
separately.

Run from the repository root with::

    ./.venv/bin/python analysis/eeg_mne/09_audit_historical_prep_and_reference.py
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from preprocessing_parameter_audit import (
    load_audit_contract,
    parse_python_assignments,
    parse_r_scalar_assignments,
    validate_audit_output_path,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = REPO_ROOT / "analysis/eeg_mne/preprocessing_parameter_audit_contract_v1.yaml"
DEFAULT_OUTPUT = REPO_ROOT / "_Data/eeg/mne_preprocessing/preprocessing_parameter_audit_v1"


def source_line_contains(path: Path, start: int, end: int, text: str) -> bool:
    """Return whether a required source fragment occurs within an exact line span.

    The check prevents the audit narrative from drifting away from the pinned
    historical source. It reads text only and has no side effects.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    return text in "\n".join(lines[start - 1 : end])


def parse_metadata_dict(path: Path) -> dict[str, Any]:
    """Recover the literal ``metadata`` dictionary without executing the converter."""

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "metadata" for target in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError(f"No literal metadata dictionary found in {path}")


def git_head(path: Path) -> str:
    """Return the checked-out historical submodule commit without changing Git state."""

    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def build_code_evidence(root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build direct historical PREP/interpolation evidence from pinned code."""

    pipeline = root / "external/DEMI_EEG_Pipeline/eeg_pipeline.py"
    settings_path = root / "_Scripts/_settings.R"
    values = parse_python_assignments(
        pipeline, ("perform_csd", "interpolate_bads", "max_interpolated", "seed")
    )
    r_values = parse_r_scalar_assignments(
        settings_path, ("drop_interpolated", "using_csd", "exclude_bad_by_emg")
    )
    required = {
        "perform_csd": True,
        "interpolate_bads": True,
        "max_interpolated": 0.25,
        "seed": 530453080,
    }
    if values != required:
        raise ValueError(f"Historical Python settings changed unexpectedly: {values}")
    if r_values.get("drop_interpolated") is not False or r_values.get("using_csd") is not True:
        raise ValueError(f"Historical R settings changed unexpectedly: {r_values}")

    checks = (
        (pipeline, 26, 29, "from pyprep.reference import Reference"),
        (pipeline, 257, 271, 'method="spectrum_fit"'),
        (pipeline, 278, 290, "reference.perform_reference()"),
        (pipeline, 292, 297, "if not interpolate_bads:"),
        (pipeline, 307, 317, "reference.interpolated_channels"),
        (pipeline, 328, 334, "prop_interpolated"),
        (pipeline, 375, 380, "raw_prepped.drop_channels(remaining_bad)"),
        (root / "_Scripts/04_import_eeg.R", 65, 88, "if (drop_interpolated)"),
    )
    for path, start, end, fragment in checks:
        if not source_line_contains(path, start, end, fragment):
            raise ValueError(f"Expected historical evidence missing: {path}:{start} {fragment}")

    rows = [
        {
            "finding_id": "historical_pyprep_api",
            "evidence_class": "direct_code_evidence",
            "finding": "The script imported PrepPipeline, NoisyChannels, removeTrend, and Reference; its active PREP path manually called Reference.perform_reference rather than PrepPipeline.fit.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "26-29;253-290",
        },
        {
            "finding_id": "historical_line_noise",
            "evidence_class": "direct_code_evidence",
            "finding": "Line removal detrended EEG in microvolts, used spectrum_fit at 60-Hz harmonics with a 10-s filter length, 2-Hz multitaper bandwidth, and p=0.01, then restored the residual trend.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "253-272",
        },
        {
            "finding_id": "historical_reference",
            "evidence_class": "direct_code_evidence",
            "finding": "Reference received all EEG channels as both ref_chs and reref_chs, enabled RANSAC, used fixed seed 530453080, and performed iterative robust average re-referencing.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "38;278-290",
        },
        {
            "finding_id": "historical_interpolation_enabled",
            "evidence_class": "direct_code_evidence",
            "finding": "interpolate_bads was True. The fallback restoring EEG_before_interpolation ran only when that setting was false, so the pinned active path retained Reference interpolation.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "52-54;292-297",
        },
        {
            "finding_id": "historical_channel_logs",
            "evidence_class": "direct_code_evidence",
            "finding": "prep_info.csv recorded noisy_channels_original['bad_all'] as initial_bad, interpolated_channels as interpolated, and still_noisy_channels as remaining_bad, with counts for each.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "307-317;420-439",
        },
        {
            "finding_id": "historical_25pct",
            "evidence_class": "direct_code_evidence",
            "finding": "More than 25% interpolated produced a warning and routed the final output to bad/too_noisy; it did not stop PREP, filtering, ICA, CSD, or export.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "328-334;337-395",
        },
        {
            "finding_id": "historical_drop_interpolated_false",
            "evidence_class": "direct_code_evidence",
            "finding": "drop_interpolated was False downstream. The R import therefore did not add already interpolated channel names to the drop list; it did not disable upstream interpolation.",
            "source": "_Scripts/_settings.R;_Scripts/04_import_eeg.R",
            "line_span": "8;65-88",
        },
        {
            "finding_id": "historical_csd_interaction",
            "evidence_class": "direct_code_evidence",
            "finding": "using_csd was True. The Python CSD branch dropped still-noisy channels, and R deliberately emptied its remaining-bad drop table because those channels were assumed absent already.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py;_Scripts/_settings.R;_Scripts/04_import_eeg.R",
            "line_span": "375-380;11;65-77",
        },
        {
            "finding_id": "historical_coupling",
            "evidence_class": "bounded_inference_from_code_order",
            "finding": "Line removal, repeated detection, robust reference, interpolation, post-interpolation detection, and downstream CSD channel dropping were operationally coupled in one subject function and one shared Reference object.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py",
            "line_span": "253-395",
        },
    ]
    return pd.DataFrame(rows), {**values, **r_values}


def build_output_recovery(root: Path) -> pd.DataFrame:
    """Inventory expected historical per-recording logs and signal outputs locally."""

    candidates = [
        root / "_Data/eeg/prep_info.csv",
        root / "external/DEMI_EEG_Pipeline/output/prep_info.csv",
        root / "legacy/dissertation/_Preprocessing/output/prep_info.csv",
    ]
    rows: list[dict[str, Any]] = []
    for path in candidates:
        exists = path.is_file()
        row_count = None
        columns = None
        if exists:
            table = pd.read_csv(path)
            row_count = len(table)
            columns = json.dumps(list(table.columns))
        rows.append(
            {
                "artifact": str(path.relative_to(root)),
                "exists": exists,
                "row_count": row_count,
                "columns_json": columns,
                "evidence_class": "recoverable_historical_output" if exists else "unavailable_historical_output",
                "interpretation": (
                    "Per-recording initial/interpolated/remaining bad lists are locally recoverable."
                    if exists
                    else "Code proves the expected schema, but participant-level channel lists are not recoverable from the current checkout."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_acquisition_reference_evidence(root: Path) -> pd.DataFrame:
    """Separate direct EDF/code facts from unsupported acquisition-reference claims."""

    manifest_path = root / "_Data/eeg/manifest/raw_eeg_file_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    channel_sets = manifest["channel_names_json"].map(json.loads)
    all_same = channel_sets.map(tuple).nunique() == 1
    channels = channel_sets.iloc[0]
    explicit_reference_labels = [
        channel for channel in channels if channel.upper() in {"REF", "REFERENCE", "A1", "A2"}
    ]
    metadata = parse_metadata_dict(root / "external/DEMI_EEG_Pipeline/edf2bids.py")
    rows = [
        {
            "finding_id": "edf_channel_surface",
            "evidence_class": "direct_edf_manifest_evidence",
            "finding": f"All {len(manifest)} EDF rows have the same {len(channels)} labels: {all_same}.",
            "source": "_Data/eeg/manifest/raw_eeg_file_manifest.csv",
        },
        {
            "finding_id": "explicit_reference_channel",
            "evidence_class": "direct_edf_manifest_evidence",
            "finding": f"No separately labelled REF/REFERENCE/A1/A2 channel is present: {explicit_reference_labels}.",
            "source": "_Data/eeg/manifest/raw_eeg_file_manifest.csv",
        },
        {
            "finding_id": "m1_m2_recorded",
            "evidence_class": "direct_edf_manifest_evidence",
            "finding": f"M1 and M2 are both recorded signal labels in every EDF: {all(name in channels for name in ('M1', 'M2'))}.",
            "source": "_Data/eeg/manifest/raw_eeg_file_manifest.csv",
        },
        {
            "finding_id": "converter_hardware_metadata",
            "evidence_class": "direct_code_evidence_with_placeholders",
            "finding": f"The converter names {metadata['Manufacturer']} {metadata['ManufacturersModelName']} and {metadata['SoftwareVersions']}; cap and placement strings explicitly require confirmation.",
            "source": "external/DEMI_EEG_Pipeline/edf2bids.py:68-93",
        },
        {
            "finding_id": "converter_reference_claim",
            "evidence_class": "unconfirmed_code_metadata_not_acquisition_proof",
            "finding": f"The converter says EEGReference={metadata['EEGReference']!r} and EEGGround={metadata['EEGGround']!r}, but neighboring fields contain placeholders and no acquisition record corroborates these values.",
            "source": "external/DEMI_EEG_Pipeline/edf2bids.py:68-93",
        },
        {
            "finding_id": "historical_rereference_behavior",
            "evidence_class": "direct_code_evidence",
            "finding": "Historical PyPREP included M1/M2 among all EEG reference and rereference channels; later GAM preparation dropped M1/M2 only from the compact feature table.",
            "source": "external/DEMI_EEG_Pipeline/eeg_pipeline.py:143-145,278-290;_Scripts/07_gam_prep.R:22-23,80-83",
        },
        {
            "finding_id": "acquisition_reference_unresolved",
            "evidence_class": "unavailable_primary_source",
            "finding": "The EDF/channel/code surface does not identify the amplifier's online acquisition reference or establish whether M1/M2 served as online, linked, ordinary mastoid, or derived channels.",
            "source": "repository-wide acquisition/reference search",
        },
    ]
    return pd.DataFrame(rows)


def write_summary(
    path: Path,
    code: pd.DataFrame,
    recovery: pd.DataFrame,
    acquisition: pd.DataFrame,
) -> None:
    """Write a concise generated summary of direct, inferred, and unavailable evidence."""

    recoverable = int(recovery["exists"].sum())
    text = f"""# Historical PREP and reference audit

Generated: {datetime.now(timezone.utc).isoformat()}

## Historical interpolation conclusion

The pinned external pipeline enabled interpolation. PyPREP ``Reference.perform_reference()`` ran with ``interpolate_bads = True``; the downstream R setting ``drop_interpolated <- FALSE`` retained those interpolated channels rather than disabling interpolation.

## Coupling and 25% behavior

Line removal, detector iterations, robust average reference, interpolation, post-interpolation detection, and CSD channel dropping occurred inside one subject workflow. More than 25% interpolated was a warning/output-routing rule, not an early processing stop.

## Recoverable participant logs

Recoverable ``prep_info.csv`` artifacts: {recoverable} of {len(recovery)} searched locations. When absent, the schema is direct code evidence but old participant/channel lists are unavailable.

## Acquisition/reference boundary

M1 and M2 are recorded signals and no separately labelled reference channel appears in the EDF channel surface. The converter's common-average/AFz strings are unconfirmed metadata, not acquisition proof. The online acquisition reference and M1/M2 acquisition roles therefore remain unresolved from available primary sources.

## Evidence tables

* Historical findings: {len(code)}
* Output-recovery locations: {len(recovery)}
* Acquisition/reference findings: {len(acquisition)}

No EEG signal was read or written, and no historical output was treated as ground truth.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    """Run the source/manifest audit and write local tabular evidence."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    contract = load_audit_contract(args.contract)
    suffixes = contract["safety"]["allowed_output_suffixes"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "code": args.output_dir / "historical_prep_code_evidence.csv",
        "recovery": args.output_dir / "historical_prep_output_recovery.csv",
        "acquisition": args.output_dir / "acquisition_reference_evidence.csv",
        "manifest": args.output_dir / "historical_prep_audit_manifest.json",
        "summary": args.output_dir / "historical_prep_audit_summary.md",
    }
    for path in paths.values():
        validate_audit_output_path(path, args.output_dir, suffixes)

    code, settings = build_code_evidence(REPO_ROOT)
    recovery = build_output_recovery(REPO_ROOT)
    acquisition = build_acquisition_reference_evidence(REPO_ROOT)
    code.to_csv(paths["code"], index=False)
    recovery.to_csv(paths["recovery"], index=False)
    acquisition.to_csv(paths["acquisition"], index=False)
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "contract_version": contract["contract_version"],
        "historical_submodule_commit": git_head(REPO_ROOT / "external/DEMI_EEG_Pipeline"),
        "settings": settings,
        "signal_samples_read": 0,
        "signal_derivatives_written": 0,
        "epochs_constructed": 0,
        "outputs": [str(path.relative_to(REPO_ROOT)) for path in paths.values()],
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_summary(paths["summary"], code, recovery, acquisition)
    print(f"Historical PREP/reference audit complete: {args.output_dir}")


if __name__ == "__main__":
    main()
