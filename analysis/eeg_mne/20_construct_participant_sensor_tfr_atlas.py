#!/usr/bin/env python3
"""Construct the descriptive DEMI participant-level sensor TFR atlas.

This is a read-only visual-QC stage.  It opens only the accepted TFR-v1
response-onset and response-end dB arrays with NumPy memory mapping, then
averages existing trial-level dB cells within each participant's assigned
condition.  Physical-group participants use accepted overt trials; imagery
participants use accepted imagery trials only.  Their 738 final-overt bridge
trials remain in the accepted source product but are deliberately excluded
from this display mean.

Inputs:
    Accepted TFR-v1 manifests, axes, row metadata, and response-onset/end dB
    arrays; accepted behavioural-lineage metadata for condition/hand/QC page
    labels; the recovered historical BESA presentation-coordinate CSV; and the
    tracked atlas configuration.

Outputs:
    Two ignored multipage PDFs (overt and imagery), a participant/page
    manifest, trial/saturation summaries, plotting/configuration and source
    provenance records, source immutability evidence, validation JSON, state,
    log, and concise Markdown under
    ``_Data/eeg/participant_sensor_tfr_atlas_v1/``.

This script explicitly does not compute TFRs, make a baseline, apply CSD,
mirror sensors, modify montage/reference/interpolation/eligibility, select a
representative participant, generate group/ROI plots, fit a model, calculate a
contrast/association/p-value, or interpret physiological effects.

Run from the repository root:

    tools/run_participant_sensor_tfr_atlas_v1.sh

Use ``--verify-current`` to reopen and validate the finished PDFs and
manifests without rewriting any output.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from typing import Any, Mapping

import matplotlib

matplotlib.use("Agg")
from matplotlib import cm, colors, pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np
import pandas as pd

from participant_tfr_atlas import (
    ATLAS_NAMESPACE,
    ATLAS_SCHEMA_VERSION,
    EXPECTED_ASSIGNED_TRIALS,
    EXPECTED_BRIDGE_EXCLUDED,
    EXPECTED_IMAGERY_PARTICIPANTS,
    EXPECTED_OVERT_PARTICIPANTS,
    EXPECTED_PARTICIPANTS,
    PRIMARY_CHANNELS,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    compare_source_snapshots,
    load_and_validate_config,
    load_layout_positions,
    mean_memory_mapped_trials,
    output_descriptor,
    participant_page_table,
    pdf_page_count,
    saturation_accounting,
    select_assigned_condition_trials,
    sha256_file,
    source_snapshot,
    stable_json_hash,
)


CONFIG_PATH = Path("analysis/eeg_mne/participant_tfr_atlas_config_v1.yaml")
OUTPUT_ROOT = Path("_Data/eeg/participant_sensor_tfr_atlas_v1")
TFR_ROOT = Path("_Data/eeg/tfr_v1")
LINEAGE_PATH = Path("_Data/eeg/features_v1/metadata/behavioural_lineage.parquet")
RUN_MANIFEST = "atlas_run_manifest.json"
RUN_STATE = "atlas_run_state.json"


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp suitable for provenance files."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def print_progress(message: str) -> None:
    """Print one timestamped progress line immediately."""

    print(f"[{utc_now()}] {message}", flush=True)


def repo_root_from_script() -> Path:
    """Return the repository root and fail if the tracked atlas config is absent."""

    root = Path(__file__).resolve().parents[2]
    if not (root / ".git").exists() or not (root / CONFIG_PATH).is_file():
        raise RuntimeError("could_not_locate_demi_repository_root")
    return root


def parse_args() -> argparse.Namespace:
    """Parse the non-destructive validation-only command-line option."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify-current",
        action="store_true",
        help="Require a complete current atlas and validate/reuse it without rewriting.",
    )
    parser.add_argument(
        "--replace-current",
        action="store_true",
        help="Archive the current atlas before an explicitly authorized rendering-only republication.",
    )
    args = parser.parse_args()
    if args.verify_current and args.replace_current:
        parser.error("--verify-current cannot be combined with --replace-current")
    return args


def read_json(path: Path) -> dict[str, Any]:
    """Read one required JSON object."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"atlas_json_not_object:{path}")
    return value


def require_ignored_output_root(repo_root: Path, output_root: Path) -> None:
    """Require the target to be an ignored descendant of this repository."""

    resolved_repo = repo_root.resolve()
    resolved_output = output_root.resolve()
    if resolved_output == resolved_repo or resolved_repo not in resolved_output.parents:
        raise RuntimeError("atlas_unsafe_output_root")
    relative = resolved_output.relative_to(resolved_repo).as_posix()
    result = subprocess.run(
        ["git", "check-ignore", "--quiet", relative], cwd=repo_root, check=False
    )
    if result.returncode != 0:
        raise RuntimeError("atlas_output_root_not_ignored")


def validate_descriptor(path: Path, descriptor: Mapping[str, Any]) -> None:
    """Require a file to match its recorded size and content digest."""

    if not path.is_file() or path.stat().st_size != int(descriptor["size_bytes"]):
        raise RuntimeError(f"atlas_output_size_mismatch:{path}")
    if sha256_file(path) != descriptor["sha256"]:
        raise RuntimeError(f"atlas_output_hash_mismatch:{path}")


def load_tfr_authority(repo_root: Path) -> dict[str, Any]:
    """Validate TFR-v1 axes/manifests and expose only current dB-array paths."""

    root = repo_root / TFR_ROOT
    run_path = root / "tfr_run_manifest.json"
    run = read_json(run_path)
    if run.get("status") != "complete" or run.get("tfr_namespace") != "tfr_v1":
        raise RuntimeError("atlas_tfr_run_manifest_mismatch")
    axes: dict[str, pd.DataFrame] = {}
    for label in ("channels", "frequencies_cycles", "task_times"):
        descriptor = run["axes"][label]
        path = root / descriptor["relative_path"]
        validate_descriptor(path, descriptor)
        axes[label] = pd.read_csv(path)
    channels = axes["channels"]["channel"].astype(str).tolist()
    frequencies = axes["frequencies_cycles"]["frequency_hz"].to_numpy(dtype=float)
    times = axes["task_times"]["time_seconds"].to_numpy(dtype=float)
    if tuple(channels) != PRIMARY_CHANNELS:
        raise RuntimeError("atlas_tfr_channel_axis_mismatch")
    if not np.array_equal(frequencies, np.arange(4.0, 41.0)):
        raise RuntimeError("atlas_tfr_frequency_axis_mismatch")
    if not np.array_equal(times, np.arange(-50, 151, dtype=float) / 100.0):
        raise RuntimeError("atlas_tfr_time_axis_mismatch")
    recordings: dict[str, dict[str, Any]] = {}
    source_paths: list[Path] = [run_path]
    source_paths.extend(root / run["axes"][key]["relative_path"] for key in run["axes"])
    for item in run["recordings"]:
        manifest_path = root / item["manifest_path"]
        manifest = read_json(manifest_path)
        if manifest.get("status") != "complete" or manifest.get("fingerprint") != item["fingerprint"]:
            raise RuntimeError(f"atlas_tfr_recording_manifest_mismatch:{manifest_path}")
        recording_root = manifest_path.parent
        metadata_path = recording_root / "row_metadata.parquet"
        if not metadata_path.is_file():
            raise RuntimeError(f"atlas_tfr_row_metadata_missing:{metadata_path}")
        arrays = {}
        # Deliberately do not open raw power or red_on baseline arrays.
        for label in ("response_onset_db", "response_end_db"):
            descriptor = manifest["arrays"][label]
            path = recording_root / descriptor["relative_path"]
            if not path.is_file() or path.stat().st_size != int(descriptor["size_bytes"]):
                raise RuntimeError(f"atlas_tfr_db_array_size_mismatch:{path}")
            arrays[label] = path
            source_paths.append(path)
        source_paths.extend((manifest_path, metadata_path))
        stem = manifest["recording"]["recording_stem"]
        recordings[stem] = {
            "manifest": manifest,
            "manifest_path": manifest_path,
            "row_metadata_path": metadata_path,
            "arrays": arrays,
        }
    if len(recordings) != EXPECTED_PARTICIPANTS:
        raise RuntimeError("atlas_tfr_recording_count_mismatch")
    return {
        "root": root,
        "run": run,
        "run_path": run_path,
        "run_sha256": sha256_file(run_path),
        "channels": channels,
        "frequencies": frequencies,
        "times": times,
        "recordings": recordings,
        "source_paths": source_paths,
    }


def load_lineage(repo_root: Path, tfr: Mapping[str, Any]) -> pd.DataFrame:
    """Load accepted participant/QC metadata and verify its TFR row alignment."""

    path = repo_root / LINEAGE_PATH
    lineage = pd.read_parquet(path)
    if len(lineage) != 8_798 or lineage["canonical_event_key"].duplicated().any():
        raise RuntimeError("atlas_lineage_identity_mismatch")
    if lineage["recording_stem"].nunique() != EXPECTED_PARTICIPANTS:
        raise RuntimeError("atlas_lineage_recording_count_mismatch")
    for stem, frame in lineage.groupby("recording_stem", sort=False, observed=True):
        if stem not in tfr["recordings"]:
            raise RuntimeError(f"atlas_lineage_unknown_recording:{stem}")
        source = pd.read_parquet(tfr["recordings"][stem]["row_metadata_path"])
        if frame["canonical_event_key"].tolist() != source["canonical_event_key"].tolist():
            raise RuntimeError(f"atlas_lineage_row_order_mismatch:{stem}")
        if not np.array_equal(frame["tfr_row_index"].to_numpy(dtype=int), source["tfr_row_index"].to_numpy(dtype=int)):
            raise RuntimeError(f"atlas_lineage_tfr_index_mismatch:{stem}")
    return lineage


def participant_metadata(frame: pd.DataFrame, saturation: Mapping[str, Any]) -> dict[str, Any]:
    """Summarize one participant's selected rows for a concise page caption."""

    def one(column: str) -> Any:
        values = frame[column].drop_duplicates().tolist()
        if len(values) != 1:
            raise RuntimeError(f"atlas_participant_metadata_not_constant:{column}:{values}")
        return values[0]

    interpolation = int(one("interpolation_count")) if "interpolation_count" in frame else 0
    continuous_warning = bool(frame["continuous_qc_warning"].astype(bool).any())
    return {
        "behavioural_id": int(one("behavioural_id")),
        "recording_stem": str(one("recording_stem")),
        "assigned_condition": str(one("atlas_assigned_condition")),
        "performed_conditions": ", ".join(sorted(frame["performed_condition"].drop_duplicates().astype(str))),
        "contributing_trial_count": int(len(frame)),
        "strict_clean_trial_count": int(frame["strict_clean_eligibility"].astype(bool).sum()),
        "duration_warning_trial_count": int(frame["duration_warning_flag"].astype(bool).sum()),
        "original_handedness": str(one("original_handedness")),
        "analysis_hand": str(one("analysis_hand")),
        "motor_laterality_mapping": str(one("motor_laterality_mapping")),
        "mapping_source": str(one("mapping_source")),
        "mapping_override": bool(one("mapping_override")),
        "interpolation_count": interpolation,
        "continuous_qc_warning": continuous_warning,
        **dict(saturation),
    }


def archive_current_rendering(output_root: Path) -> Path:
    """Recoverably move a stale rendering aside before authorized republication.

    This never changes a TFR or participant summary. The old ignored atlas is
    preserved beside the atomic destination so a failed new rendering cannot
    overwrite it in place.
    """

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    archived = output_root.parent / f".{output_root.name}.previous-rendering-{stamp}"
    os.replace(output_root, archived)
    return archived


def render_participant_page(
    *,
    onset: np.ndarray,
    end: np.ndarray,
    metadata: Mapping[str, Any],
    channels: list[str],
    frequencies: np.ndarray,
    times: np.ndarray,
    positions: Mapping[str, Any],
    config: Mapping[str, Any],
) -> plt.Figure:
    """Render one paired-onset/end all-sensor page in the recovered visual style."""

    display = config["display"]
    minimum = float(display["db_minimum"])
    maximum = float(display["db_maximum"])
    width = float(display["panel_width_fraction"])
    height = float(display["panel_height_fraction"])
    gap = float(display["onset_end_gap_seconds"])
    end_start = float(times[-1]) + gap
    end_zero = end_start + 0.5
    end_stop = end_start + float(times[-1] - times[0])
    cmap = plt.get_cmap(display["palette"])
    norm = colors.Normalize(vmin=minimum, vmax=maximum, clip=True)
    figure = plt.figure(
        figsize=(float(display["figure_width_inches"]), float(display["figure_height_inches"])),
        facecolor="white",
    )
    figure.text(
        0.5, 0.985,
        f"Participant {metadata['behavioural_id']} — assigned {metadata['assigned_condition']} trials",
        ha="center", va="top", fontsize=15, fontweight="bold",
    )
    page_line = (
        f"n={metadata['contributing_trial_count']} (strict-clean {metadata['strict_clean_trial_count']}; "
        f"warnings {metadata['duration_warning_trial_count']})  |  performed: {metadata['performed_conditions']}  |  "
        f"original handedness: {metadata['original_handedness']}; analysis hand: {metadata['analysis_hand']}"
    )
    qc_line = (
        f"physical, unmirrored channels; interpolation={metadata['interpolation_count']}/30; "
        f"recording QC warning={'yes' if metadata['continuous_qc_warning'] else 'no'}  |  "
        f"current average-referenced TFR-v1 dB (not historical CSD)  |  "
        f"scale {minimum:g} to +{maximum:g} dB; saturation {metadata['saturated_percent']:.3f}%"
    )
    figure.text(0.5, 0.963, page_line, ha="center", va="top", fontsize=8.2)
    figure.text(0.5, 0.948, qc_line, ha="center", va="top", fontsize=7.5)

    guide_y = [float(item) for item in display["frequency_guides_hz"]]
    for index, channel in enumerate(channels):
        position = positions[channel]
        axis = figure.add_axes([
            position.x_fraction - width / 2.0,
            position.y_fraction - height / 2.0,
            width,
            height,
        ])
        axis.imshow(
            onset[index], origin="lower", aspect="auto", interpolation="nearest", cmap=cmap,
            norm=norm, extent=[float(times[0]), float(times[-1]), float(frequencies[0]), float(frequencies[-1])],
        )
        axis.imshow(
            end[index], origin="lower", aspect="auto", interpolation="nearest", cmap=cmap,
            norm=norm, extent=[end_start, end_stop, float(frequencies[0]), float(frequencies[-1])],
        )
        axis.set_xlim(float(times[0]), end_stop)
        axis.set_ylim(float(frequencies[0]), float(frequencies[-1]))
        axis.axvline(0.0, color="white", linewidth=0.75)
        axis.axvline(end_zero, color="white", linewidth=0.75)
        axis.axvspan(float(times[-1]), end_start, color="white", alpha=1.0, linewidth=0.0)
        for frequency in guide_y:
            axis.axhline(frequency, color="white", linewidth=0.65)
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.text(-0.08, 0.98, channel, transform=axis.transAxes, ha="right", va="top", fontsize=8)

    # A compact historical-style key mirrors the paired sensor-panel geometry:
    # two same-width epoch boxes with independent within-epoch t=0 markers.
    key = figure.add_axes([0.075, 0.064, 0.235, 0.105], label="atlas_explanatory_key")
    key.set_xlim(float(times[0]) - 0.24, end_stop + 0.04)
    key.set_ylim(float(frequencies[0]), float(frequencies[-1]))
    epoch_width = float(times[-1] - times[0])
    for epoch_start in (float(times[0]), end_start):
        key.add_patch(plt.Rectangle(
            (epoch_start, float(frequencies[0])), epoch_width,
            float(frequencies[-1] - frequencies[0]), fill=False,
            edgecolor="black", linewidth=0.8,
        ))
        for frequency in guide_y:
            key.plot([epoch_start, epoch_start + epoch_width], [frequency, frequency], color="black", linewidth=0.6)
    key.axvline(0.0, color="black", linewidth=0.8)
    key.axvline(end_zero, color="black", linewidth=0.8)
    label_centres = display["frequency_label_centres_hz"]
    key.set_yticks(
        [
            float(label_centres["theta"]), float(label_centres["alpha"]),
            float(label_centres["beta"]), float(label_centres["diagnostic_exploratory"]),
        ],
        ["theta\n4–8 Hz", "alpha\n9–12 Hz", "beta\n13–30 Hz", "diagnostic/\nexploratory\n31–40 Hz"],
        fontsize=5.8,
    )
    key.set_xticks([])
    key.tick_params(axis="y", length=0, pad=2)
    for spine in key.spines.values():
        spine.set_visible(False)
    key.text(0.0, -0.08, "t = 0", transform=key.get_xaxis_transform(), ha="center", va="top", fontsize=6.5, clip_on=False)
    key.text(end_zero, -0.08, "t = 0", transform=key.get_xaxis_transform(), ha="center", va="top", fontsize=6.5, clip_on=False)
    key.text(0.5, -0.27, "During", transform=key.get_xaxis_transform(), ha="center", va="top", fontsize=7.2, clip_on=False)
    key.text(end_start + 1.0, -0.27, "After", transform=key.get_xaxis_transform(), ha="center", va="top", fontsize=7.2, clip_on=False)

    colour_axis = figure.add_axes([0.76, 0.028, 0.18, 0.016])
    colourbar = figure.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=colour_axis, orientation="horizontal")
    colourbar.set_label("dB relative power", fontsize=8, labelpad=3)
    colourbar.ax.tick_params(labelsize=7, length=2)
    figure.text(0.76, 0.074, "white separators: 8.5, 12.5, 30.5 Hz", ha="left", va="bottom", fontsize=6.7)
    return figure


def source_paths(repo_root: Path, tfr: Mapping[str, Any], config: Mapping[str, Any]) -> list[Path]:
    """Return exact read-only sources whose before/after state is recorded."""

    paths = list(tfr["source_paths"])
    paths.extend(
        [
            repo_root / LINEAGE_PATH,
            repo_root / CONFIG_PATH,
            repo_root / config["paths"]["historical_layout"],
            repo_root / config["paths"]["historical_audit"],
            repo_root / "_Data/eeg/features_v1/feature_run_manifest.json",
            repo_root / "_Data/eeg/model_tables_v1/model_table_run_manifest.json",
        ]
    )
    return [path for path in paths if path.is_file()]


def authority_seed(
    *, repo_root: Path, tfr: Mapping[str, Any], config_sha256: str, source_before: list[dict[str, Any]]
) -> dict[str, Any]:
    """Create a deterministic source/config/code identity for cache-safe reuse."""

    code_paths = [
        repo_root / CONFIG_PATH,
        repo_root / "analysis/eeg_mne/participant_tfr_atlas.py",
        repo_root / "analysis/eeg_mne/20_construct_participant_sensor_tfr_atlas.py",
    ]
    return {
        "atlas_namespace": ATLAS_NAMESPACE,
        "atlas_config_sha256": config_sha256,
        "tfr_run_manifest_sha256": tfr["run_sha256"],
        "source_size_mtime_fingerprint": stable_json_hash(source_before),
        "code_sha256": {path.relative_to(repo_root).as_posix(): sha256_file(path) for path in code_paths},
    }


def validate_existing(output_root: Path, expected_authority: Mapping[str, Any]) -> dict[str, Any] | None:
    """Hash/reopen a matching finished atlas and return it without writing anything."""

    manifest_path = output_root / RUN_MANIFEST
    if not manifest_path.is_file():
        return None
    manifest = read_json(manifest_path)
    if manifest.get("status") != "complete" or manifest.get("authority") != dict(expected_authority):
        return None
    for descriptor in manifest["outputs"].values():
        validate_descriptor(output_root / descriptor["relative_path"], descriptor)
    pages = pd.read_csv(output_root / "participant_page_manifest.csv")
    if len(pages) != EXPECTED_PARTICIPANTS or pages["behavioural_id"].duplicated().any():
        raise RuntimeError("atlas_existing_page_manifest_mismatch")
    for group, expected in (("overt", EXPECTED_OVERT_PARTICIPANTS), ("imagery", EXPECTED_IMAGERY_PARTICIPANTS)):
        pdf = output_root / manifest["pdfs"][group]["relative_path"]
        if pdf_page_count(pdf) != expected:
            raise RuntimeError(f"atlas_existing_pdf_page_count_mismatch:{group}")
    validation = read_json(output_root / "validation/atlas_validation.json")
    if validation.get("status") != "pass":
        raise RuntimeError("atlas_existing_validation_not_pass")
    return manifest


def write_summary_markdown(summary: Mapping[str, Any]) -> str:
    """Return a concise local summary without scientific interpretation."""

    return (
        "# Participant sensor TFR atlas v1\n\n"
        f"Completed {summary['completed_at']} with {summary['participant_pages']} participant pages "
        f"and {summary['assigned_condition_trial_count']} contributing assigned-condition trials. "
        "This is a descriptive visual-QC atlas, not an inferential result.\n\n"
        "The pages use physical, unmirrored 30-channel labels and historical BESA-derived "
        "presentation coordinates only. They show current average-referenced TFR-v1 dB power, "
        "not the historical CSD surface. Every panel retains the stored 4–40 Hz and −0.5 to +1.5 s "
        "axes, with paired response-onset and response-end anchors.\n"
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Build or validation-reuse the atomically published atlas namespace."""

    started = time.perf_counter()
    repo_root = repo_root_from_script()
    config, config_sha256 = load_and_validate_config(repo_root / CONFIG_PATH)
    final_root = repo_root / OUTPUT_ROOT
    require_ignored_output_root(repo_root, final_root)
    tfr = load_tfr_authority(repo_root)
    sources = source_paths(repo_root, tfr, config)
    before = source_snapshot(sources, root=repo_root)
    authority = authority_seed(repo_root=repo_root, tfr=tfr, config_sha256=config_sha256, source_before=before)
    existing = validate_existing(final_root, authority)
    if existing is not None:
        print_progress("PASS participant_sensor_tfr_atlas_v1 reused validated unchanged output")
        return existing
    if args.verify_current:
        raise RuntimeError("atlas_verify_current_requires_matching_complete_output")
    if final_root.exists():
        if not args.replace_current:
            raise RuntimeError("atlas_existing_output_is_stale_or_invalid_refusing_to_overwrite")
        archived = archive_current_rendering(final_root)
        print_progress(f"archived prior rendering={archived.name}")

    lineage = load_lineage(repo_root, tfr)
    selection = select_assigned_condition_trials(lineage, config)
    pages = participant_page_table(selection)
    positions = load_layout_positions(repo_root / config["paths"]["historical_layout"], tfr["channels"])
    temporary_root = final_root.parent / f".{final_root.name}.tmp-{os.getpid()}-{time.time_ns()}"
    temporary_root.mkdir(parents=True, exist_ok=False)
    try:
        atomic_write_json(temporary_root / RUN_STATE, {
            "schema_version": ATLAS_SCHEMA_VERSION, "atlas_namespace": ATLAS_NAMESPACE,
            "status": "running", "started_at": utc_now(), "authority": authority,
        })
        (temporary_root / "configuration").mkdir(parents=True, exist_ok=True)
        shutil.copyfile(repo_root / CONFIG_PATH, temporary_root / "configuration/participant_tfr_atlas_config_v1.yaml")
        shutil.copyfile(repo_root / config["paths"]["historical_layout"], temporary_root / "configuration/historical_besa_layout_coordinates.csv")
        pdf_paths = {group: temporary_root / filename for group, filename in config["display"]["pdfs"].items()}
        page_rows: list[dict[str, Any]] = []
        saturation_rows: list[dict[str, Any]] = []
        selected = selection.loc[selection["atlas_selected_for_mean"]].copy()
        with PdfPages(pdf_paths["overt"]) as overt_pdf, PdfPages(pdf_paths["imagery"]) as imagery_pdf:
            writers = {"overt": overt_pdf, "imagery": imagery_pdf}
            for page in pages.itertuples(index=False):
                participant = selected.loc[selected["behavioural_id"].eq(page.behavioural_id)].copy()
                recording = tfr["recordings"][page.recording_stem]
                rows = participant["tfr_row_index"].to_numpy(dtype=int)
                onset = mean_memory_mapped_trials(recording["arrays"]["response_onset_db"], rows, (30, 37, 201))
                end = mean_memory_mapped_trials(recording["arrays"]["response_end_db"], rows, (30, 37, 201))
                saturation = saturation_accounting(onset, end, -8.0, 8.0)
                metadata = participant_metadata(participant, saturation)
                figure = render_participant_page(
                    onset=onset, end=end, metadata=metadata, channels=tfr["channels"],
                    frequencies=tfr["frequencies"], times=tfr["times"], positions=positions, config=config,
                )
                writers[page.pdf_group].savefig(figure)
                plt.close(figure)
                page_rows.append({
                    **metadata, "pdf_group": page.pdf_group,
                    "pdf_filename": config["display"]["pdfs"][page.pdf_group],
                    "page_number": int(page.page_number), "physical_unmirrored_positions": True,
                    "layout_coordinate_role": "historical_besa_derived_presentation_only",
                })
                saturation_rows.append({
                    "behavioural_id": metadata["behavioural_id"], "pdf_group": page.pdf_group,
                    "page_number": int(page.page_number), **saturation,
                })
                print_progress(f"rendered participant={page.behavioural_id} pdf={page.pdf_group} page={page.page_number}")
        page_manifest = pd.DataFrame(page_rows).sort_values(["pdf_group", "page_number"]).reset_index(drop=True)
        saturation_frame = pd.DataFrame(saturation_rows).sort_values(["pdf_group", "page_number"]).reset_index(drop=True)
        if len(page_manifest) != EXPECTED_PARTICIPANTS or page_manifest["behavioural_id"].duplicated().any():
            raise RuntimeError("atlas_rendered_page_identity_mismatch")
        atomic_write_csv(temporary_root / "participant_page_manifest.csv", page_manifest)
        atomic_write_csv(temporary_root / "summaries/participant_trial_and_saturation.csv", saturation_frame)
        atomic_write_csv(temporary_root / "summaries/assigned_condition_trial_selection.csv", selection.loc[:, [
            "canonical_event_key", "behavioural_id", "recording_stem", "group", "performed_condition",
            "condition_semantics", "atlas_pdf_group", "atlas_assigned_condition", "atlas_selected_for_mean",
            "atlas_bridge_excluded_from_mean",
        ]])

        pdf_counts = {group: pdf_page_count(path) for group, path in pdf_paths.items()}
        after = source_snapshot(sources, root=repo_root)
        compare_source_snapshots(before, after)
        saturation_total = int(saturation_frame["saturated_cell_count"].sum())
        rendered_cells = int(saturation_frame["rendered_source_cell_count"].sum())
        validation = {
            "status": "pass",
            "participant_page_count": len(page_manifest),
            "one_page_per_accepted_participant": not page_manifest["behavioural_id"].duplicated().any(),
            "pdf_page_counts": pdf_counts,
            "correct_overt_imagery_assignment": (
                int((page_manifest["pdf_group"] == "overt").sum()) == EXPECTED_OVERT_PARTICIPANTS
                and int((page_manifest["pdf_group"] == "imagery").sum()) == EXPECTED_IMAGERY_PARTICIPANTS
            ),
            "contributing_assigned_condition_trials": int(page_manifest["contributing_trial_count"].sum()),
            "imagery_final_overt_bridge_trials_excluded_from_averaging": int(selection["atlas_bridge_excluded_from_mean"].sum()),
            "all_30_sensors_once_per_page": len(tfr["channels"]) == 30 and len(set(tfr["channels"])) == 30,
            "onset_end_axes_match_tfr_v1": bool(np.array_equal(tfr["times"], np.arange(-50, 151) / 100.0)),
            "frequency_axis_matches_tfr_v1": bool(np.array_equal(tfr["frequencies"], np.arange(4.0, 41.0))),
            "fixed_db_limits": {"minimum": -8.0, "maximum": 8.0, "palette": "plasma"},
            "saturation_reproduced_from_rendered_source_cells": (
                saturation_total == int(saturation_frame["below_minimum_count"].sum() + saturation_frame["above_maximum_count"].sum())
            ),
            "saturation": {"cells": saturation_total, "total_cells": rendered_cells, "percent": 100.0 * saturation_total / rendered_cells},
            "participant_70": page_manifest.loc[page_manifest["behavioural_id"].eq(70), [
                "original_handedness", "analysis_hand", "mapping_source", "physical_unmirrored_positions"
            ]].to_dict(orient="records"),
            "participant_4_included": bool(page_manifest["behavioural_id"].eq(4).any()),
            "participant_60_included": bool(page_manifest["behavioural_id"].eq(60).any()),
            "file49_included": bool(page_manifest["recording_stem"].eq("demi_49_data").any()),
            "file54_1_absent": not page_manifest["recording_stem"].str.contains("54_1", regex=False).any(),
            "participant_86_absent": not page_manifest["behavioural_id"].eq(86).any(),
            "current_surface_reference": "average_referenced_tfr_v1_not_historical_csd",
            "source_tfr_and_downstream_products_unchanged": True,
        }
        expected_true = (
            validation["participant_page_count"] == EXPECTED_PARTICIPANTS,
            validation["one_page_per_accepted_participant"], validation["correct_overt_imagery_assignment"],
            validation["contributing_assigned_condition_trials"] == EXPECTED_ASSIGNED_TRIALS,
            validation["imagery_final_overt_bridge_trials_excluded_from_averaging"] == EXPECTED_BRIDGE_EXCLUDED,
            validation["all_30_sensors_once_per_page"], validation["onset_end_axes_match_tfr_v1"],
            validation["frequency_axis_matches_tfr_v1"], validation["saturation_reproduced_from_rendered_source_cells"],
            validation["participant_70"] == [{
                "original_handedness": "a", "analysis_hand": "right",
                "mapping_source": "owner_recollection", "physical_unmirrored_positions": True,
            }],
            validation["participant_4_included"], validation["participant_60_included"], validation["file49_included"],
            validation["file54_1_absent"], validation["participant_86_absent"],
            validation["source_tfr_and_downstream_products_unchanged"],
        )
        if not all(expected_true) or pdf_counts != {"overt": 41, "imagery": 40}:
            raise RuntimeError(f"atlas_validation_failed:{validation}")
        atomic_write_json(temporary_root / "validation/atlas_validation.json", validation)
        provenance = {
            "atlas_namespace": ATLAS_NAMESPACE,
            "source_tfr_v1_run_manifest_sha256": tfr["run_sha256"],
            "source_arrays_read": "response_onset_db and response_end_db only",
            "array_access": "numpy mmap_mode='r'; one participant at a time",
            "normalization": "existing trial-matched red_on -0.5 to -0.2 s dB values averaged as stored",
            "reference_and_csd_caveat": config["presentation"]["reference_statement"],
            "layout_caveat": config["presentation"]["layout_coordinate_statement"],
            "historical_audit": config["paths"]["historical_audit"],
            "selected_condition_rule": config["selection"]["selection_description"],
        }
        atomic_write_json(temporary_root / "source_provenance_manifest.json", provenance)
        immutability = {
            "status": "unchanged", "source_count": len(before),
            "source_size_mtime_fingerprint": stable_json_hash(before), "before": before, "after": after,
        }
        atomic_write_json(temporary_root / "source_immutability_evidence.json", immutability)
        summary = {
            "status": "complete", "completed_at": utc_now(), "participant_pages": EXPECTED_PARTICIPANTS,
            "assigned_condition_trial_count": EXPECTED_ASSIGNED_TRIALS,
            "imagery_bridge_trials_excluded_from_means": EXPECTED_BRIDGE_EXCLUDED,
            "pdf_page_counts": pdf_counts, "saturation": validation["saturation"],
            "runtime_seconds": time.perf_counter() - started,
        }
        atomic_write_json(temporary_root / "atlas_summary.json", summary)
        atomic_write_text(temporary_root / "atlas_summary.md", write_summary_markdown(summary))
        atomic_write_text(temporary_root / "logs/atlas_run.log", "\n".join([
            f"started={utc_now()}", "status=complete", f"participant_pages={EXPECTED_PARTICIPANTS}",
            f"assigned_condition_trials={EXPECTED_ASSIGNED_TRIALS}",
            f"saturated_rendered_mean_cells={saturation_total}",
        ]) + "\n")
        output_files = [
            pdf_paths["overt"], pdf_paths["imagery"], temporary_root / "participant_page_manifest.csv",
            temporary_root / "summaries/participant_trial_and_saturation.csv",
            temporary_root / "summaries/assigned_condition_trial_selection.csv",
            temporary_root / "configuration/participant_tfr_atlas_config_v1.yaml",
            temporary_root / "configuration/historical_besa_layout_coordinates.csv",
            temporary_root / "source_provenance_manifest.json", temporary_root / "source_immutability_evidence.json",
            temporary_root / "validation/atlas_validation.json", temporary_root / "atlas_summary.json",
            temporary_root / "atlas_summary.md", temporary_root / "logs/atlas_run.log",
        ]
        descriptors = {path.relative_to(temporary_root).as_posix(): output_descriptor(path, root=temporary_root) for path in output_files}
        manifest = {
            "schema_version": ATLAS_SCHEMA_VERSION, "atlas_namespace": ATLAS_NAMESPACE, "status": "complete",
            "started_at": utc_now(), "completed_at": utc_now(), "authority": authority,
            "pdfs": {group: output_descriptor(path, root=temporary_root) for group, path in pdf_paths.items()},
            "outputs": descriptors, "summary": summary,
            "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
        }
        atomic_write_json(temporary_root / RUN_MANIFEST, manifest)
        atomic_write_json(temporary_root / RUN_STATE, {
            "schema_version": ATLAS_SCHEMA_VERSION, "atlas_namespace": ATLAS_NAMESPACE, "status": "complete",
            "completed_at": utc_now(), "elapsed_seconds": summary["runtime_seconds"],
        })
        os.replace(temporary_root, final_root)
        print_progress(f"PASS participant_sensor_tfr_atlas_v1 pages=81 trials=8060 runtime={summary['runtime_seconds']:.1f}s")
        return manifest
    except Exception:
        # Keep incomplete output separate from the accepted destination for inspection.
        if temporary_root.exists():
            failed = final_root.parent / f"{temporary_root.name}.failed"
            os.replace(temporary_root, failed)
        raise


def main() -> None:
    """Run atlas construction, preserving the original exception on failure."""

    run(parse_args())


if __name__ == "__main__":
    main()
