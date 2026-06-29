"""Create a local manifest of raw DEMI EEG EDF files.

This script inspects the raw EDF files stored under ``_Data/eeg/raw/`` and
writes local manifest tables under ``_Data/eeg/manifest/``. The manifest exists
to give Tony a transparent first-pass inventory of observable raw-file facts
before any new MNE preprocessing or event reconstruction is attempted.

The script records facts that can be measured directly from the EDF filename or
from MNE's header/annotation view of the EDF file:

- participant ID parsed from the filename;
- source filename and local path;
- whether the file appears to be a standard single EDF, a split-file part, or
  an already concatenated EDF;
- sampling frequency, channel count, channel names, channel-type labels, and
  approximate recording duration;
- annotation counts and per-description annotation counts;
- compact historical-note categories recovered in the private inventory
  reports, kept separate from the raw EDF facts.

The script does not preload the EEG signal, preprocess the recording, repair
events, concatenate split files, write BIDS data, modify raw files, or create
new participant-retention rules. Historical categories are reported as
historical notes only. They are not treated as new MNE decisions.

Expected inputs:

- raw EDF files in ``_Data/eeg/raw/``;
- filenames that begin with ``demi_<ID>``, optionally followed by a split-part
  number such as ``demi_54_1 Data.edf``;
- an installed Python environment with ``mne`` and ``pandas`` available.

Generated local outputs:

- ``_Data/eeg/manifest/raw_eeg_file_manifest.csv``: one row per EDF file;
- ``_Data/eeg/manifest/raw_eeg_annotation_counts.csv``: one row per EDF file
  and annotation description;
- ``_Data/eeg/manifest/raw_eeg_manifest_summary.md``: a short local summary of
  file counts, annotation-count flags, split/concatenated files, and read
  errors.

Safety boundaries:

- EDFs are opened with ``mne.io.read_raw_edf(..., preload=False)``.
- No raw signal array is loaded into memory.
- No raw file is edited.
- One unreadable EDF produces an error row and does not stop the manifest run.
- The generated manifest outputs live under ``_Data/`` and should remain
  local-only.

Run from the repository root with:

    python analysis/eeg_mne/00_make_raw_eeg_manifest.py

The script also works when launched from another directory because paths are
resolved relative to this file's location in the repository.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import mne
import pandas as pd


# The old Python pipeline used fewer than 600 annotations as a completeness
# warning. Here the value is used only as a factual low-count flag, not as a
# participant-retention rule.
LOW_ANNOTATION_COUNT_THRESHOLD = 600

# These IDs were named in the private decision-source report as cases needing
# special attention or conflict tracking. The flag created from this set is a
# manifest inspection flag only.
KNOWN_CONFLICT_OR_SPECIAL_CASE_IDS = {5, 11, 54, 56, 65, 86, 88, 99}

# Historical behavioural/session categories recovered from old R-script
# comments and summarized in the private inventory reports. These labels are
# deliberately compact so they do not copy private free-text notes into code.
HISTORICAL_BEHAVIOURAL_SESSION_CATEGORIES = {
    9: "experimenter_error_or_wrong_protocol",
    10: "experimenter_error_or_wrong_protocol",
    12: "experimenter_error_or_wrong_protocol",
    14: "experimenter_error_or_wrong_protocol",
    20: "experimenter_error_or_wrong_protocol",
    89: "technical_noncompletion",
    96: "technical_noncompletion",
    100: "technical_noncompletion",
    13: "id_skipped",
    26: "id_skipped",
    36: "id_skipped",
    78: "id_skipped",
}

# Historical EEG categories recovered from the old dissertation EEG path and
# summarized in the private decision-source report. These are historical notes,
# not current MNE decisions.
HISTORICAL_EEG_CATEGORIES = {
    6: "eeg_not_recorded",
    7: "eeg_not_recorded",
    11: "eeg_not_recorded",
    24: "very_bad_eeg_note",
    54: "triggers_missing",
    56: "triggers_missing",
    65: "software_crash",
    86: "veol_or_veo_l_did_not_record",
}

# Special handling categories summarize historical workflow limitations or
# source conflicts that should be visible during raw EDF review.
SPECIAL_HANDLING_CATEGORIES = {
    5: "concatenated_or_duplicate_file_review",
    11: "legacy_note_conflicts_with_local_raw_edf",
    54: "split_raw_edf_trigger_review",
    56: "split_raw_edf_trigger_review",
    65: "split_raw_edf_software_crash_note",
    86: "auxiliary_channel_recording_note",
    88: "incomplete_task_marker_requires_review",
    99: "incomplete_task_marker_requires_review",
}

# These IDs were explicitly described as conflicting or ambiguous in the
# private reports. This field keeps ambiguity separate from measured EDF facts.
CONFLICT_OR_AMBIGUITY_FLAGS = {
    11: "legacy_eeg_not_recorded_note_conflicts_with_local_raw_edf",
    54: "split_raw_edf_and_historical_trigger_note",
    56: "split_raw_edf_and_historical_trigger_note",
    65: "split_raw_edf_and_historical_crash_note",
    88: "incomplete_task_marker_not_resolved_by_raw_edf_manifest",
    99: "incomplete_task_marker_not_resolved_by_raw_edf_manifest",
}

# Name-based channel labels borrowed from the old EDF-to-BIDS helper. MNE reads
# many EDF channels as EEG by default, so these labels make auxiliary channels
# easier to audit without mutating the Raw object.
CHANNEL_TYPE_OVERRIDES = {
    "HEO": "eog",
    "VEO": "eog",
    "VEOL": "eog",
    "VEO-L": "eog",
    "EMG-L": "emg",
    "EMG-A": "emg",
}

EDF_FILENAME_RE = re.compile(
    r"^demi_(?P<participant_id>\d{1,3})(?:_(?P<split_part>\d+))?(?P<label>.*)\.edf$",
    flags=re.IGNORECASE,
)


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Inputs/arguments:
        None.

    Return value:
        A ``Path`` pointing to the repository root. The script is expected to
        live at ``analysis/eeg_mne/00_make_raw_eeg_manifest.py``.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def discover_edf_files(raw_dir: Path) -> list[Path]:
    """Find EDF files in the raw EEG directory.

    Inputs/arguments:
        raw_dir: Directory expected to contain raw DEMI EDF files.

    Return value:
        A sorted list of EDF paths directly inside ``raw_dir``. The search is
        case-insensitive for the ``.edf`` extension and is not recursive.

    Side effects:
        None.
    """

    if not raw_dir.exists():
        return []

    return sorted(path for path in raw_dir.iterdir() if path.is_file() and path.suffix.lower() == ".edf")


def parse_edf_filename(edf_path: Path) -> dict[str, Any]:
    """Parse DEMI participant and file-role information from an EDF filename.

    Inputs/arguments:
        edf_path: Path to an EDF file. Only ``edf_path.name`` is parsed.

    Return value:
        A dictionary with the parsed participant ID, zero-padded participant
        label, file role, split-part number if present, and a parse-warning
        string. If the filename cannot be parsed, participant fields are null
        and the role is ``unknown``.

    Side effects:
        None.
    """

    match = EDF_FILENAME_RE.match(edf_path.name)
    if match is None:
        return {
            "participant_id": None,
            "participant_id_padded": "",
            "file_role": "unknown",
            "split_part": None,
            "filename_parse_warning": "filename_does_not_match_expected_demi_pattern",
        }

    participant_id = int(match.group("participant_id"))
    split_part_text = match.group("split_part")
    split_part = int(split_part_text) if split_part_text is not None else None
    label_text = match.group("label").lower()

    # Split-part filenames carry an extra numeric token after the participant
    # ID, for example "demi_54_1 Data.edf". Concatenated files are labelled by
    # text in the filename. Plain files such as "demi_27.edf" are treated as
    # single raw EDFs because the participant ID still parses cleanly.
    if split_part is not None:
        file_role = "split_part"
    elif "concatenated" in label_text:
        file_role = "concatenated"
    else:
        file_role = "single"

    return {
        "participant_id": participant_id,
        "participant_id_padded": f"{participant_id:03d}",
        "file_role": file_role,
        "split_part": split_part,
        "filename_parse_warning": "",
    }


def historical_notes_for_participant(participant_id: int | None) -> dict[str, Any]:
    """Return compact historical-note fields for one participant ID.

    Inputs/arguments:
        participant_id: Integer DEMI participant ID, or ``None`` when a
        filename could not be parsed.

    Return value:
        A dictionary containing historical behavioural/session category,
        historical EEG category, special handling category, and
        conflict/ambiguity fields. Empty strings mean no category was recovered
        from the private reports for that ID.

    Side effects:
        None.
    """

    if participant_id is None:
        return {
            "historical_behavioural_session_exclusion_category": "",
            "historical_eeg_exclusion_category": "",
            "special_handling_category": "",
            "conflict_ambiguity_flag": False,
            "conflict_ambiguity_note": "",
            "known_conflict_or_special_case_id": False,
        }

    conflict_note = CONFLICT_OR_AMBIGUITY_FLAGS.get(participant_id, "")
    return {
        "historical_behavioural_session_exclusion_category": HISTORICAL_BEHAVIOURAL_SESSION_CATEGORIES.get(
            participant_id, ""
        ),
        "historical_eeg_exclusion_category": HISTORICAL_EEG_CATEGORIES.get(participant_id, ""),
        "special_handling_category": SPECIAL_HANDLING_CATEGORIES.get(participant_id, ""),
        "conflict_ambiguity_flag": bool(conflict_note),
        "conflict_ambiguity_note": conflict_note,
        "known_conflict_or_special_case_id": participant_id in KNOWN_CONFLICT_OR_SPECIAL_CASE_IDS,
    }


def infer_channel_types(channel_names: list[str], mne_channel_types: list[str]) -> list[str]:
    """Create auditable channel-type labels for manifest output.

    Inputs/arguments:
        channel_names: Channel names in the order returned by MNE.
        mne_channel_types: Channel-type labels returned by
        ``raw.get_channel_types()`` in the same order.

    Return value:
        A list of channel-type labels. Name-based overrides are applied for
        known auxiliary channels such as HEO, VEO, and EMG channels; all other
        labels keep MNE's reported type.

    Side effects:
        None.
    """

    inferred_types: list[str] = []
    for channel_name, mne_type in zip(channel_names, mne_channel_types):
        # EDF readers often lack rich channel typing. Use conservative,
        # explicit name matches for known auxiliary channels and otherwise keep
        # MNE's label untouched.
        normalized_name = channel_name.strip().upper()
        inferred_types.append(CHANNEL_TYPE_OVERRIDES.get(normalized_name, mne_type))

    return inferred_types


def json_dump(value: Any) -> str:
    """Serialize a manifest value into stable JSON text.

    Inputs/arguments:
        value: A JSON-serializable Python value.

    Return value:
        A JSON string with deterministic key ordering where applicable.

    Side effects:
        None.
    """

    return json.dumps(value, sort_keys=True)


def summarize_annotations(raw: mne.io.BaseRaw) -> tuple[int, dict[str, int]]:
    """Summarize EDF annotations already attached to an MNE Raw object.

    Inputs/arguments:
        raw: An MNE Raw object opened from an EDF file with ``preload=False``.

    Return value:
        A tuple of ``(annotation_count, description_counts)`` where
        ``description_counts`` maps each annotation description to its count.

    Side effects:
        None.
    """

    descriptions = [str(description) for description in raw.annotations.description]
    return len(descriptions), dict(sorted(Counter(descriptions).items()))


def extract_mne_metadata(edf_path: Path) -> dict[str, Any]:
    """Read observable EDF metadata through MNE without preloading data.

    Inputs/arguments:
        edf_path: Path to one raw EDF file.

    Return value:
        A dictionary containing read status, sampling frequency, channel count,
        channel names, channel-type labels, recording duration, and annotation
        summaries. If MNE raises an error, the dictionary records
        ``read_status='mne_read_error'`` and keeps metadata fields empty.

    Side effects:
        Opens the EDF through MNE. The file handle is closed before returning
        when the Raw object exposes ``close()``.
    """

    raw: mne.io.BaseRaw | None = None
    try:
        raw = mne.io.read_raw_edf(edf_path, preload=False, verbose="ERROR")

        channel_names = list(raw.ch_names)
        mne_channel_types = raw.get_channel_types()
        inferred_channel_types = infer_channel_types(channel_names, mne_channel_types)
        annotation_count, annotation_description_counts = summarize_annotations(raw)

        # ``raw.n_times`` and ``info['sfreq']`` are header-derived values here;
        # no signal samples are loaded because preload=False was used above.
        sfreq = float(raw.info["sfreq"])
        duration_seconds = float(raw.n_times / sfreq) if sfreq else None

        return {
            "read_status": "ok",
            "mne_read_error_type": "",
            "mne_read_error_message": "",
            "sampling_frequency_hz": sfreq,
            "n_channels": len(channel_names),
            "channel_names_json": json_dump(channel_names),
            "channel_types_json": json_dump(inferred_channel_types),
            "channel_type_counts_json": json_dump(dict(sorted(Counter(inferred_channel_types).items()))),
            "duration_seconds": duration_seconds,
            "annotation_count": annotation_count,
            "annotation_description_counts_json": json_dump(annotation_description_counts),
        }
    except Exception as error:  # noqa: BLE001 - a bad EDF should not stop the manifest.
        # Keep a row for unreadable EDFs so the manifest shows which file needs
        # manual follow-up. The exact exception class and message are preserved.
        return {
            "read_status": "mne_read_error",
            "mne_read_error_type": type(error).__name__,
            "mne_read_error_message": str(error),
            "sampling_frequency_hz": None,
            "n_channels": None,
            "channel_names_json": "[]",
            "channel_types_json": "[]",
            "channel_type_counts_json": "{}",
            "duration_seconds": None,
            "annotation_count": None,
            "annotation_description_counts_json": "{}",
        }
    finally:
        if raw is not None and hasattr(raw, "close"):
            raw.close()


def inspection_issue_codes(row: dict[str, Any]) -> list[str]:
    """Create factual inspection flags for a manifest row.

    Inputs/arguments:
        row: A manifest-row dictionary after filename parsing, historical-note
        lookup, and MNE metadata extraction.

    Return value:
        A list of short issue codes. Codes describe observable read/filename/
        annotation facts or known historical-conflict IDs. They are not
        participant-retention rules.

    Side effects:
        None.
    """

    issue_codes: list[str] = []

    if row["read_status"] == "mne_read_error":
        issue_codes.append("mne_read_error")

    annotation_count = row.get("annotation_count")
    if annotation_count == 0:
        issue_codes.append("zero_annotations")
    elif annotation_count is not None and annotation_count < LOW_ANNOTATION_COUNT_THRESHOLD:
        issue_codes.append("low_annotation_count")

    if row["file_role"] == "split_part":
        issue_codes.append("split_file")
    elif row["file_role"] == "concatenated":
        issue_codes.append("concatenated_file")

    if row["known_conflict_or_special_case_id"]:
        issue_codes.append("known_conflict_or_special_case_id")

    if row["filename_parse_warning"]:
        issue_codes.append("filename_parse_warning")

    return issue_codes


def manifest_row_for_file(edf_path: Path, repo_root: Path) -> dict[str, Any]:
    """Build one manifest row for one EDF file.

    Inputs/arguments:
        edf_path: Path to the raw EDF file.
        repo_root: Repository root used to store a readable relative path.

    Return value:
        A dictionary containing filename-derived fields, historical-note
        fields, MNE metadata fields, and semicolon-separated inspection issues.

    Side effects:
        Opens the EDF through MNE via ``extract_mne_metadata``. Does not write
        files or modify the EDF.
    """

    parsed = parse_edf_filename(edf_path)
    historical_notes = historical_notes_for_participant(parsed["participant_id"])
    mne_metadata = extract_mne_metadata(edf_path)

    row = {
        **parsed,
        "source_filename": edf_path.name,
        "file_path": edf_path.relative_to(repo_root).as_posix(),
        **historical_notes,
        **mne_metadata,
    }
    row["inspection_issue_codes"] = ";".join(inspection_issue_codes(row))
    return row


def build_file_manifest(edf_files: list[Path], repo_root: Path) -> pd.DataFrame:
    """Build the one-row-per-EDF manifest table.

    Inputs/arguments:
        edf_files: Sorted list of raw EDF paths.
        repo_root: Repository root used for relative file paths.

    Return value:
        A pandas ``DataFrame`` with one row per EDF file.

    Side effects:
        Opens each EDF through MNE without preloading signal data.
    """

    rows = [manifest_row_for_file(edf_path, repo_root) for edf_path in edf_files]
    return pd.DataFrame(rows)


def build_annotation_counts(manifest_df: pd.DataFrame) -> pd.DataFrame:
    """Build a long annotation-count table from the file manifest.

    Inputs/arguments:
        manifest_df: The one-row-per-EDF manifest table.

    Return value:
        A pandas ``DataFrame`` with one row per EDF file and annotation
        description. Files with no annotations or read errors are represented by
        a row with an empty annotation description and count 0.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []

    for _, row in manifest_df.iterrows():
        # The manifest stores nested annotation counts as JSON so the main CSV
        # remains one row per EDF file. Decode it here for a tidy count table.
        description_counts = json.loads(row["annotation_description_counts_json"])
        if not description_counts:
            rows.append(
                {
                    "participant_id": row["participant_id"],
                    "participant_id_padded": row["participant_id_padded"],
                    "source_filename": row["source_filename"],
                    "file_role": row["file_role"],
                    "split_part": row["split_part"],
                    "annotation_description": "",
                    "annotation_count": 0,
                }
            )
            continue

        for description, count in description_counts.items():
            rows.append(
                {
                    "participant_id": row["participant_id"],
                    "participant_id_padded": row["participant_id_padded"],
                    "source_filename": row["source_filename"],
                    "file_role": row["file_role"],
                    "split_part": row["split_part"],
                    "annotation_description": description,
                    "annotation_count": count,
                }
            )

    return pd.DataFrame(rows)


def write_manifest_summary(manifest_df: pd.DataFrame, output_path: Path, raw_dir: Path) -> None:
    """Write a short Markdown summary of the raw EDF manifest.

    Inputs/arguments:
        manifest_df: The one-row-per-EDF manifest table.
        output_path: Destination Markdown path.
        raw_dir: Raw EDF directory summarized by this run.

    Return value:
        ``None``.

    Side effects:
        Writes ``output_path``.
    """

    issue_counter: Counter[str] = Counter()
    for issue_text in manifest_df["inspection_issue_codes"].fillna(""):
        if not issue_text:
            continue
        issue_counter.update(issue_text.split(";"))

    read_error_df = manifest_df[manifest_df["read_status"] == "mne_read_error"]
    special_df = manifest_df[manifest_df["known_conflict_or_special_case_id"]]

    role_counts = manifest_df["file_role"].value_counts(dropna=False).sort_index()
    status_counts = manifest_df["read_status"].value_counts(dropna=False).sort_index()

    lines = [
        "# Raw EEG File Manifest Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Raw EDF directory: `{raw_dir.as_posix()}`",
        f"EDF files inspected: {len(manifest_df)}",
        f"Unique parsed participant IDs: {manifest_df['participant_id'].dropna().nunique()}",
        "",
        "## File roles",
        "",
    ]

    for role, count in role_counts.items():
        lines.append(f"- {role}: {count}")

    lines.extend(["", "## MNE read status", ""])
    for status, count in status_counts.items():
        lines.append(f"- {status}: {count}")

    lines.extend(["", "## Inspection issue counts", ""])
    if issue_counter:
        for issue, count in sorted(issue_counter.items()):
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- none")

    lines.extend(["", "## Known conflict or special-case IDs present", ""])
    if special_df.empty:
        lines.append("- none")
    else:
        for _, row in special_df.sort_values(["participant_id", "source_filename"]).iterrows():
            note_parts = [
                part
                for part in [
                    row["historical_behavioural_session_exclusion_category"],
                    row["historical_eeg_exclusion_category"],
                    row["special_handling_category"],
                    row["conflict_ambiguity_note"],
                ]
                if part
            ]
            note_text = "; ".join(note_parts) if note_parts else "historical category present"
            lines.append(f"- {int(row['participant_id'])}: `{row['source_filename']}` - {note_text}")

    lines.extend(["", "## MNE read errors", ""])
    if read_error_df.empty:
        lines.append("- none")
    else:
        for _, row in read_error_df.sort_values("source_filename").iterrows():
            lines.append(
                f"- `{row['source_filename']}`: {row['mne_read_error_type']} - {row['mne_read_error_message']}"
            )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Annotation-count flags are factual count flags only.",
            "- Historical fields are recovered notes, not new MNE decisions.",
            "- This script does not preload signal data and does not alter EDF files.",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(manifest_df: pd.DataFrame, annotation_counts_df: pd.DataFrame, output_dir: Path, raw_dir: Path) -> None:
    """Write manifest CSVs and the local Markdown summary.

    Inputs/arguments:
        manifest_df: One-row-per-EDF manifest table.
        annotation_counts_df: Long annotation-count table.
        output_dir: Directory where manifest outputs should be written.
        raw_dir: Raw EDF directory summarized by this run.

    Return value:
        ``None``.

    Side effects:
        Creates ``output_dir`` if needed and writes three local-only manifest
        files below it.
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = output_dir / "raw_eeg_file_manifest.csv"
    annotation_counts_path = output_dir / "raw_eeg_annotation_counts.csv"
    summary_path = output_dir / "raw_eeg_manifest_summary.md"

    manifest_df.to_csv(manifest_path, index=False)
    annotation_counts_df.to_csv(annotation_counts_path, index=False)
    write_manifest_summary(manifest_df, summary_path, raw_dir)


def main() -> None:
    """Run the raw EEG manifest workflow.

    Inputs/arguments:
        None. Paths are resolved relative to the repository root.

    Return value:
        ``None``.

    Side effects:
        Reads EDF headers and annotations from ``_Data/eeg/raw/`` without
        preloading signal data, then writes local outputs under
        ``_Data/eeg/manifest/``.
    """

    repo_root = repo_root_from_script()
    raw_dir = repo_root / "_Data" / "eeg" / "raw"
    output_dir = repo_root / "_Data" / "eeg" / "manifest"

    edf_files = discover_edf_files(raw_dir)
    manifest_df = build_file_manifest(edf_files, repo_root)
    annotation_counts_df = build_annotation_counts(manifest_df)
    write_outputs(manifest_df, annotation_counts_df, output_dir, raw_dir)

    print(f"Inspected {len(manifest_df)} EDF files from {raw_dir}")
    print(f"Wrote manifest outputs to {output_dir}")


if __name__ == "__main__":
    main()
