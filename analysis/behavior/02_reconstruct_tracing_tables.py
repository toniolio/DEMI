"""Reconstruct DEMI TraceLab figure and tracing tables from zip contents.

This script opens local TraceLab figure/tracing zip files below
``_Data/figure/`` and reconstructs the first old-style intermediate tables used
by the historical DEMI behavioural/tracing pipeline. It is part of the EEG
reanalysis infrastructure because later EEG event reconstruction depends on
knowing the participant/session/block/trial structure of task rows, TraceLab
figures, physical tracings, and their timing offsets before raw EEG annotations
are compared to those behavioural/tracing records.

The old R import script, ``_Scripts/00_import.R``, already solved the basic
TraceLab zip-content import. It recursively found figure zip files, set aside
basenames containing ``learned``, opened each regular zip, read four expected
internal TraceLab text files, parsed participant/session/block/trial/date from
the zip basename, and created table-like objects called ``points``,
``segments``, ``frames``, and ``tracings``. This Python script deliberately
ports that narrow import step in a more explicit and auditable form.

Old logic ported here:

- recursive discovery of ``.zip`` files below ``_Data/figure/``;
- zip-basename parsing for participant ID, session, block, trial, and date;
- read-only inspection of expected TraceLab text members:
  ``.tlfp`` for figure vertex points, ``.tlfs`` for figure segments,
  ``.tlf`` for animation frames, and ``.tlt`` for physical tracing samples;
- reconstruction of separate points, segments, frames, and tracings tables;
- manifest recording for missing or malformed zip/internal-file contents.

Old logic intentionally deferred:

- task-to-tracing joins from ``_Scripts/01_preprocessing.R``;
- tracing filters for no-shape trials, failed trial ends, glitches, false
  starts, hand noise, incomplete tracings, excessive gaps, or edge hits;
- movement-time summaries and accuracy/error metrics;
- behavioural trial cleanup;
- EEG event-offset construction from ``_Scripts/04_import_eeg.R``;
- EEG event updating from ``_Scripts/_functions/eeg.R``;
- raw EEG annotation inspection or EEG preprocessing.

Expected inputs:

- ``_Data/figure/`` containing local TraceLab figure/tracing zip files;
- zip basenames like ``p100_s1_b1_t10_2019-09-16.zip``;
- internal TraceLab text files named from the zip stem with extensions
  ``.tlfp``, ``.tlfs``, ``.tlf``, and ``.tlt``;
- a Python environment with ``pandas`` and ``pyarrow`` available.

Generated local outputs:

- ``_Data/behavior/tracing_tables/points.parquet``;
- ``_Data/behavior/tracing_tables/segments.parquet``;
- ``_Data/behavior/tracing_tables/frames.parquet``;
- ``_Data/behavior/tracing_tables/tracings.parquet``;
- ``_Data/behavior/tracing_tables/zip_parse_manifest.csv``;
- ``_Data/behavior/tracing_tables/tracing_reconstruction_summary.md``.

Safety boundaries:

- source zip files are opened read-only through Python's ``zipfile`` module;
- source data are never modified;
- no TraceLab filter or behavioural outcome is computed here;
- one malformed zip or internal file is recorded in the manifest and does not
  stop reconstruction for other zip files;
- outputs are written only below ``_Data/behavior/tracing_tables/`` and should
  remain local-only.

Run from the repository root with:

    PATH="$(pwd)/.venv/bin:$PATH" python3 analysis/behavior/02_reconstruct_tracing_tables.py

The script also works when launched from another directory because paths are
resolved relative to this file's location in the repository.
"""

from __future__ import annotations

import ast
import json
import re
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


EXPECTED_INTERNAL_EXTENSIONS = {
    "points": ".tlfp",
    "segments": ".tlfs",
    "frames": ".tlf",
    "tracings": ".tlt",
}

TABLE_OUTPUT_FILENAMES = {
    "points": "points.parquet",
    "segments": "segments.parquet",
    "frames": "frames.parquet",
    "tracings": "tracings.parquet",
}

BATCH_ROW_LIMIT = 50_000

FIGURE_ZIP_FILENAME_RE = re.compile(
    r"^p(?P<participant_id>\d{1,3})_s(?P<session>\d+)_b(?P<block>\d+)_t(?P<trial>\d+)_"
    r"(?P<date_token>\d{4}-\d{2}-\d{2})(?P<extra>.*)\.zip$",
    flags=re.IGNORECASE,
)
DATE_TOKEN_RE = re.compile(r"\d{4}-\d{2}-\d{2}")
TIME_TOKEN_RE = re.compile(r"(?<!\d)\d{2}:\d{2}:\d{2}(?!\d)")


POINTS_COLUMNS = [
    "participant_id",
    "session",
    "block",
    "trial",
    "zip_date_token",
    "source_zip_path",
    "source_zip_basename",
    "internal_source_file_name",
    "learned_flag",
    "point_index",
    "x",
    "y",
]

SEGMENTS_COLUMNS = [
    "participant_id",
    "session",
    "block",
    "trial",
    "zip_date_token",
    "source_zip_path",
    "source_zip_basename",
    "internal_source_file_name",
    "learned_flag",
    "segment_index",
    "start.x",
    "start.y",
    "end.x",
    "end.y",
    "ctrl.x",
    "ctrl.y",
]

TIMED_POINT_COLUMNS = [
    "participant_id",
    "session",
    "block",
    "trial",
    "zip_date_token",
    "source_zip_path",
    "source_zip_basename",
    "internal_source_file_name",
    "learned_flag",
    "sample_index",
    "x",
    "y",
    "time",
]

TABLE_COLUMNS = {
    "points": POINTS_COLUMNS,
    "segments": SEGMENTS_COLUMNS,
    "frames": TIMED_POINT_COLUMNS,
    "tracings": TIMED_POINT_COLUMNS,
}


def repo_root_from_script() -> Path:
    """Return the repository root inferred from this script location.

    Args:
        None.

    Returns:
        A ``Path`` pointing to the repository root. The script is expected to
        live at ``analysis/behavior/02_reconstruct_tracing_tables.py``.

    Side effects:
        None.
    """

    return Path(__file__).resolve().parents[2]


def require_pyarrow_modules() -> tuple[Any, Any]:
    """Import pyarrow modules needed for parquet output.

    Args:
        None.

    Returns:
        ``(pa, pq)`` where ``pa`` is ``pyarrow`` and ``pq`` is
        ``pyarrow.parquet``.

    Side effects:
        Imports modules. Raises a clear ``RuntimeError`` if parquet support is
        not available in the active Python environment.
    """

    try:
        import pyarrow as pa  # type: ignore[import-not-found]
        import pyarrow.parquet as pq  # type: ignore[import-not-found]
    except ImportError as error:
        raise RuntimeError(
            "pyarrow is required to write tracing-table parquet outputs. "
            "Run this script in the project environment with pyarrow installed."
        ) from error

    return pa, pq


def json_dump(value: Any) -> str:
    """Serialize a value into stable JSON text for manifest cells.

    Args:
        value: JSON-serializable Python value.

    Returns:
        Deterministic JSON text.

    Side effects:
        None.
    """

    return json.dumps(value, sort_keys=True)


def discover_figure_zip_files(figure_dir: Path) -> list[Path]:
    """Find TraceLab zip files below the local figure directory.

    Args:
        figure_dir: Directory expected to contain participant figure folders.

    Returns:
        A sorted list of paths whose suffix is ``.zip`` case-insensitively.

    Side effects:
        Reads directory entries only. Zip contents are not opened here.
    """

    if not figure_dir.exists():
        return []

    return sorted(path for path in figure_dir.rglob("*") if path.is_file() and path.suffix.lower() == ".zip")


def relative_path(path: Path, repo_root: Path) -> str:
    """Return a repository-relative POSIX path when possible.

    Args:
        path: Path to represent.
        repo_root: Repository root used for relative paths.

    Returns:
        A POSIX-style relative path if ``path`` is below ``repo_root``;
        otherwise an absolute POSIX-style path.

    Side effects:
        None.
    """

    try:
        return path.relative_to(repo_root).as_posix()
    except ValueError:
        return path.as_posix()


def parse_zip_filename(zip_path: Path) -> dict[str, Any]:
    """Parse participant/session/block/trial/date tokens from a zip basename.

    Args:
        zip_path: Path to one TraceLab zip file.

    Returns:
        A dictionary with parsed tokens, date/time tokens, learned flag, and a
        filename parse status.

    Side effects:
        None.
    """

    filename = zip_path.name
    match = FIGURE_ZIP_FILENAME_RE.match(filename)
    date_tokens = DATE_TOKEN_RE.findall(filename)
    time_tokens = TIME_TOKEN_RE.findall(filename)
    learned_flag = "learned" in filename.lower()

    if match is None:
        return {
            "participant_id": None,
            "session": None,
            "block": None,
            "trial": None,
            "zip_date_token": date_tokens[0] if date_tokens else "",
            "filename_date_tokens": ";".join(date_tokens),
            "filename_time_tokens": ";".join(time_tokens),
            "zip_filename_parse_status": "zip_filename_parse_error",
            "zip_filename_extra_text": "",
            "learned_flag": learned_flag,
        }

    return {
        "participant_id": int(match.group("participant_id")),
        "session": int(match.group("session")),
        "block": int(match.group("block")),
        "trial": int(match.group("trial")),
        "zip_date_token": match.group("date_token"),
        "filename_date_tokens": ";".join(date_tokens),
        "filename_time_tokens": ";".join(time_tokens),
        "zip_filename_parse_status": "ok",
        "zip_filename_extra_text": match.group("extra"),
        "learned_flag": learned_flag,
    }


def find_internal_trace_file(zip_names: list[str], zip_stem: str, extension: str) -> tuple[str | None, str]:
    """Find the internal TraceLab member corresponding to one expected suffix.

    Args:
        zip_names: Internal member names reported by ``ZipFile.namelist()``.
        zip_stem: Stem of the source zip basename.
        extension: Expected TraceLab extension such as ``.tlfp``.

    Returns:
        ``(internal_name, issue_code)``. ``internal_name`` is ``None`` when no
        suitable member is found. ``issue_code`` is empty for an exact stem
        match and describes missing or nonstandard cases otherwise.

    Side effects:
        None.
    """

    expected_basename = f"{zip_stem}{extension}"
    exact_matches = sorted(name for name in zip_names if Path(name).name == expected_basename)
    if len(exact_matches) == 1:
        return exact_matches[0], ""
    if len(exact_matches) > 1:
        return exact_matches[0], f"multiple_exact_internal_files_{extension.lstrip('.')}"

    # The old R code used unz(f, paste0(fname, extension)), which expects the
    # exact basename. This fallback keeps reconstruction moving if a zip stores
    # the expected TraceLab type under a different internal path or basename,
    # while the manifest keeps that naming irregularity visible.
    suffix_matches = sorted(name for name in zip_names if Path(name).suffix.lower() == extension)
    if len(suffix_matches) == 1:
        return suffix_matches[0], f"nonstandard_internal_filename_{extension.lstrip('.')}"
    if len(suffix_matches) > 1:
        return suffix_matches[0], f"multiple_internal_files_{extension.lstrip('.')}"

    return None, f"missing_internal_file_{extension.lstrip('.')}"


def read_internal_text(zip_handle: zipfile.ZipFile, internal_name: str) -> str:
    """Read one internal TraceLab text file from an open zip.

    Args:
        zip_handle: Open ``ZipFile`` object.
        internal_name: Internal member name to read.

    Returns:
        Decoded text, using replacement characters for unexpected byte
        sequences.

    Side effects:
        Reads one internal zip member. The source zip is not modified.
    """

    return zip_handle.read(internal_name).decode("utf-8", errors="replace").strip()


def parse_literal_sequence(text: str) -> list[Any]:
    """Parse a TraceLab tuple/list text payload into Python objects.

    Args:
        text: TraceLab text payload, usually formatted like
            ``[(x, y), ...]`` or ``[(x, y, time), ...]``.

    Returns:
        A list of parsed records. ``NA`` and blank payloads return an empty
        list because the old tracing path treated ``NA`` as no physical tracing
        samples for that trial.

    Side effects:
        None.
    """

    cleaned = text.strip()
    if cleaned in {"", "NA"}:
        return []

    parsed = ast.literal_eval(cleaned)
    if isinstance(parsed, tuple):
        return list(parsed)
    if isinstance(parsed, list):
        return parsed

    raise ValueError(f"expected list-like TraceLab payload, got {type(parsed).__name__}")


def numeric_value(value: Any) -> float:
    """Coerce a parsed TraceLab coordinate or timestamp to float.

    Args:
        value: Parsed scalar value from a TraceLab payload.

    Returns:
        Numeric value as ``float``.

    Side effects:
        None.
    """

    return float(value)


def base_output_row(metadata: dict[str, Any], zip_path: Path, repo_root: Path, internal_name: str) -> dict[str, Any]:
    """Build source metadata shared by every reconstructed output row.

    Args:
        metadata: Parsed zip filename metadata.
        zip_path: Source zip path.
        repo_root: Repository root used for relative source paths.
        internal_name: Internal TraceLab member that supplied the row.

    Returns:
        A dictionary containing participant/session/block/trial and source
        provenance columns.

    Side effects:
        None.
    """

    return {
        "participant_id": metadata["participant_id"],
        "session": metadata["session"],
        "block": metadata["block"],
        "trial": metadata["trial"],
        "zip_date_token": metadata["zip_date_token"],
        "source_zip_path": relative_path(zip_path, repo_root),
        "source_zip_basename": zip_path.name,
        "internal_source_file_name": internal_name,
        "learned_flag": bool(metadata["learned_flag"]),
    }


def parse_points_rows(text: str, metadata: dict[str, Any], zip_path: Path, repo_root: Path, internal_name: str) -> list[dict[str, Any]]:
    """Parse TraceLab ``.tlfp`` figure vertex points into rows.

    Args:
        text: Internal ``.tlfp`` text payload.
        metadata: Parsed zip filename metadata.
        zip_path: Source zip path.
        repo_root: Repository root used for relative source paths.
        internal_name: Internal member name that supplied the text.

    Returns:
        A list of rows with ``point_index``, ``x``, and ``y`` columns plus
        source metadata.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    for index, point in enumerate(parse_literal_sequence(text), start=1):
        if len(point) != 2:
            raise ValueError(f"point {index} has {len(point)} values; expected 2")
        row = base_output_row(metadata, zip_path, repo_root, internal_name)
        row.update({"point_index": index, "x": numeric_value(point[0]), "y": numeric_value(point[1])})
        rows.append(row)
    return rows


def parse_segment_rows(text: str, metadata: dict[str, Any], zip_path: Path, repo_root: Path, internal_name: str) -> list[dict[str, Any]]:
    """Parse TraceLab ``.tlfs`` figure segments into rows.

    Args:
        text: Internal ``.tlfs`` text payload.
        metadata: Parsed zip filename metadata.
        zip_path: Source zip path.
        repo_root: Repository root used for relative source paths.
        internal_name: Internal member name that supplied the text.

    Returns:
        A list of rows with start, end, and control-point coordinate columns
        plus source metadata.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    for index, segment in enumerate(parse_literal_sequence(text), start=1):
        if len(segment) != 3:
            raise ValueError(f"segment {index} has {len(segment)} points; expected 3")
        start, end, control = segment
        if len(start) != 2 or len(end) != 2 or len(control) != 2:
            raise ValueError(f"segment {index} contains a point without x/y coordinates")
        row = base_output_row(metadata, zip_path, repo_root, internal_name)
        row.update(
            {
                "segment_index": index,
                "start.x": numeric_value(start[0]),
                "start.y": numeric_value(start[1]),
                "end.x": numeric_value(end[0]),
                "end.y": numeric_value(end[1]),
                "ctrl.x": numeric_value(control[0]),
                "ctrl.y": numeric_value(control[1]),
            }
        )
        rows.append(row)
    return rows


def parse_timed_point_rows(
    text: str,
    metadata: dict[str, Any],
    zip_path: Path,
    repo_root: Path,
    internal_name: str,
) -> list[dict[str, Any]]:
    """Parse TraceLab ``.tlf`` or ``.tlt`` timed coordinate rows.

    Args:
        text: Internal TraceLab text payload.
        metadata: Parsed zip filename metadata.
        zip_path: Source zip path.
        repo_root: Repository root used for relative source paths.
        internal_name: Internal member name that supplied the text.

    Returns:
        A list of rows with ``sample_index``, ``x``, ``y``, and ``time`` plus
        source metadata. ``NA`` payloads produce an empty list.

    Side effects:
        None.
    """

    rows: list[dict[str, Any]] = []
    for index, point in enumerate(parse_literal_sequence(text), start=1):
        if len(point) != 3:
            raise ValueError(f"timed point {index} has {len(point)} values; expected 3")
        row = base_output_row(metadata, zip_path, repo_root, internal_name)
        row.update(
            {
                "sample_index": index,
                "x": numeric_value(point[0]),
                "y": numeric_value(point[1]),
                "time": numeric_value(point[2]),
            }
        )
        rows.append(row)
    return rows


def parse_internal_table(
    table_name: str,
    text: str,
    metadata: dict[str, Any],
    zip_path: Path,
    repo_root: Path,
    internal_name: str,
) -> list[dict[str, Any]]:
    """Dispatch one internal TraceLab payload to the correct table parser.

    Args:
        table_name: One of ``points``, ``segments``, ``frames``, or
            ``tracings``.
        text: Internal TraceLab text payload.
        metadata: Parsed zip filename metadata.
        zip_path: Source zip path.
        repo_root: Repository root used for relative source paths.
        internal_name: Internal member name that supplied the text.

    Returns:
        Parsed output rows for the requested table.

    Side effects:
        None.
    """

    if table_name == "points":
        return parse_points_rows(text, metadata, zip_path, repo_root, internal_name)
    if table_name == "segments":
        return parse_segment_rows(text, metadata, zip_path, repo_root, internal_name)
    if table_name in {"frames", "tracings"}:
        return parse_timed_point_rows(text, metadata, zip_path, repo_root, internal_name)

    raise ValueError(f"unknown TraceLab table name: {table_name}")


class ParquetBatchWriter:
    """Write one reconstructed table to parquet in row batches."""

    def __init__(self, output_path: Path, columns: list[str], pa_module: Any, pq_module: Any) -> None:
        """Create a parquet batch writer for one output table.

        Args:
            output_path: Parquet path to write.
            columns: Stable output column order.
            pa_module: Imported ``pyarrow`` module.
            pq_module: Imported ``pyarrow.parquet`` module.

        Returns:
            None.

        Side effects:
            Stores writer configuration. The parquet file is opened lazily on
            the first non-empty batch.
        """

        self.output_path = output_path
        self.columns = columns
        self.pa = pa_module
        self.pq = pq_module
        self.schema = self._schema_for_columns(columns)
        self.writer: Any | None = None
        self.row_count = 0

    def _schema_for_columns(self, columns: list[str]) -> Any:
        """Build a pyarrow schema for one output table.

        Args:
            columns: Stable output column order.

        Returns:
            A ``pyarrow.Schema`` matching the table columns.

        Side effects:
            None.
        """

        int_columns = {"participant_id", "session", "block", "trial", "point_index", "segment_index", "sample_index"}
        float_columns = {"x", "y", "time", "start.x", "start.y", "end.x", "end.y", "ctrl.x", "ctrl.y"}
        bool_columns = {"learned_flag"}

        fields = []
        for column in columns:
            if column in int_columns:
                fields.append(self.pa.field(column, self.pa.int64()))
            elif column in float_columns:
                fields.append(self.pa.field(column, self.pa.float64()))
            elif column in bool_columns:
                fields.append(self.pa.field(column, self.pa.bool_()))
            else:
                fields.append(self.pa.field(column, self.pa.string()))
        return self.pa.schema(fields)

    def write_rows(self, rows: list[dict[str, Any]]) -> None:
        """Write a batch of rows to parquet.

        Args:
            rows: Output rows for this table.

        Returns:
            None.

        Side effects:
            Opens or appends to the parquet file at ``self.output_path``.
            Raises ``RuntimeError`` with path context if parquet writing fails.
        """

        if not rows:
            return

        try:
            frame = pd.DataFrame(rows)
            for column in self.columns:
                if column not in frame.columns:
                    frame[column] = None
            frame = frame[self.columns]
            table = self.pa.Table.from_pandas(frame, schema=self.schema, preserve_index=False)
            if self.writer is None:
                self.writer = self.pq.ParquetWriter(self.output_path, self.schema)
            self.writer.write_table(table)
            self.row_count += len(rows)
        except Exception as error:  # noqa: BLE001 - add output path context for parquet failures.
            raise RuntimeError(f"failed to write parquet output {self.output_path}: {error}") from error

    def close(self) -> None:
        """Close the parquet writer and create an empty parquet if needed.

        Args:
            None.

        Returns:
            None.

        Side effects:
            Finalizes the parquet file. If no rows were written, writes an
            empty parquet table with the expected schema.
        """

        if self.writer is not None:
            self.writer.close()
            return

        try:
            empty_table = self.pa.Table.from_pylist([], schema=self.schema)
            self.pq.write_table(empty_table, self.output_path)
        except Exception as error:  # noqa: BLE001 - add output path context for parquet failures.
            raise RuntimeError(f"failed to write empty parquet output {self.output_path}: {error}") from error


def issue_string(issue_codes: list[str]) -> str:
    """Join issue codes in stable order for a manifest cell.

    Args:
        issue_codes: Raw issue-code list.

    Returns:
        Semicolon-separated unique issue codes.

    Side effects:
        None.
    """

    return ";".join(sorted({code for code in issue_codes if code}))


def reconstruct_one_zip(zip_path: Path, repo_root: Path) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Reconstruct all available TraceLab tables from one zip file.

    Args:
        zip_path: Source TraceLab zip path.
        repo_root: Repository root used for relative manifest paths.

    Returns:
        ``(table_rows, manifest_row)`` where ``table_rows`` maps table names to
        parsed rows and ``manifest_row`` records filename, zip, internal-file,
        and parse statuses.

    Side effects:
        Opens the source zip read-only. Source data are not modified.
    """

    metadata = parse_zip_filename(zip_path)
    zip_stem = zip_path.stem
    table_rows: dict[str, list[dict[str, Any]]] = {table_name: [] for table_name in EXPECTED_INTERNAL_EXTENSIONS}
    issue_codes: list[str] = []

    if metadata["zip_filename_parse_status"] != "ok":
        issue_codes.append("zip_filename_parse_error")

    manifest_row: dict[str, Any] = {
        **metadata,
        "source_zip_path": relative_path(zip_path, repo_root),
        "source_zip_basename": zip_path.name,
        "zip_open_status": "not_opened",
        "internal_file_names_json": "[]",
    }

    for table_name in EXPECTED_INTERNAL_EXTENSIONS:
        manifest_row[f"{table_name}_internal_name"] = ""
        manifest_row[f"{table_name}_internal_present"] = False
        manifest_row[f"{table_name}_parse_status"] = "not_parsed"
        manifest_row[f"{table_name}_row_count"] = 0
        manifest_row[f"{table_name}_parse_error_message"] = ""

    try:
        with zipfile.ZipFile(zip_path, mode="r") as zip_handle:
            internal_names = zip_handle.namelist()
            manifest_row["zip_open_status"] = "ok"
            manifest_row["internal_file_names_json"] = json_dump(internal_names)

            for table_name, extension in EXPECTED_INTERNAL_EXTENSIONS.items():
                internal_name, internal_issue = find_internal_trace_file(internal_names, zip_stem, extension)
                if internal_issue:
                    issue_codes.append(internal_issue)
                if internal_name is None:
                    manifest_row[f"{table_name}_parse_status"] = "missing_internal_file"
                    continue

                manifest_row[f"{table_name}_internal_name"] = internal_name
                manifest_row[f"{table_name}_internal_present"] = True

                try:
                    text = read_internal_text(zip_handle, internal_name)
                    parsed_rows = parse_internal_table(table_name, text, metadata, zip_path, repo_root, internal_name)
                    table_rows[table_name] = parsed_rows
                    if table_name == "tracings" and text.strip() == "NA":
                        manifest_row[f"{table_name}_parse_status"] = "na_payload"
                    else:
                        manifest_row[f"{table_name}_parse_status"] = "ok"
                    manifest_row[f"{table_name}_row_count"] = len(parsed_rows)
                except Exception as error:  # noqa: BLE001 - record one table failure and keep processing.
                    issue_code = f"{table_name}_parse_error"
                    issue_codes.append(issue_code)
                    manifest_row[f"{table_name}_parse_status"] = "parse_error"
                    manifest_row[f"{table_name}_parse_error_message"] = f"{type(error).__name__}: {error}"

    except Exception as error:  # noqa: BLE001 - one bad zip should not stop the full run.
        issue_codes.append("zip_open_error")
        manifest_row["zip_open_status"] = "zip_open_error"
        manifest_row["zip_open_error_message"] = f"{type(error).__name__}: {error}"
    else:
        manifest_row["zip_open_error_message"] = ""

    manifest_row["issue_codes"] = issue_string(issue_codes)
    return table_rows, manifest_row


def flush_table_buffers(
    buffers: dict[str, list[dict[str, Any]]],
    writers: dict[str, ParquetBatchWriter],
    force: bool = False,
) -> None:
    """Write buffered table rows that have reached the batch threshold.

    Args:
        buffers: Pending rows by table name.
        writers: Parquet writers by table name.
        force: Whether to write all pending rows regardless of batch size.

    Returns:
        None.

    Side effects:
        Writes parquet batches and clears buffers that were written.
    """

    for table_name, rows in buffers.items():
        if rows and (force or len(rows) >= BATCH_ROW_LIMIT):
            writers[table_name].write_rows(rows)
            rows.clear()


def reconstruct_all_zips(zip_files: list[Path], repo_root: Path, output_dir: Path) -> pd.DataFrame:
    """Reconstruct TraceLab tables for all discovered zip files.

    Args:
        zip_files: Sorted source zip paths.
        repo_root: Repository root used for relative source paths.
        output_dir: Directory where parquet outputs are written.

    Returns:
        Manifest ``DataFrame`` with one row per source zip.

    Side effects:
        Opens zip files read-only. Writes four parquet files below
        ``output_dir``.
    """

    pa_module, pq_module = require_pyarrow_modules()
    output_dir.mkdir(parents=True, exist_ok=True)

    writers = {
        table_name: ParquetBatchWriter(output_dir / output_filename, TABLE_COLUMNS[table_name], pa_module, pq_module)
        for table_name, output_filename in TABLE_OUTPUT_FILENAMES.items()
    }
    buffers: dict[str, list[dict[str, Any]]] = {table_name: [] for table_name in TABLE_OUTPUT_FILENAMES}
    manifest_rows: list[dict[str, Any]] = []

    try:
        for index, zip_path in enumerate(zip_files, start=1):
            if index == 1 or index % 500 == 0 or index == len(zip_files):
                print(f"Reconstructing TraceLab zip {index:,} of {len(zip_files):,}: {zip_path.name}")

            table_rows, manifest_row = reconstruct_one_zip(zip_path, repo_root)
            manifest_rows.append(manifest_row)
            for table_name, rows in table_rows.items():
                buffers[table_name].extend(rows)
            flush_table_buffers(buffers, writers)

        flush_table_buffers(buffers, writers, force=True)
    finally:
        for writer in writers.values():
            writer.close()

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.attrs["table_row_counts"] = {table_name: writer.row_count for table_name, writer in writers.items()}
    return manifest_df


def issue_counts_from_manifest(manifest_df: pd.DataFrame) -> Counter[str]:
    """Count issue codes from the zip parse manifest.

    Args:
        manifest_df: One-row-per-zip manifest.

    Returns:
        A ``Counter`` keyed by issue code.

    Side effects:
        None.
    """

    counts: Counter[str] = Counter()
    if manifest_df.empty or "issue_codes" not in manifest_df.columns:
        return counts

    for value in manifest_df["issue_codes"].fillna(""):
        for code in str(value).split(";"):
            if code:
                counts[code] += 1
    return counts


def write_manifest(manifest_df: pd.DataFrame, output_path: Path) -> None:
    """Write the zip parse manifest as CSV.

    Args:
        manifest_df: One-row-per-zip manifest.
        output_path: CSV path to write.

    Returns:
        None.

    Side effects:
        Writes ``output_path``.
    """

    manifest_df.to_csv(output_path, index=False)


def write_summary(manifest_df: pd.DataFrame, output_path: Path, table_row_counts: dict[str, int], started_at: datetime) -> None:
    """Write a Markdown summary of TraceLab reconstruction outputs.

    Args:
        manifest_df: One-row-per-zip manifest.
        output_path: Markdown path to write.
        table_row_counts: Number of rows written to each parquet table.
        started_at: Timestamp captured when the run began.

    Returns:
        None.

    Side effects:
        Writes ``output_path``.
    """

    finished_at = datetime.now()
    issue_counts = issue_counts_from_manifest(manifest_df)
    zip_count = len(manifest_df)
    learned_count = int(manifest_df["learned_flag"].fillna(False).sum()) if zip_count else 0
    zip_open_errors = int((manifest_df["zip_open_status"] != "ok").sum()) if zip_count else 0
    filename_parse_errors = (
        int((manifest_df["zip_filename_parse_status"] != "ok").sum()) if zip_count else 0
    )

    parse_status_lines: list[str] = []
    for table_name in EXPECTED_INTERNAL_EXTENSIONS:
        status_counts = Counter(manifest_df[f"{table_name}_parse_status"].fillna("")) if zip_count else Counter()
        status_text = ", ".join(f"{status}: {count}" for status, count in sorted(status_counts.items()))
        parse_status_lines.append(f"- {table_name}: {status_text if status_text else 'no rows'}")

    issue_lines = [f"- {code}: {count}" for code, count in sorted(issue_counts.items())]
    if not issue_lines:
        issue_lines = ["- none"]

    row_count_lines = [f"- {table_name}: {count:,}" for table_name, count in sorted(table_row_counts.items())]

    summary = f"""# TraceLab Reconstruction Summary

Generated: {finished_at.isoformat(timespec="seconds")}

Started: {started_at.isoformat(timespec="seconds")}

## Scope

This local summary describes read-only reconstruction of TraceLab zip contents
from `_Data/figure/` into parquet tables under `_Data/behavior/tracing_tables/`.
It reports structural import facts only. No tracing filters, behavioural
outcome calculations, raw EEG annotation checks, or EEG preprocessing were run.

## Zip Files

- Zip files discovered: {zip_count:,}
- Zip files with `learned` in the basename: {learned_count:,}
- Zip open errors: {zip_open_errors:,}
- Filename parse errors: {filename_parse_errors:,}

## Output Row Counts

{chr(10).join(row_count_lines)}

## Internal Parse Status Counts

{chr(10).join(parse_status_lines)}

## Issue Codes

{chr(10).join(issue_lines)}

## Next Step

The next linkage-recovery step is to port the old tracing filters from
`_Scripts/01_preprocessing.R` using these reconstructed tables as input. Raw
EEG annotations should be compared only after tracing filters and event-offset
tables have been reconstructed.
"""

    output_path.write_text(summary, encoding="utf-8")


def main() -> int:
    """Run TraceLab zip-content reconstruction.

    Args:
        None.

    Returns:
        ``0`` when reconstruction completes.

    Side effects:
        Reads TraceLab zip files below ``_Data/figure/`` and writes local
        parquet/CSV/Markdown outputs below ``_Data/behavior/tracing_tables/``.
    """

    started_at = datetime.now()
    repo_root = repo_root_from_script()
    figure_dir = repo_root / "_Data" / "figure"
    output_dir = repo_root / "_Data" / "behavior" / "tracing_tables"

    zip_files = discover_figure_zip_files(figure_dir)
    manifest_df = reconstruct_all_zips(zip_files, repo_root, output_dir)

    table_row_counts = manifest_df.attrs.get("table_row_counts", {})
    manifest_path = output_dir / "zip_parse_manifest.csv"
    summary_path = output_dir / "tracing_reconstruction_summary.md"
    write_manifest(manifest_df, manifest_path)
    write_summary(manifest_df, summary_path, table_row_counts, started_at)

    print(f"Wrote {manifest_path.relative_to(repo_root)}")
    print(f"Wrote {summary_path.relative_to(repo_root)}")
    for table_name, output_filename in TABLE_OUTPUT_FILENAMES.items():
        print(f"Wrote {(output_dir / output_filename).relative_to(repo_root)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
