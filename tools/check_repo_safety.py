"""Check that private notes and local DEMI data are not tracked by git.

This script is a lightweight repository safety check for the DEMI EEG
reanalysis repository. It is meant to catch two common mistakes before a commit:

- private workflow files, participant-facing notes, local raw data, or generated
  manifests becoming tracked by git;
- root ignore rules being weakened so private/local paths are no longer ignored.

The check is intentionally conservative and boring. It uses git itself as the
source of truth rather than trying to reimplement ignore matching in Python:

- ``git ls-files`` lists files already tracked by the repository;
- ``git check-ignore`` confirms that expected local/private paths are ignored;
- ``git rev-parse --show-toplevel`` makes the script work from the repository
  root or from any subdirectory.

The only tracked file allowed under ``_Data/`` is ``_Data/eeg/BESA-81.csv``.
That file is a public coordinate/layout resource used by the analysis code, not
raw participant data or a generated local manifest.

Run from anywhere inside the repository with:

    python tools/check_repo_safety.py

The script prints a concise PASS/FAIL summary and exits nonzero on failure.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


# These paths should be ignored in a clean checkout. They include private notes,
# local raw data, generated manifests, local virtual environments, and OS files.
EXPECTED_IGNORED_PATHS = (
    "AGENTS.md",
    "_Private/",
    ".private/",
    ".venv/",
    ".DS_Store",
    "_Data/eeg/raw/",
    "_Data/eeg/manifest/",
    "_Data/behavior/manifest/",
    "_Data/task/",
    "_Data/figure/",
)

# This is the one known public/tracked file that legitimately lives below
# _Data/. The tracked-file audit treats every other _Data path as suspicious.
ALLOWED_TRACKED_DATA_FILES = {"_Data/eeg/BESA-81.csv"}

# Exact tracked paths that should never be committed. Directory-like cases are
# handled by DISALLOWED_TRACKED_PREFIXES below.
DISALLOWED_TRACKED_EXACT = {"AGENTS.md"}

# Directory prefixes that should never contain tracked files. Prefixes include
# a trailing slash so similarly named files do not match accidentally.
DISALLOWED_TRACKED_PREFIXES = (
    "_Private/",
    ".private/",
    ".venv/",
    "_Data/eeg/raw/",
    "_Data/eeg/manifest/",
    "_Data/behavior/manifest/",
    "_Data/task/",
    "_Data/figure/",
)


def run_git(args: list[str], repo_root: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    """Run a git command and return the completed process.

    Args:
        args: Git arguments excluding the leading ``git`` executable.
        repo_root: Directory from which to run the command. ``None`` means use
            the current process working directory.
        check: Whether to raise ``CalledProcessError`` for nonzero exit codes.

    Returns:
        A ``CompletedProcess`` with text-mode stdout and stderr.

    Side effects:
        Spawns a subprocess. The repository is only inspected, not modified.
    """

    return subprocess.run(
        ["git", *args],
        cwd=repo_root,
        check=check,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def find_repo_root() -> Path:
    """Find the git repository root for the current working directory.

    Args:
        None.

    Returns:
        The absolute path reported by ``git rev-parse --show-toplevel``.

    Side effects:
        Runs ``git rev-parse``. Exits through the caller if the current
        directory is not inside a git repository.
    """

    result = run_git(["rev-parse", "--show-toplevel"])
    return Path(result.stdout.strip()).resolve()


def list_tracked_files(repo_root: Path) -> list[str]:
    """Return all files tracked by git, using repository-relative paths.

    Args:
        repo_root: Absolute path to the repository root.

    Returns:
        A sorted list of tracked file paths with POSIX-style separators.

    Side effects:
        Runs ``git ls-files``. It does not inspect untracked files directly.
    """

    result = run_git(["ls-files"], repo_root=repo_root)
    return sorted(path for path in result.stdout.splitlines() if path)


def tracked_safety_violations(tracked_files: list[str]) -> list[str]:
    """Identify tracked files that violate private/local-data safety rules.

    Args:
        tracked_files: Repository-relative paths returned by ``git ls-files``.

    Returns:
        A list of human-readable violation messages.

    Side effects:
        None.
    """

    violations: list[str] = []

    for path in tracked_files:
        if path in ALLOWED_TRACKED_DATA_FILES:
            continue

        if path in DISALLOWED_TRACKED_EXACT:
            violations.append(f"tracked private/local file: {path}")
            continue

        if Path(path).name == ".DS_Store":
            violations.append(f"tracked OS metadata file: {path}")
            continue

        if path.startswith("_Data/"):
            violations.append(f"tracked local data file outside allowlist: {path}")
            continue

        for prefix in DISALLOWED_TRACKED_PREFIXES:
            if path.startswith(prefix):
                violations.append(f"tracked file under ignored prefix {prefix}: {path}")
                break

    return violations


def check_expected_ignored_paths(repo_root: Path) -> list[str]:
    """Verify that expected private/local paths are ignored by git.

    Args:
        repo_root: Absolute path to the repository root.

    Returns:
        A list of human-readable violation messages for paths that are not
        ignored according to ``git check-ignore``.

    Side effects:
        Runs ``git check-ignore`` once per expected path. This avoids relying on
        path existence, because some safety paths may not exist in every clone.
    """

    violations: list[str] = []

    for path in EXPECTED_IGNORED_PATHS:
        result = run_git(["check-ignore", "--quiet", path], repo_root=repo_root, check=False)
        if result.returncode != 0:
            violations.append(f"expected ignored path is not ignored: {path}")

    return violations


def check_allowed_data_file(tracked_files: list[str], repo_root: Path) -> list[str]:
    """Ensure the public BESA coordinate file remains tracked and unignored.

    Args:
        tracked_files: Repository-relative paths returned by ``git ls-files``.
        repo_root: Absolute path to the repository root.

    Returns:
        A list of violation messages if the BESA file is missing from tracking
        or is matched by ignore rules.

    Side effects:
        Runs ``git check-ignore`` for ``_Data/eeg/BESA-81.csv``.
    """

    besa_path = "_Data/eeg/BESA-81.csv"
    violations: list[str] = []

    if besa_path not in tracked_files:
        violations.append(f"allowed coordinate file is not tracked: {besa_path}")

    # Tracked files can remain tracked even if a later ignore rule matches them.
    # This explicit check catches accidental ignore-rule changes before someone
    # removes or tries to re-add the coordinate file.
    result = run_git(["check-ignore", "--quiet", besa_path], repo_root=repo_root, check=False)
    if result.returncode == 0:
        violations.append(f"allowed coordinate file is unexpectedly ignored: {besa_path}")

    return violations


def print_summary(violations: list[str]) -> None:
    """Print a concise PASS/FAIL summary for humans and commit hooks.

    Args:
        violations: Human-readable safety violations collected by the checks.

    Returns:
        None.

    Side effects:
        Writes the summary to stdout.
    """

    if not violations:
        print("PASS repo safety: private/local data paths are untracked and ignored as expected.")
        return

    print(f"FAIL repo safety: {len(violations)} issue(s) found.")
    for violation in violations:
        print(f"- {violation}")


def main() -> int:
    """Run all repository safety checks.

    Args:
        None.

    Returns:
        ``0`` when all checks pass, otherwise ``1``.

    Side effects:
        Runs read-only git commands and prints a PASS/FAIL summary.
    """

    try:
        repo_root = find_repo_root()
        tracked_files = list_tracked_files(repo_root)
    except subprocess.CalledProcessError as error:
        stderr = error.stderr.strip()
        print(f"FAIL repo safety: unable to inspect git repository. {stderr}")
        return 1

    violations: list[str] = []
    violations.extend(tracked_safety_violations(tracked_files))
    violations.extend(check_expected_ignored_paths(repo_root))
    violations.extend(check_allowed_data_file(tracked_files, repo_root))

    print_summary(violations)
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
