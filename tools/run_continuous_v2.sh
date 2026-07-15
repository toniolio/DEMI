#!/bin/bash

# Launch the authoritative DEMI continuous-v2 preprocessing surface.
#
# The script is intentionally self-locating and contains no continued shell
# lines, so it is safe to invoke from any directory and cannot be broken by
# whitespace introduced while pasting a multiline command. The Python driver
# owns per-recording atomic publication, interruption recovery, cache checks,
# and full raw-source immutability snapshots.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${REPO_ROOT}/_Data/eeg/mne_preprocessing/continuous_v2"
LOG_PATH="${OUTPUT_ROOT}/terminal_run_$(date -u +%Y%m%dT%H%M%SZ).log"
PYTHON="${REPO_ROOT}/.venv/bin/python"
DRIVER="${REPO_ROOT}/analysis/eeg_mne/13_run_continuous_preprocessing_validation.py"

cd "${REPO_ROOT}"
mkdir -p "${OUTPUT_ROOT}"

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: EEG Python environment is missing or not executable: ${PYTHON}" >&2
    exit 1
fi

if ! command -v caffeinate >/dev/null 2>&1; then
    echo "ERROR: caffeinate is required for the unattended macOS run." >&2
    exit 1
fi

echo "Repository root: ${REPO_ROOT}"
echo "Target output root: ${OUTPUT_ROOT}"
echo "Terminal log: ${LOG_PATH}"
printf 'Command:'
printf ' %q' caffeinate -dimsu "${PYTHON}" -u "${DRIVER}" --all-recordings
printf '\n'

caffeinate -dimsu "${PYTHON}" -u "${DRIVER}" --all-recordings 2>&1 | tee "${LOG_PATH}"
