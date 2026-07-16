#!/usr/bin/env bash
# Run the accepted DEMI epochs_v1 construction with durable progress logging.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${REPO_ROOT}/_Data/eeg/epochs_v1"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${OUTPUT_ROOT}/terminal_run_${TIMESTAMP}.log"
MIN_FREE_KIB=$((30 * 1024 * 1024))

cd "${REPO_ROOT}"
mkdir -p "${OUTPUT_ROOT}"

if ! git check-ignore --quiet "_Data/eeg/epochs_v1"; then
  echo "FAIL: _Data/eeg/epochs_v1 is not ignored by git." >&2
  exit 1
fi

FREE_KIB="$(df -Pk "${OUTPUT_ROOT}" | awk 'NR == 2 {print $4}')"
if (( FREE_KIB < MIN_FREE_KIB )); then
  echo "FAIL: epochs_v1 requires at least 30 GiB free; observed ${FREE_KIB} KiB." >&2
  exit 1
fi

export PYTHONUNBUFFERED=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/demi-mpl-cache"
mkdir -p "${MPLCONFIGDIR}"

PYTHON="${REPO_ROOT}/.venv/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
  echo "FAIL: expected EEG environment is absent: ${PYTHON}" >&2
  exit 1
fi

echo "DEMI epochs_v1 launcher"
echo "Repository: ${REPO_ROOT}"
echo "Output: ${OUTPUT_ROOT}"
echo "Log: ${LOG_PATH}"
echo "Free space: ${FREE_KIB} KiB"

COMMAND=("${PYTHON}" "analysis/eeg_mne/15_construct_accepted_epochs.py" "$@")
if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -dimsu "${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
else
  "${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
fi
