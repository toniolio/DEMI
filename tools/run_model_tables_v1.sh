#!/usr/bin/env bash
# Run accepted DEMI model_tables_v1 construction with durable launcher logging.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${REPO_ROOT}/_Data/eeg/model_tables_v1"
LAUNCHER_LOG_ROOT="${REPO_ROOT}/_Data/eeg/model_tables_v1_launcher_logs"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LAUNCHER_LOG_ROOT}/terminal_run_${TIMESTAMP}.log"
MIN_FREE_KIB=$((1 * 1024 * 1024))

cd "${REPO_ROOT}"

if ! git check-ignore --quiet "_Data/eeg/model_tables_v1"; then
  echo "FAIL: _Data/eeg/model_tables_v1 is not ignored by git." >&2
  exit 1
fi

mkdir -p "${LAUNCHER_LOG_ROOT}"
FREE_KIB="$(df -Pk "${LAUNCHER_LOG_ROOT}" | awk 'NR == 2 {print $4}')"
if (( FREE_KIB < MIN_FREE_KIB )); then
  echo "FAIL: model_tables_v1 requires at least 1 GiB free; observed ${FREE_KIB} KiB." >&2
  exit 1
fi

PYTHON="${REPO_ROOT}/.venv/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
  echo "FAIL: expected EEG environment is absent: ${PYTHON}" >&2
  exit 1
fi

export PYTHONUNBUFFERED=1

echo "DEMI model_tables_v1 launcher"
echo "Repository: ${REPO_ROOT}"
echo "Output: ${OUTPUT_ROOT}"
echo "Launcher log: ${LOG_PATH}"
echo "Free space: ${FREE_KIB} KiB"

COMMAND=("${PYTHON}" "analysis/eeg_mne/18_construct_model_ready_tables.py" "$@")
"${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
