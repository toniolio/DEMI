#!/usr/bin/env bash
# Run accepted DEMI features_v1 construction with durable progress logging.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${REPO_ROOT}/_Data/eeg/features_v1"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${OUTPUT_ROOT}/terminal_run_${TIMESTAMP}.log"
MIN_FREE_KIB=$((5 * 1024 * 1024))

cd "${REPO_ROOT}"
mkdir -p "${OUTPUT_ROOT}"

if ! git check-ignore --quiet "_Data/eeg/features_v1"; then
  echo "FAIL: _Data/eeg/features_v1 is not ignored by git." >&2
  exit 1
fi

FREE_KIB="$(df -Pk "${OUTPUT_ROOT}" | awk 'NR == 2 {print $4}')"
if (( FREE_KIB < MIN_FREE_KIB )); then
  echo "FAIL: features_v1 requires at least 5 GiB free; observed ${FREE_KIB} KiB." >&2
  exit 1
fi

PYTHON="${REPO_ROOT}/.venv/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
  echo "FAIL: expected EEG environment is absent: ${PYTHON}" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "FAIL: Rscript is required to read the frozen bdat2.rds authority." >&2
  exit 1
fi

export PYTHONUNBUFFERED=1
export R_PROFILE_USER=/dev/null
export MPLCONFIGDIR="${TMPDIR:-/tmp}/demi-features-mpl-cache"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/demi-features-xdg-cache"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

echo "DEMI features_v1 launcher"
echo "Repository: ${REPO_ROOT}"
echo "Output: ${OUTPUT_ROOT}"
echo "Log: ${LOG_PATH}"
echo "Free space: ${FREE_KIB} KiB"

COMMAND=("${PYTHON}" "analysis/eeg_mne/17_construct_trial_level_features.py" "$@")
if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -dimsu "${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
else
  "${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
fi
