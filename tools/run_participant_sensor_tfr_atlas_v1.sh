#!/usr/bin/env bash
# Run the read-only participant TFR visual-QC atlas with durable terminal logging.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${REPO_ROOT}/_Data/eeg/participant_sensor_tfr_atlas_v1"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${TMPDIR:-/tmp}/demi-participant-atlas-terminal-${TIMESTAMP}.log"

cd "${REPO_ROOT}"
if ! git check-ignore --quiet "_Data/eeg/participant_sensor_tfr_atlas_v1"; then
  echo "FAIL: atlas output root is not ignored by git." >&2
  exit 1
fi

PYTHON="${REPO_ROOT}/.venv/bin/python"
if [[ ! -x "${PYTHON}" ]]; then
  echo "FAIL: expected EEG environment is absent: ${PYTHON}" >&2
  exit 1
fi

export PYTHONUNBUFFERED=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/demi-participant-atlas-mpl-cache"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/demi-participant-atlas-xdg-cache"
mkdir -p "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}"

echo "DEMI participant sensor TFR atlas v1"
echo "Repository: ${REPO_ROOT}"
echo "Output: ${OUTPUT_ROOT}"
echo "Log: ${LOG_PATH}"

COMMAND=("${PYTHON}" "analysis/eeg_mne/20_construct_participant_sensor_tfr_atlas.py" "$@")
if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -dimsu "${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
else
  "${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"
fi
