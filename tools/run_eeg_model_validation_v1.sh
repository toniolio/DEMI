#!/usr/bin/env bash
# Run DEMI model_validation_v1 with durable launcher logging.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_ROOT="${REPO_ROOT}/_Data/eeg/model_validation_v1"
LAUNCHER_LOG_ROOT="${REPO_ROOT}/_Data/eeg/model_validation_v1_launcher_logs"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LAUNCHER_LOG_ROOT}/terminal_run_${TIMESTAMP}.log"
MIN_FREE_KIB=$((2 * 1024 * 1024))

cd "${REPO_ROOT}"

if ! git check-ignore --quiet "_Data/eeg/model_validation_v1"; then
  echo "FAIL: _Data/eeg/model_validation_v1 is not ignored by git." >&2
  exit 1
fi

mkdir -p "${LAUNCHER_LOG_ROOT}"
FREE_KIB="$(df -Pk "${LAUNCHER_LOG_ROOT}" | awk 'NR == 2 {print $4}')"
if (( FREE_KIB < MIN_FREE_KIB )); then
  echo "FAIL: model_validation_v1 requires at least 2 GiB free; observed ${FREE_KIB} KiB." >&2
  exit 1
fi

if [[ ! -x "${REPO_ROOT}/.venv/bin/python" ]]; then
  echo "FAIL: expected locked Python environment is absent." >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "FAIL: Rscript is unavailable." >&2
  exit 1
fi

echo "DEMI model_validation_v1 launcher"
echo "Repository: ${REPO_ROOT}"
echo "Output: ${OUTPUT_ROOT}"
echo "Launcher log: ${LOG_PATH}"
echo "Free space: ${FREE_KIB} KiB"

COMMAND=(
  Rscript --no-init-file
  analysis/eeg_models/19_validate_synthetic_models_and_priors.R
  "$@"
)
"${COMMAND[@]}" 2>&1 | tee "${LOG_PATH}"

