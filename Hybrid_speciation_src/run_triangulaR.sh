#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${BASE_DIR:-.}"
TRIANGULAR_SCRIPT="${TRIANGULAR_SCRIPT:-run_triangulaR.R}"

cd "${BASE_DIR}"

if [[ ! -f "${TRIANGULAR_SCRIPT}" ]]; then
    echo "Missing triangulaR script: ${TRIANGULAR_SCRIPT}" >&2
    exit 1
fi

Rscript "${TRIANGULAR_SCRIPT}"
