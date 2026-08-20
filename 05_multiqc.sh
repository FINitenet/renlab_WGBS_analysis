#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib.sh"

REPORT_DIR="${OUTDIR}/06_multiqc"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${REPORT_DIR}" "${LOG_DIR}"
require_command multiqc

log "Aggregate FastQC, Trim Galore, Bismark alignment/dedup/extraction reports"
inputs=("${OUTDIR}")
[[ -z "${RAW_QC_DIR:-}" ]] || inputs+=("${RAW_QC_DIR}")
multiqc --force --outdir "${REPORT_DIR}" "${inputs[@]}" \
    >"${LOG_DIR}/06_multiqc.log" 2>&1
log "Combined report: ${REPORT_DIR}/multiqc_report.html"
