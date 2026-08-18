#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib.sh"

QC_DIR="${OUTDIR}/01_qc"
FASTQC_DIR="${QC_DIR}/fastqc"
MULTIQC_DIR="${QC_DIR}/multiqc"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${FASTQC_DIR}" "${MULTIQC_DIR}" "${LOG_DIR}"
require_command fastqc multiqc

qc_sample() {
    local sample="$1" r1="$2" r2="$3" r1_zip r2_zip="" marker
    marker="${FASTQC_DIR}/.${sample}.complete"
    r1_zip="${FASTQC_DIR}/$(basename "${r1}" | sed -E 's/\.(fastq|fq)(\.gz)?$/_fastqc.zip/i')"
    [[ -z "${r2}" ]] || r2_zip="${FASTQC_DIR}/$(basename "${r2}" | sed -E 's/\.(fastq|fq)(\.gz)?$/_fastqc.zip/i')"
    if [[ -s "${marker}" && -s "${r1_zip}" && ( -z "${r2}" || -s "${r2_zip}" ) ]]; then
        log "FastQC skip (complete): ${sample}"
        return
    fi
    log "FastQC: ${sample}"
    if [[ -n "${r2}" ]]; then
        fastqc --threads "${FASTQC_THREADS}" --outdir "${FASTQC_DIR}" "${r1}" "${r2}" \
            >"${LOG_DIR}/${sample}.fastqc.log" 2>&1
    else
        fastqc --threads "${FASTQC_THREADS}" --outdir "${FASTQC_DIR}" "${r1}" \
            >"${LOG_DIR}/${sample}.fastqc.log" 2>&1
    fi
    [[ -s "${r1_zip}" && ( -z "${r2}" || -s "${r2_zip}" ) ]] || die "FastQC output incomplete: ${sample}"
    printf 'complete\n' > "${marker}"
}

log "Raw-read QC from ${SAMPLES_TSV}"
for_each_sample qc_sample
log "MultiQC summary"
multiqc --force --outdir "${MULTIQC_DIR}" "${FASTQC_DIR}" \
    >"${LOG_DIR}/01_multiqc.log" 2>&1
log "QC complete: ${MULTIQC_DIR}/multiqc_report.html"
