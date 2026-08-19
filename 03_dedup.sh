#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib.sh"

MAP_DIR="${OUTDIR}/03_align"
DEDUP_DIR="${OUTDIR}/04_dedup"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${DEDUP_DIR}" "${LOG_DIR}"
require_command deduplicate_bismark samtools
[[ "${ENABLE_DEDUP}" == "true" ]] || die \
    "ENABLE_DEDUP=false. This is correct for RRBS; do not run coordinate deduplication."

dedup_sample() {
    local sample="$1" r1="$2" r2="$3"
    local input_bam="${MAP_DIR}/${sample}/${sample}.bam" sample_dir="${DEDUP_DIR}/${sample}"
    local dedup_bam="${sample_dir}/${sample}.deduplicated.bam"
    local sorted_bam="${sample_dir}/${sample}.deduplicated.sorted.bam"
    local namesorted_bam="${sample_dir}/${sample}.deduplicated.namesorted.bam"
    local marker="${sample_dir}/.complete"
    require_file "${input_bam}"
    mkdir -p "${sample_dir}"

    if [[ ! -s "${dedup_bam}" ]] || ! samtools quickcheck "${dedup_bam}"; then
        log "Deduplicate: ${sample}"
        local cmd=(deduplicate_bismark --bam --output_dir "${sample_dir}" --outfile "${sample}")
        [[ -z "${r2}" ]] || cmd+=(--paired)
        cmd+=("${input_bam}")
        "${cmd[@]}" >"${LOG_DIR}/${sample}.deduplicate.log" 2>&1
    else
        log "Dedup skip (complete): ${sample}"
    fi
    require_file "${dedup_bam}"

    if [[ ! -s "${sorted_bam}" || ! -s "${sorted_bam}.bai" ]]; then
        log "Coordinate sort/index: ${sample}"
        samtools sort --threads "${SAMTOOLS_THREADS}" --output-fmt BAM \
            -o "${sorted_bam}.tmp" "${dedup_bam}"
        mv "${sorted_bam}.tmp" "${sorted_bam}"
        samtools index -@ "${SAMTOOLS_THREADS}" "${sorted_bam}"
    else
        samtools quickcheck -v "${sorted_bam}"
        log "Sort/index skip (complete): ${sample}"
    fi

    # bismark_methylation_extractor requires paired mates to be adjacent.
    # Keep the coordinate-sorted/indexed BAM for genome browsers and downstream
    # interval tools, and create a separate query-name-sorted BAM for extraction.
    if [[ ! -s "${namesorted_bam}" ]] || ! samtools quickcheck "${namesorted_bam}"; then
        log "Query-name sort for methylation extraction: ${sample}"
        samtools sort -n --threads "${SAMTOOLS_THREADS}" --output-fmt BAM \
            -o "${namesorted_bam}.tmp" "${dedup_bam}"
        mv "${namesorted_bam}.tmp" "${namesorted_bam}"
    else
        log "Query-name sort skip (complete): ${sample}"
    fi
    samtools quickcheck -v "${namesorted_bam}"
    printf 'complete\n' > "${marker}"
}

for_each_sample dedup_sample
log "Deduplication complete: ${DEDUP_DIR}"
