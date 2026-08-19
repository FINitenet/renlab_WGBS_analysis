#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib.sh"

MODE="${1:-extract}"
[[ "${MODE}" == "mbias" || "${MODE}" == "extract" ]] || die "Usage: $0 [mbias|extract]"
DEDUP_DIR="${OUTDIR}/04_dedup"
METH_DIR="${OUTDIR}/05_methylation"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${METH_DIR}" "${LOG_DIR}"
require_command bismark_methylation_extractor
[[ -n "${BISMARK_GENOME_DIR}" ]] || die "BISMARK_GENOME_DIR is not configured"
has_bismark_bt2_index || die "Incomplete Bismark index: ${BISMARK_GENOME_DIR}"

if [[ "${MODE}" == "extract" && "${MBIAS_REVIEWED}" != "true" ]]; then
    die "M-bias gate not approved. Run '$0 mbias', inspect reports, set IGNORE_* in config.sh, then set MBIAS_REVIEWED=true."
fi

extract_sample() {
    local sample="$1" r1="$2" r2="$3"
    # Paired-end extraction requires mates to remain adjacent. A coordinate-
    # sorted BAM is retained separately for browsing and interval operations.
    local input_bam="${DEDUP_DIR}/${sample}/${sample}.deduplicated.namesorted.bam"
    local sample_dir="${METH_DIR}/${sample}" report marker
    require_file "${input_bam}"
    mkdir -p "${sample_dir}"

    local cmd=(bismark_methylation_extractor --comprehensive
               --parallel "${EXTRACT_PARALLEL}" --buffer_size "${EXTRACT_BUFFER}"
               --output "${sample_dir}" --genome_folder "${BISMARK_GENOME_DIR}")
    [[ -z "${r2}" ]] || cmd+=(--paired --no_overlap)

    if [[ "${MODE}" == "mbias" ]]; then
        report="${sample_dir}/${sample}.deduplicated.namesorted.M-bias.txt"
        marker="${sample_dir}/.mbias.complete"
        if [[ -s "${marker}" && -s "${report}" ]]; then
            log "M-bias skip (complete): ${sample}"
            return
        fi
        log "M-bias diagnostic: ${sample}"
        "${cmd[@]}" --mbias_only "${input_bam}" >"${LOG_DIR}/${sample}.mbias.log" 2>&1
        require_file "${report}"
        printf 'complete\n' > "${marker}"
        return
    fi

    report="${sample_dir}/${sample}.deduplicated.namesorted_splitting_report.txt"
    marker="${sample_dir}/.extract.complete"
    if [[ -s "${marker}" && -s "${report}" ]]; then
        log "Methylation extraction skip (complete): ${sample}"
        return
    fi
    log "Plant methylation extraction (CpG/CHG/CHH): ${sample}"
    local clips=(--ignore "${IGNORE_R1_5}" --ignore_3prime "${IGNORE_R1_3}")
    if [[ -n "${r2}" ]]; then
        clips+=(--ignore_r2 "${IGNORE_R2_5}" --ignore_3prime_r2 "${IGNORE_R2_3}")
    fi
    "${cmd[@]}" --bedGraph --cytosine_report --CX --gzip --report \
        "${clips[@]}" "${input_bam}" >"${LOG_DIR}/${sample}.methylation_extractor.log" 2>&1
    require_file "${report}"
    printf 'complete\n' > "${marker}"
}

for_each_sample extract_sample
if [[ "${MODE}" == "mbias" ]]; then
    log "M-bias diagnostics complete. Review ${METH_DIR}/*/*M-bias* before extraction."
else
    log "Methylation extraction complete: ${METH_DIR}"
fi
