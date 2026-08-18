#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib.sh"

TRIM_DIR="${OUTDIR}/02_trim"
MAP_DIR="${OUTDIR}/03_align"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${TRIM_DIR}" "${MAP_DIR}" "${LOG_DIR}"
require_command trim_galore bismark bismark_genome_preparation bowtie2 samtools
[[ -n "${GENOME_FA}" ]] || die "GENOME_FA is not configured"
[[ -n "${BISMARK_GENOME_DIR}" ]] || die "BISMARK_GENOME_DIR is not configured"
require_file "${GENOME_FA}"

[[ "${ENABLE_DECONTAM}" == "false" ]] || die \
    "ENABLE_DECONTAM=true is unsupported: ordinary Bowtie2 is not bisulfite-aware. Build a Bismark contaminant index instead."

case "${LIBRARY_TYPE}" in
    directional) BM_DIRECTION=() ;;
    non_directional) BM_DIRECTION=(--non_directional) ;;
    pbat) BM_DIRECTION=(--pbat) ;;
    *) die "LIBRARY_TYPE must be directional, non_directional, or pbat" ;;
esac

if ! has_bismark_bt2_index; then
    log "Complete Bismark Bowtie2 index not found; preparing ${BISMARK_GENOME_DIR}"
    mkdir -p "${BISMARK_GENOME_DIR}"
    if [[ ! -e "${BISMARK_GENOME_DIR}/$(basename "${GENOME_FA}")" ]]; then
        ln -s "${GENOME_FA}" "${BISMARK_GENOME_DIR}/$(basename "${GENOME_FA}")"
    fi
    bismark_genome_preparation --verbose --bowtie2 "${BISMARK_GENOME_DIR}" \
        >"${LOG_DIR}/bismark_genome_preparation.log" 2>&1
    has_bismark_bt2_index || die "Bismark genome preparation did not create a complete Bowtie2 index"
fi

process_sample() {
    local sample="$1" r1="$2" r2="$3"
    local trim_sample_dir="${TRIM_DIR}/${sample}" map_sample_dir="${MAP_DIR}/${sample}"
    local r1_trim r2_trim bam_out="${map_sample_dir}/${sample}.bam"
    local trim_marker="${trim_sample_dir}/.complete" map_marker="${map_sample_dir}/.complete"
    mkdir -p "${trim_sample_dir}" "${map_sample_dir}"

    if [[ -n "${r2}" ]]; then
        r1_trim="${trim_sample_dir}/$(basename "${r1}" | sed -E 's/\.(fastq|fq)(\.gz)?$/_val_1.fq.gz/i')"
        r2_trim="${trim_sample_dir}/$(basename "${r2}" | sed -E 's/\.(fastq|fq)(\.gz)?$/_val_2.fq.gz/i')"
        if [[ ! -s "${trim_marker}" || ! -s "${r1_trim}" || ! -s "${r2_trim}" ]]; then
            log "Trim Galore PE: ${sample}"
            trim_galore --paired --quality "${TRIM_Q}" --length "${TRIM_MINLEN}" \
                --cores "${TRIM_CORES}" --gzip --output_dir "${trim_sample_dir}" "${r1}" "${r2}" \
                >"${LOG_DIR}/${sample}.trim_galore.log" 2>&1
        else
            log "Trim skip (complete): ${sample}"
        fi
        require_file "${r1_trim}"; require_file "${r2_trim}"
        printf 'complete\n' > "${trim_marker}"
    else
        r1_trim="${trim_sample_dir}/$(basename "${r1}" | sed -E 's/\.(fastq|fq)(\.gz)?$/_trimmed.fq.gz/i')"
        r2_trim=""
        if [[ ! -s "${trim_marker}" || ! -s "${r1_trim}" ]]; then
            log "Trim Galore SE: ${sample}"
            trim_galore --quality "${TRIM_Q}" --length "${TRIM_MINLEN}" \
                --cores "${TRIM_CORES}" --gzip --output_dir "${trim_sample_dir}" "${r1}" \
                >"${LOG_DIR}/${sample}.trim_galore.log" 2>&1
        else
            log "Trim skip (complete): ${sample}"
        fi
        require_file "${r1_trim}"
        printf 'complete\n' > "${trim_marker}"
    fi

    if [[ -s "${map_marker}" && -s "${bam_out}" ]] && samtools quickcheck "${bam_out}"; then
        log "Alignment skip (complete): ${sample}"
        return
    fi
    log "Bismark ${LIBRARY_TYPE} alignment: ${sample}"
    # Bismark 0.25.1 rejects --basename together with --parallel >1. Use a
    # sample prefix, then normalize the single resulting BAM name below.
    local bm=(bismark --bowtie2 --genome "${BISMARK_GENOME_DIR}"
              --parallel "${BISMARK_PARALLEL}" --prefix "${sample}."
              --rg_tag --rg_id "${sample}" --rg_sample "${sample}"
              --output_dir "${map_sample_dir}" "${BM_DIRECTION[@]}")
    if [[ -n "${r2_trim}" ]]; then
        bm+=(-1 "${r1_trim}" -2 "${r2_trim}")
    else
        bm+=("${r1_trim}")
    fi
    "${bm[@]}" >"${LOG_DIR}/${sample}.bismark.log" 2>&1
    local produced_bams=()
    mapfile -t produced_bams < <(find "${map_sample_dir}" -maxdepth 1 -type f \
        -name '*.bam' ! -name "${sample}.bam" -print)
    [[ ${#produced_bams[@]} -eq 1 ]] || die \
        "${sample}: expected one Bismark BAM, found ${#produced_bams[@]} in ${map_sample_dir}"
    mv "${produced_bams[0]}" "${bam_out}"
    require_file "${bam_out}"
    samtools quickcheck -v "${bam_out}"
    printf 'complete\n' > "${map_marker}"
}

log "Trim and align samples from ${SAMPLES_TSV}"
for_each_sample process_sample
log "Trim/alignment complete: ${MAP_DIR}"
