#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib.sh"

require_command fastqc multiqc trim_galore cutadapt bowtie2 bismark \
    bismark_genome_preparation deduplicate_bismark \
    bismark_methylation_extractor samtools
[[ -n "${GENOME_FA}" ]] || die "GENOME_FA is not configured"
[[ -n "${BISMARK_GENOME_DIR}" ]] || die "BISMARK_GENOME_DIR is not configured"
require_file "${GENOME_FA}"
has_bismark_bt2_index || die "Incomplete Bismark Bowtie2 index: ${BISMARK_GENOME_DIR}"

case "${LIBRARY_TYPE}" in
    directional|non_directional|pbat) ;;
    *) die "LIBRARY_TYPE must be directional, non_directional, or pbat" ;;
esac

sample_count=0
pe_count=0
declare -A seen_samples=()
validate_row() {
    local sample="$1" r1="$2" r2="$3"
    [[ -z "${seen_samples[${sample}]:-}" ]] || die "Duplicate sample in ${SAMPLES_TSV}: ${sample}"
    seen_samples["${sample}"]=1
    sample_count=$((sample_count + 1))
    [[ -z "${r2}" ]] || pe_count=$((pe_count + 1))
}
for_each_sample validate_row
(( sample_count > 0 )) || die "No samples found in ${SAMPLES_TSV}"

log "Validation passed"
printf 'Samples: %d (%d paired-end)\n' "${sample_count}" "${pe_count}"
printf 'Reference: %s\n' "${GENOME_FA}"
printf 'Bismark index: %s\n' "${BISMARK_GENOME_DIR}"
printf 'Output: %s\n' "${OUTDIR}"
printf 'Library type: %s\n' "${LIBRARY_TYPE}"
bismark --version 2>&1 | grep 'Bismark Version' | sed 's/^ */Tool: /'
multiqc --version
samtools --version | head -n 1
