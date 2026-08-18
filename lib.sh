#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
if [[ -n "${TOOL_BIN}" ]]; then
    export PATH="${TOOL_BIN}:${PATH}"
fi

timestamp() { date '+%F %T'; }
log() { printf '[%s] %s\n' "$(timestamp)" "$*"; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

require_command() {
    local tool
    for tool in "$@"; do
        command -v "${tool}" >/dev/null 2>&1 || die "Required command not found: ${tool}"
    done
}

require_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }

validate_sample_name() {
    [[ "$1" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
        die "Unsafe sample name '$1' (allowed: letters, digits, dot, underscore, hyphen)"
}

# Calls a function once for each non-header samples.tsv row.
for_each_sample() {
    local callback="$1" sample r1 r2 extra line_no=0
    require_file "${SAMPLES_TSV}"
    while IFS=$'\t' read -r sample r1 r2 extra; do
        line_no=$((line_no + 1))
        sample="${sample%$'\r'}"; r1="${r1%$'\r'}"; r2="${r2%$'\r'}"
        [[ -z "${sample}" ]] && continue
        [[ ${line_no} -eq 1 && "${sample}" == "sample" ]] && continue
        [[ -z "${extra:-}" ]] || die "${SAMPLES_TSV}:${line_no}: expected exactly 3 tab-separated columns"
        validate_sample_name "${sample}"
        require_file "${r1}"
        [[ -z "${r2}" ]] || require_file "${r2}"
        "${callback}" "${sample}" "${r1}" "${r2}"
    done < "${SAMPLES_TSV}"
}

has_bismark_bt2_index() {
    local conv prefix suffix
    for conv in CT_conversion GA_conversion; do
        prefix="BS_${conv%%_*}"
        for suffix in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
            [[ -s "${BISMARK_GENOME_DIR}/Bisulfite_Genome/${conv}/${prefix}.${suffix}" ]] || return 1
        done
    done
}
