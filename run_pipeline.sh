#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
stage="${1:-preprocess}"

case "${stage}" in
    validate) "${SCRIPT_DIR}/00_validate.sh" ;;
    qc) "${SCRIPT_DIR}/01_qc.sh" ;;
    preprocess)
        "${SCRIPT_DIR}/01_qc.sh"
        "${SCRIPT_DIR}/02_trim_decontam_map.sh"
        "${SCRIPT_DIR}/03_dedup.sh"
        "${SCRIPT_DIR}/04_meth_extractor.sh" mbias
        ;;
    extract)
        "${SCRIPT_DIR}/04_meth_extractor.sh" extract
        "${SCRIPT_DIR}/05_multiqc.sh"
        ;;
    *)
        printf 'Usage: %s {validate|qc|preprocess|extract}\n' "$0" >&2
        exit 2
        ;;
esac
