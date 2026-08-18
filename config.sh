#!/usr/bin/env bash
# WGBS pipeline configuration. Every value can be overridden from the environment.

PROJECT_DIR="${PROJECT_DIR:-${PWD}}"
SAMPLES_TSV="${SAMPLES_TSV:-${PROJECT_DIR}/samples.tsv}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/wgbs_results}"
TOOL_BIN="${TOOL_BIN:-${CONDA_PREFIX:+${CONDA_PREFIX}/bin}}"

# Required project-specific paths. Export them in the shell or set them here.
GENOME_FA="${GENOME_FA:-}"
BISMARK_GENOME_DIR="${BISMARK_GENOME_DIR:-}"

# Resource limits. Bismark --parallel multiplies CPU and RAM. Trim Galore
# --cores 4 may use ~15 processes; extractor --parallel 4 uses ~12 cores.
THREADS="${THREADS:-24}"
FASTQC_THREADS="${FASTQC_THREADS:-12}"
TRIM_CORES="${TRIM_CORES:-4}"
BISMARK_PARALLEL="${BISMARK_PARALLEL:-4}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-8}"
EXTRACT_PARALLEL="${EXTRACT_PARALLEL:-4}"
EXTRACT_BUFFER="${EXTRACT_BUFFER:-10G}"
TRIM_Q="${TRIM_Q:-20}"
TRIM_MINLEN="${TRIM_MINLEN:-20}"

# directional | non_directional | pbat
LIBRARY_TYPE="${LIBRARY_TYPE:-directional}"

# Ordinary Bowtie2 is generally inappropriate for bisulfite-converted reads.
ENABLE_DECONTAM="${ENABLE_DECONTAM:-false}"

# First run 04_meth_extractor.sh mbias. Inspect the reports, set these four
# values, then set MBIAS_REVIEWED=true before making final calls.
MBIAS_REVIEWED="${MBIAS_REVIEWED:-false}"
IGNORE_R1_5="${IGNORE_R1_5:-0}"
IGNORE_R1_3="${IGNORE_R1_3:-0}"
IGNORE_R2_5="${IGNORE_R2_5:-0}"
IGNORE_R2_3="${IGNORE_R2_3:-0}"

# WGBS/EM-seq: true. RRBS: false (coordinate dedup is invalid for RRBS).
ENABLE_DEDUP="${ENABLE_DEDUP:-true}"
