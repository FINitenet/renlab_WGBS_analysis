# renlab WGBS analysis

A restartable WGBS workflow built around FastQC, MultiQC, Trim Galore and
Bismark. It provides two entry points:

- `wgbs_qc.py`: input FASTQ folder to sample-level raw-read QC figures.
- `run_pipeline.sh`: trimming, alignment, deduplication, M-bias review,
  methylation extraction and final MultiQC aggregation.

The workflow supports directional, non-directional and PBAT libraries. WGBS
coordinate deduplication is enabled by default and must be disabled for RRBS.

## Requirements

- Python 3.8+
- matplotlib (only for the static PNG/PDF overview)
- FastQC
- MultiQC
- Trim Galore / Cutadapt
- Bismark 0.25+
- Bowtie2
- samtools

Activate an environment containing these tools, or set `TOOL_BIN` to its `bin`
directory.

## Input folder to QC figures

Only an input FASTQ directory is required. The default output directory is
`<input parent>/wgbs_qc`:

```bash
./wgbs_qc.py /path/to/input_fastq
```

Choose the output directory or inspect sample/lane inference without writing:

```bash
./wgbs_qc.py /path/to/input_fastq -o /path/to/wgbs_qc --threads 12
./wgbs_qc.py /path/to/input_fastq --dry-run
```

The Python workflow:

1. Recursively discovers paired `.fastq[.gz]` / `.fq[.gz]` files.
2. Removes `_L001`-style lane tokens to form sample IDs.
3. Links single-lane gzip inputs and concatenates true multi-lane samples.
4. Checks leading mate names and verifies complete R1/R2 counts via FastQC.
5. Runs FastQC and a FastQC-scoped MultiQC report.
6. Writes:
   - `samples.tsv`
   - `02_qc/qc_summary.tsv`
   - `02_qc/fastqc_overview.png`
   - `02_qc/fastqc_overview.pdf`
   - `02_qc/multiqc/wgbs_qc_multiqc.html`

State files make reruns incremental. Use `--force` to rebuild outputs.

> FastQC traffic lights are generic WGS heuristics. WGBS commonly fails the
> per-base sequence-content module because bisulfite conversion depletes C and
> enriches T. Interpret the plot and cohort outliers rather than that flag alone.

## Full WGBS workflow

The shell workflow consumes a three-column, tab-separated sample sheet:

```text
sample	R1	R2
sample_1	/path/sample_1.R1.fastq.gz	/path/sample_1.R2.fastq.gz
```

The QC Python entry point produces this compatible `samples.tsv` directly;
lane counts are retained in `qc_summary.tsv`.

Configure project paths through environment variables:

```bash
export PROJECT_DIR=/path/to/project
export SAMPLES_TSV=/path/to/project/samples.tsv
export OUTDIR=/path/to/project/wgbs_results
export GENOME_FA=/path/to/reference/genome.fa
export BISMARK_GENOME_DIR=/path/to/bismark_genome
export LIBRARY_TYPE=directional   # directional | non_directional | pbat
export TOOL_BIN=/path/to/conda/env/bin   # optional when tools are already on PATH
```

Validate without processing reads, then run through M-bias diagnostics:

```bash
./run_pipeline.sh validate
./run_pipeline.sh preprocess
```

Inspect each `05_methylation/<sample>/*M-bias*` report. Set the four end-clipping
values in `config.sh` (or export them) and explicitly approve the QC gate:

```bash
export IGNORE_R1_5=0
export IGNORE_R1_3=0
export IGNORE_R2_5=0
export IGNORE_R2_3=0
export MBIAS_REVIEWED=true
./run_pipeline.sh extract
```

The extractor retains CpG, CHG and CHH contexts for plant methylomes. In plants,
CHH is biological and cannot be used as a bisulfite conversion proxy; use an
unmethylated spike-in or the library's independent conversion control.

## Important design choices

- Bismark uses its default `-N 0` seed policy instead of relaxing to `-N 1`.
- Bismark index completeness is checked using both converted Bowtie2 indexes.
- Read groups use `--rg_tag --rg_id --rg_sample` syntax supported by Bismark
  0.25.x.
- Bismark 0.25.x cannot combine `--basename` with multicore mapping, so output
  uses a prefix and is normalized only after successful completion.
- Ordinary Bowtie2 contaminant filtering is rejected because it is not
  bisulfite-aware.
- Paired overlap is counted once during methylation extraction.
- M-bias clipping is a required review gate, not a fixed arbitrary number.
- RRBS must set `ENABLE_DEDUP=false`; coordinate deduplication is invalid for
  restriction-enzyme-defined fragments.

## DMR analysis

DMR calling is intentionally not hard-coded because it requires an explicit
biological design, contrasts, covariates and raw methylated/unmethylated counts.
For small replicated WGBS studies, DSS with smoothing is a reasonable model;
with sufficient replication, dmrseq provides selection-aware region FDR.
Single-CpG DMCs must not be reported as DMRs.
