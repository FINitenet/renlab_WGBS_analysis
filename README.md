# renlab WGBS analysis

A restartable WGBS workflow built around FastQC, MultiQC, Trim Galore and
Bismark. The preferred entry point is:

- `wgbs_oneclick.py`: FASTQ to QC, mapping, methylation calls, summaries and
  all publication figures in one restartable command.

The smaller Python/shell programs remain internal, independently testable
stages. `wgbs_qc.py` and `run_pipeline.sh` are retained for compatibility.

The workflow supports directional, non-directional and PBAT libraries. WGBS
coordinate deduplication is enabled by default and must be disabled for RRBS.

## Requirements

- Python 3.8+
- numpy, pandas, matplotlib, seaborn and pyBigWig
- deepTools `computeMatrix`, pigz and awk
- FastQC
- MultiQC
- Trim Galore / Cutadapt
- Bismark 0.25+
- Bowtie2
- samtools

Activate an environment containing these tools, or set `TOOL_BIN` to its `bin`
directory.

## One-command complete analysis

The analysis directory should contain `1_rawdata_merged/` (preferred) or
`1_rawdata/`. The program searches upward for
`reference/TAIR10_plus_transgene_bismark_bt2`, uses its combined FASTA/index,
and writes new-reference results to `wgbs_joint_results/` so old TAIR10-only
completion markers cannot be reused accidentally.

Run the complete workflow with one command:

```bash
./wgbs_oneclick.py \
  /path/to/analysis \
  --tool-bin /path/to/conda/env/bin
```

This single command performs:

1. FASTQ discovery/lane staging, pair validation, FastQC and raw MultiQC.
2. Configuration/index validation, Trim Galore and Bismark alignment.
3. WGBS deduplication, coordinate/name sorting and M-bias reports.
4. Final extraction with the reviewed cohort clipping values
   `R1=10/0 bp, R2=10/3 bp`, then final MultiQC and cohort summaries.
5. Gene/TE CG/CHG/CHH metaprofiles, the requested mC/mCG/mCHG/mCHH boxplot
   panel, and whole-transgene locus/feature-bar figures.

Every expensive stage is restartable. The file `.oneclick_config.json` binds an
output directory to the FASTA checksum, Bismark genome and sample sheet; the
program refuses to reuse results if those analysis-defining inputs change.
Inspect a complete command plan without launching computation:

```bash
./wgbs_oneclick.py /path/to/analysis --dry-run
```

Resume or run only a portion with `--from-stage` / `--to-stage`, for example:

```bash
./wgbs_oneclick.py /path/to/analysis --from-stage figures
```

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

Aggregate a cohort's Bismark M-bias reports into a weighted TSV and PNG/PDF:

```bash
./summarize_mbias.py /path/to/wgbs_results/05_methylation \
  -o /path/to/wgbs_results/05_methylation/mbias_summary
```

After extraction, build per-sample/group TSVs, a cohort figure and a Markdown
summary (the default genome size is TAIR10):

```bash
./summarize_results.py /path/to/wgbs_results /path/to/samples.tsv \
  /path/to/wgbs_qc/02_qc/qc_summary.tsv \
  -o /path/to/wgbs_results/07_summary
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
- Deduplicated BAMs are written in two forms: coordinate-sorted/indexed for
  browsing and interval tools, and query-name-sorted for paired-end Bismark
  methylation extraction so mates remain adjacent.
- M-bias clipping is a required review gate, not a fixed arbitrary number.
- RRBS must set `ENABLE_DEDUP=false`; coordinate deduplication is invalid for
  restriction-enzyme-defined fragments.

## DMR analysis

DMR calling is intentionally not hard-coded because it requires an explicit
biological design, contrasts, covariates and raw methylated/unmethylated counts.
For small replicated WGBS studies, DSS with smoothing is a reasonable model;
with sufficient replication, dmrseq provides selection-aware region FDR.
Single-CpG DMCs must not be reported as DMRs.

## Reference-style downstream figures

`wgbs_reference_figures.py` turns completed Bismark CX reports into reusable
context-specific bigWig tracks and count-weighted windows, then exports the
following as both PNG and editable-text PDF:

- CG/CHG/CHH profiles across scaled genes and transposable elements, including
  2 kb flanks;
- mC/mCG/mCHG/mCHH distributions across 1 kb genomic windows;
- CG/CHG/CHH profiles for any requested locus.

```bash
./wgbs_reference_figures.py prepare \
  --methylation-dir /path/to/wgbs_results/05_methylation \
  --gff /path/to/TAIR10_GFF3_genes_transposons.gff \
  --outdir /path/to/wgbs_results/08_reference_figures

./wgbs_reference_figures.py locus \
  --outdir /path/to/wgbs_results/08_reference_figures \
  --region chr1:100000-110000
```

The default sample groups are inferred by removing the replicate number and
UDI suffix (for example, `F1A1_UDI9077` becomes `F1A`). Review the emitted
`sample_metadata.inferred.tsv` before biological interpretation. Window
boxplots require both replicates and use summed methylated/unmethylated read
counts; they deliberately omit significance stars because genomic windows are
not biological replicates.
