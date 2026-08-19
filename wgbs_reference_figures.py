#!/usr/bin/env python3
"""Create publication-style downstream WGBS figures from Bismark CX reports.

The script prepares context-specific bigWig tracks and count-weighted genomic
windows, then draws (1) gene/TE metaprofiles, (2) regional methylation
distribution boxplots, and (3) a locus view when a genomic interval is given.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-wgbs-figures")

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns

CONTEXTS = ("CG", "CHG", "CHH")
COLORS = {
    "C": "#D55E00",
    "CMDR": "#0072B2",
    "D": "#009E73",
    "F1A": "#CC79A7",
    "F1B": "#E69F00",
    "YE": "#56B4E9",
}


def natural_key(text: str):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", text)]


def normalize_chrom(chrom: str) -> str:
    return chrom.lower()


def infer_group(sample: str) -> str:
    core = sample.split("_UDI", 1)[0]
    return re.sub(r"\d+$", "", core)


def parse_gff(gff: Path, outdir: Path):
    chrom_sizes = {}
    genes, tes = [], []
    with gff.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, feature = normalize_chrom(fields[0]), fields[2]
            start, end, strand = int(fields[3]), int(fields[4]), fields[6]
            attrs = fields[8]
            if feature == "chromosome":
                chrom_sizes[chrom] = max(chrom_sizes.get(chrom, 0), end)
            if feature not in {"gene", "transposable_element"}:
                continue
            # deepTools scale-regions cannot scale features shorter than one bin.
            if end - start + 1 < 100:
                continue
            match = re.search(r"(?:^|;)ID=([^;]+)", attrs)
            name = match.group(1) if match else f"{feature}_{len(genes) + len(tes) + 1}"
            row = (chrom, start - 1, end, name, 0, strand)
            (genes if feature == "gene" else tes).append(row)

    if not chrom_sizes:
        for row in genes + tes:
            chrom_sizes[row[0]] = max(chrom_sizes.get(row[0], 0), row[2])
    order = {c: i for i, c in enumerate(chrom_sizes)}
    for name, rows in (("genes", genes), ("transposable_elements", tes)):
        rows.sort(key=lambda x: (order.get(x[0], 10**9), x[1], x[2]))
        with (outdir / f"{name}.bed").open("w") as handle:
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")
    with (outdir / "chrom.sizes").open("w") as handle:
        for chrom, length in chrom_sizes.items():
            handle.write(f"{chrom}\t{length}\n")
    return chrom_sizes


def bedgraph_to_bigwig(bedgraph: Path, bigwig: Path, chrom_sizes):
    bw = pyBigWig.open(str(bigwig), "w")
    bw.addHeader(list(chrom_sizes.items()))
    chroms, starts, ends, values = [], [], [], []
    with bedgraph.open() as handle:
        for line in handle:
            chrom, start, end, value = line.rstrip("\n").split("\t")
            chroms.append(chrom)
            starts.append(int(start))
            ends.append(int(end))
            values.append(float(value))
            if len(chroms) >= 500000:
                bw.addEntries(chroms, starts, ends=ends, values=values)
                chroms, starts, ends, values = [], [], [], []
    if chroms:
        bw.addEntries(chroms, starts, ends=ends, values=values)
    bw.close()


def prepare_sample(report_s: str, outdir_s: str, chrom_sizes, window_size: int,
                   min_site_coverage: int):
    report, outdir = Path(report_s), Path(outdir_s)
    sample = report.parent.name
    track_dir = outdir / "tracks"
    window_dir = outdir / "windows"
    track_dir.mkdir(parents=True, exist_ok=True)
    window_dir.mkdir(parents=True, exist_ok=True)
    window_path = window_dir / f"{sample}.windows.tsv.gz"
    temp_dir = outdir / "tmp" / sample
    temp_dir.mkdir(parents=True, exist_ok=True)
    prefix = temp_dir / sample
    window_tsv = temp_dir / f"{sample}.windows.tsv"
    for path in [window_tsv, *(Path(f"{prefix}.{ctx}.bedGraph") for ctx in CONTEXTS)]:
        path.unlink(missing_ok=True)
    awk_script = Path(__file__).with_name("cx_to_bedgraphs.awk")
    pigz = subprocess.Popen(["pigz", "-dc", str(report)], stdout=subprocess.PIPE)
    awk_cmd = [
        "awk", "-v", f"prefix={prefix}", "-v", f"window_file={window_tsv}",
        "-v", f"sizes_file={outdir / 'chrom.sizes'}", "-v", f"sample={sample}",
        "-v", f"group={infer_group(sample)}", "-v", f"window_size={window_size}",
        "-v", f"min_site_coverage={min_site_coverage}", "-f", str(awk_script),
    ]
    awk = subprocess.run(awk_cmd, stdin=pigz.stdout)
    pigz.stdout.close()
    pigz_rc = pigz.wait()
    if pigz_rc != 0 or awk.returncode != 0:
        raise RuntimeError(f"CX conversion failed for {sample}: pigz={pigz_rc}, awk={awk.returncode}")

    counts = {}
    for ctx in CONTEXTS:
        bedgraph = Path(f"{prefix}.{ctx}.bedGraph")
        bedgraph_to_bigwig(bedgraph, track_dir / f"{sample}.{ctx}.bw", chrom_sizes)
        with bedgraph.open() as handle:
            counts[ctx] = sum(1 for _ in handle)
        bedgraph.unlink()
    with window_tsv.open("rb") as source, gzip.open(window_path, "wb") as target:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            target.write(chunk)
    window_tsv.unlink()
    temp_dir.rmdir()
    return sample, counts


def discover_reports(methylation_dir: Path):
    reports = sorted(methylation_dir.glob("*/*.CX_report.txt.gz"), key=lambda p: natural_key(p.parent.name))
    if not reports:
        raise SystemExit(f"No *.CX_report.txt.gz files found under {methylation_dir}")
    return reports


def write_metadata(reports, path: Path):
    rows = []
    replicate_no = defaultdict(int)
    for report in reports:
        sample = report.parent.name
        group = infer_group(sample)
        replicate_no[group] += 1
        rows.append((sample, group, replicate_no[group], "inferred_from_sample_prefix"))
    pd.DataFrame(rows, columns=["sample", "group", "replicate", "metadata_status"]).to_csv(
        path, sep="\t", index=False
    )


def run_compute_matrix(outdir: Path, metadata: pd.DataFrame, threads: int):
    matrix_dir = outdir / "matrices"
    matrix_dir.mkdir(exist_ok=True)
    samples = metadata["sample"].tolist()
    for ctx in CONTEXTS:
        tracks = [str(outdir / "tracks" / f"{sample}.{ctx}.bw") for sample in samples]
        matrix = matrix_dir / f"genes_TEs.{ctx}.matrix.gz"
        cmd = [
            "computeMatrix", "scale-regions", "-S", *tracks,
            "-R", str(outdir / "genes.bed"), str(outdir / "transposable_elements.bed"),
            "--beforeRegionStartLength", "2000", "--regionBodyLength", "4000",
            "--afterRegionStartLength", "2000", "--binSize", "100",
            "--averageTypeBins", "mean", "--numberOfProcessors", str(threads),
            "--samplesLabel", *samples,
            "-o", str(matrix),
        ]
        subprocess.run(cmd, check=True)


def read_matrix_profiles(matrix: Path, metadata: pd.DataFrame):
    with gzip.open(matrix, "rt") as handle:
        header = json.loads(handle.readline()[1:])
        values = np.loadtxt(handle, usecols=range(6, 6 + header["sample_boundaries"][-1]))
    sample_bounds = header["sample_boundaries"]
    region_bounds = header["group_boundaries"]
    result = {}
    canonical_regions = ("Genes", "TEs")
    for region_i, _region in enumerate(header["group_labels"]):
        region = canonical_regions[region_i]
        block = values[region_bounds[region_i]:region_bounds[region_i + 1], :]
        for sample_i, sample in enumerate(metadata["sample"]):
            profile = np.nanmean(block[:, sample_bounds[sample_i]:sample_bounds[sample_i + 1]], axis=0)
            result[(region, sample)] = profile
    return result


def save_figure(fig, stem: Path):
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def plot_metaprofiles(outdir: Path, metadata: pd.DataFrame):
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 5.8), constrained_layout=True)
    for col, ctx in enumerate(CONTEXTS):
        profiles = read_matrix_profiles(outdir / "matrices" / f"genes_TEs.{ctx}.matrix.gz", metadata)
        for row, region in enumerate(("Genes", "TEs")):
            ax = axes[row, col]
            for group, group_df in metadata.groupby("group", sort=False):
                vals = np.vstack([profiles[(region, sample)] for sample in group_df["sample"]])
                mean = np.nanmean(vals, axis=0)
                x = np.arange(len(mean))
                ax.plot(x, mean, lw=1.5, color=COLORS.get(group), label=group)
                if len(vals) > 1:
                    sem = np.nanstd(vals, axis=0, ddof=1) / np.sqrt(len(vals))
                    ax.fill_between(x, mean - sem, mean + sem, color=COLORS.get(group), alpha=0.15, lw=0)
            ax.axvline(20, color="#777777", ls="--", lw=0.7)
            ax.axvline(60, color="#777777", ls="--", lw=0.7)
            ax.set_xticks([0, 20, 60, 79], ["-2 kb", "TSS", "TES", "+2 kb"])
            ax.set_title(f"{region}: {ctx} methylation")
            ax.set_ylabel("Methylation level")
            sns.despine(ax=ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=len(labels), frameon=False)
    save_figure(fig, outdir / "gene_TE_methylation_profiles")


def load_group_windows(outdir: Path, metadata: pd.DataFrame, min_coverage: int, min_sites: int):
    usecols = ["sample", "group", "chrom", "start", "end", "context", "methylated", "unmethylated", "covered_sites"]
    frames = [pd.read_csv(outdir / "windows" / f"{sample}.windows.tsv.gz", sep="\t", usecols=usecols)
              for sample in metadata["sample"]]
    all_data = pd.concat(frames, ignore_index=True)
    keys = ["group", "chrom", "start", "end", "context"]
    agg = all_data.groupby(keys, observed=True, sort=False).agg(
        methylated=("methylated", "sum"), unmethylated=("unmethylated", "sum"),
        covered_sites=("covered_sites", "sum"), replicates=("sample", "nunique")
    ).reset_index()
    expected = metadata.groupby("group")["sample"].nunique().to_dict()
    agg = agg[agg.apply(lambda r: r["replicates"] == expected[r["group"]], axis=1)]
    agg["coverage"] = agg["methylated"] + agg["unmethylated"]
    agg = agg[(agg["coverage"] >= min_coverage) & (agg["covered_sites"] >= min_sites)]
    agg["methylation_percent"] = 100 * agg["methylated"] / agg["coverage"]
    return agg


def plot_boxplots(outdir: Path, metadata: pd.DataFrame, min_coverage: int, min_sites: int):
    data = load_group_windows(outdir, metadata, min_coverage, min_sites)
    groups = metadata["group"].drop_duplicates().tolist()
    palette = {g: COLORS.get(g, "#777777") for g in groups}
    fig, axes = plt.subplots(4, 1, figsize=(7.2, 8.6), sharex=True, constrained_layout=True)
    for ax, ctx in zip(axes, ("mC", *CONTEXTS)):
        subset = data[data["context"] == ctx]
        sns.boxplot(data=subset, x="group", y="methylation_percent", order=groups,
                    hue="group", palette=palette, legend=False, showfliers=False,
                    width=0.68, linewidth=0.8, ax=ax)
        display_context = ctx if ctx == "mC" else f"m{ctx}"
        ax.set_ylabel(f"{display_context}\nmethylation (%)")
        ax.set_xlabel("")
        ax.text(0.995, 0.95, f"n={len(subset):,} windows", transform=ax.transAxes,
                ha="right", va="top", fontsize=7, color="#555555")
        sns.despine(ax=ax)
    axes[-1].set_xlabel("Sample group")
    save_figure(fig, outdir / "genome_1kb_window_methylation_boxplots")
    data.to_csv(outdir / "genome_1kb_window_group_methylation.tsv.gz", sep="\t", index=False)


def parse_region(text: str):
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", text.replace(",", ""))
    if not match:
        raise argparse.ArgumentTypeError("region must look like chr1:1000-5000")
    return normalize_chrom(match.group(1)), int(match.group(2)) - 1, int(match.group(3))


def plot_locus(outdir: Path, metadata: pd.DataFrame, region_text: str, bins: int):
    chrom, start, end = parse_region(region_text)
    edges = np.linspace(start, end, bins + 1, dtype=int)
    centers = (edges[:-1] + edges[1:]) / 2
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 6.0), sharex=True, constrained_layout=True)
    for ax, ctx in zip(axes, CONTEXTS):
        for group, group_df in metadata.groupby("group", sort=False):
            replicate_profiles = []
            for sample in group_df["sample"]:
                bw = pyBigWig.open(str(outdir / "tracks" / f"{sample}.{ctx}.bw"))
                values = np.asarray(bw.stats(chrom, start, end, nBins=bins, type="mean"), dtype=float)
                bw.close()
                replicate_profiles.append(values)
            profile = np.nanmean(np.vstack(replicate_profiles), axis=0)
            ax.plot(centers, profile, lw=1.4, color=COLORS.get(group), label=group)
        ax.set_ylabel(f"{ctx}\nmethylation")
        ax.set_ylim(bottom=0)
        sns.despine(ax=ax)
    axes[-1].set_xlabel(f"{chrom} coordinate (bp)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=len(labels), frameon=False)
    save_figure(fig, outdir / f"locus_{chrom}_{start+1}_{end}")


def command_prepare(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chrom_sizes = parse_gff(Path(args.gff), outdir)
    reports = discover_reports(Path(args.methylation_dir))
    metadata_path = outdir / "sample_metadata.inferred.tsv"
    write_metadata(reports, metadata_path)
    futures = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for report in reports:
            futures.append(pool.submit(prepare_sample, str(report), str(outdir), chrom_sizes,
                                       args.window_size, args.min_site_coverage))
        for future in as_completed(futures):
            sample, counts = future.result()
            print(f"prepared {sample}: {counts}", flush=True)
    metadata = pd.read_csv(metadata_path, sep="\t")
    run_compute_matrix(outdir, metadata, args.threads)
    plot_metaprofiles(outdir, metadata)
    plot_boxplots(outdir, metadata, args.min_window_coverage, args.min_window_sites)
    params = vars(args).copy()
    params.pop("func", None)
    params["methylation_dir"] = str(Path(args.methylation_dir).resolve())
    params["gff"] = str(Path(args.gff).resolve())
    with (outdir / "figure_parameters.json").open("w") as handle:
        json.dump(params, handle, indent=2, sort_keys=True)


def command_locus(args):
    outdir = Path(args.outdir)
    metadata = pd.read_csv(outdir / "sample_metadata.inferred.tsv", sep="\t")
    plot_locus(outdir, metadata, args.region, args.bins)


def command_plot(args):
    outdir = Path(args.outdir)
    parse_gff(Path(args.gff), outdir)
    metadata = pd.read_csv(outdir / "sample_metadata.inferred.tsv", sep="\t")
    run_compute_matrix(outdir, metadata, args.threads)
    plot_metaprofiles(outdir, metadata)
    plot_boxplots(outdir, metadata, args.min_window_coverage, args.min_window_sites)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare", help="prepare tracks/windows and draw cohort figures")
    prep.add_argument("--methylation-dir", required=True)
    prep.add_argument("--gff", required=True)
    prep.add_argument("--outdir", required=True)
    prep.add_argument("--workers", type=int, default=2, help="CX reports processed concurrently")
    prep.add_argument("--threads", type=int, default=8, help="computeMatrix worker count")
    prep.add_argument("--window-size", type=int, default=1000)
    prep.add_argument("--min-site-coverage", type=int, default=1)
    prep.add_argument("--min-window-coverage", type=int, default=20)
    prep.add_argument("--min-window-sites", type=int, default=5)
    prep.set_defaults(func=command_prepare)
    plot = sub.add_parser("plot", help="reuse prepared tracks/windows and redraw cohort figures")
    plot.add_argument("--gff", required=True)
    plot.add_argument("--outdir", required=True)
    plot.add_argument("--threads", type=int, default=8)
    plot.add_argument("--min-window-coverage", type=int, default=20)
    plot.add_argument("--min-window-sites", type=int, default=5)
    plot.set_defaults(func=command_plot)
    locus = sub.add_parser("locus", help="draw a context-specific locus profile")
    locus.add_argument("--outdir", required=True)
    locus.add_argument("--region", required=True, help="e.g. chr1:10000-20000")
    locus.add_argument("--bins", type=int, default=100)
    locus.set_defaults(func=command_locus)
    return parser


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
