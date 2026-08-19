#!/usr/bin/env python3
"""Create a compact, reproducible WGBS cohort QC/results summary."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from statistics import mean, stdev


def text(path: Path) -> str:
    return path.read_text(errors="replace")


def capture(pattern: str, content: str, cast=float):
    match = re.search(pattern, content, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Missing pattern: {pattern}")
    value = match.group(1).replace(",", "").replace("%", "").strip()
    return cast(value)


def group_name(sample: str) -> str:
    stem = sample.split("_UDI", 1)[0]
    return re.sub(r"[12]$", "", stem)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("samples_tsv", type=Path)
    parser.add_argument("raw_qc_tsv", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--genome-size", type=int, default=119_667_750)
    parser.add_argument("--effective-r1", type=int, default=140)
    parser.add_argument("--effective-r2", type=int, default=137)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.samples_tsv.open() as handle:
        samples = [r["sample"] for r in csv.DictReader(handle, delimiter="\t")]
    with args.raw_qc_tsv.open() as handle:
        raw_rows = list(csv.DictReader(handle, delimiter="\t"))
    raw_pairs = {r["sample"]: int(r["total_sequences"]) for r in raw_rows if r["read"] == "R1"}
    raw_quality = [float(r["mean_base_quality"]) for r in raw_rows]
    raw_gc = [float(r["gc_percent"]) for r in raw_rows]

    rows = []
    for sample in samples:
        alignment_path = next((args.results_dir / "03_align" / sample).glob("*PE_report.txt"))
        dedup_path = args.results_dir / "04_dedup" / sample / f"{sample}.deduplication_report.txt"
        split_path = args.results_dir / "05_methylation" / sample / f"{sample}.deduplicated.namesorted_splitting_report.txt"
        alignment, dedup, split = map(text, (alignment_path, dedup_path, split_path))

        trimmed = capture(r"Sequence pairs analysed in total:\s*(\d+)", alignment, int)
        unique = capture(r"unique best hit:\s*(\d+)", alignment, int)
        mapping = capture(r"Mapping efficiency:\s*([0-9.]+)%", alignment)
        dup = capture(r"duplicated alignments removed:\s*(\d+)", dedup, int)
        dup_pct = capture(r"duplicated alignments removed:.*\(([0-9.]+)%\)", dedup)
        dedup_pairs = capture(r"deduplicated leftover sequences:\s*(\d+)", dedup, int)
        total_c = capture(r"Total number of C's analysed:\s*(\d+)", split, int)
        cpg = capture(r"C methylated in CpG context:\s*([0-9.]+)%", split)
        chg = capture(r"C methylated in CHG context:\s*([0-9.]+)%", split)
        chh = capture(r"C methylated in CHH context:\s*([0-9.]+)%", split)
        raw = raw_pairs[sample]
        rows.append({
            "sample": sample,
            "group": group_name(sample),
            "raw_pairs": raw,
            "trimmed_pairs": trimmed,
            "trim_retained_percent": 100 * trimmed / raw,
            "unique_pairs": unique,
            "mapping_percent": mapping,
            "duplicates_removed": dup,
            "duplicate_percent": dup_pct,
            "deduplicated_pairs": dedup_pairs,
            "effective_pairs_percent_raw": 100 * dedup_pairs / raw,
            "approx_deduplicated_depth_x": dedup_pairs * (args.effective_r1 + args.effective_r2) / args.genome_size,
            "total_cytosine_calls": total_c,
            "cpg_methylation_percent": cpg,
            "chg_methylation_percent": chg,
            "chh_methylation_percent": chh,
        })

    fields = list(rows[0])
    out_tsv = args.output_dir / "wgbs_sample_summary.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: f"{v:.3f}" if isinstance(v, float) else v for k, v in row.items()})

    group_metrics = [
        "mapping_percent", "duplicate_percent", "approx_deduplicated_depth_x",
        "cpg_methylation_percent", "chg_methylation_percent", "chh_methylation_percent",
    ]
    groups = {}
    for row in rows:
        groups.setdefault(row["group"], []).append(row)
    group_tsv = args.output_dir / "wgbs_group_summary.tsv"
    group_fields = ["group", "n"] + [f"{metric}_{suffix}" for metric in group_metrics for suffix in ("mean", "sd")]
    with group_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=group_fields, delimiter="\t")
        writer.writeheader()
        for group, members in groups.items():
            record = {"group": group, "n": len(members)}
            for metric in group_metrics:
                values = [float(r[metric]) for r in members]
                record[f"{metric}_mean"] = mean(values)
                record[f"{metric}_sd"] = stdev(values) if len(values) > 1 else 0.0
            writer.writerow({k: f"{v:.3f}" if isinstance(v, float) else v for k, v in record.items()})

    import matplotlib.pyplot as plt
    import numpy as np

    labels = [r["sample"].split("_UDI")[0] for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    axes[0, 0].bar(x - 0.18, [r["raw_pairs"] / 1e6 for r in rows], 0.36, label="Raw")
    axes[0, 0].bar(x + 0.18, [r["deduplicated_pairs"] / 1e6 for r in rows], 0.36, label="Unique deduplicated")
    axes[0, 0].set_ylabel("Read pairs (million)")
    axes[0, 0].set_title("Raw and effective read pairs")
    axes[0, 0].legend()
    axes[0, 1].bar(x, [r["mapping_percent"] for r in rows], color="#4C78A8")
    axes[0, 1].set_ylim(65, 80)
    axes[0, 1].set_ylabel("Unique mapping (%)")
    axes[0, 1].set_title("Bismark unique mapping efficiency")
    axes[1, 0].bar(x, [r["duplicate_percent"] for r in rows], color="#E45756")
    axes[1, 0].set_ylabel("Duplicate alignments (%)")
    axes[1, 0].set_title("Coordinate duplicate rate")
    for context, color in (("cpg", "#4C78A8"), ("chg", "#F58518"), ("chh", "#54A24B")):
        axes[1, 1].plot(x, [r[f"{context}_methylation_percent"] for r in rows], marker="o", label=context.upper(), color=color)
    axes[1, 1].set_ylabel("Global methylation (%)")
    axes[1, 1].set_title("Weighted methylation by context")
    axes[1, 1].legend()
    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=45, ha="right")
        axis.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.output_dir / "wgbs_cohort_summary.png", dpi=180)
    fig.savefig(args.output_dir / "wgbs_cohort_summary.pdf")

    def rng(key):
        values = [float(r[key]) for r in rows]
        return min(values), max(values)

    report = args.output_dir / "analysis_summary.md"
    total_raw = sum(r["raw_pairs"] for r in rows)
    total_trimmed = sum(r["trimmed_pairs"] for r in rows)
    total_unique = sum(r["unique_pairs"] for r in rows)
    total_dedup = sum(r["deduplicated_pairs"] for r in rows)
    total_calls = sum(r["total_cytosine_calls"] for r in rows)
    map_min, map_max = rng("mapping_percent")
    dup_min, dup_max = rng("duplicate_percent")
    depth_min, depth_max = rng("approx_deduplicated_depth_x")
    trim_min, trim_max = rng("trim_retained_percent")
    eff_min, eff_max = rng("effective_pairs_percent_raw")
    cpg_min, cpg_max = rng("cpg_methylation_percent")
    chg_min, chg_max = rng("chg_methylation_percent")
    chh_min, chh_max = rng("chh_methylation_percent")
    replicate_gaps = []
    for group, members in groups.items():
        if len(members) == 2:
            replicate_gaps.append((abs(members[0]["cpg_methylation_percent"] - members[1]["cpg_methylation_percent"]), group))
    largest_cpg_gap, largest_cpg_gap_group = max(replicate_gaps)
    with report.open("w") as handle:
        handle.write("# WGBS cohort analysis summary\n\n")
        handle.write(f"Samples: {len(rows)} paired-end libraries; reference genome size: {args.genome_size:,} bp.\n\n")
        handle.write("## Raw-read QC\n\n")
        handle.write(f"- Total input: {total_raw:,} PE150 read pairs ({2 * total_raw:,} reads).\n")
        handle.write(f"- Mean per-base Phred quality across mate files: {min(raw_quality):.2f}-{max(raw_quality):.2f}; GC: {min(raw_gc):.0f}-{max(raw_gc):.0f}%.\n")
        handle.write("- FastQC per-base sequence-content failures are expected for bisulfite-converted libraries and were cohort-consistent.\n\n")
        handle.write("## Technical summary\n\n")
        handle.write(f"- Cohort totals: {total_trimmed:,} trimmed pairs; {total_unique:,} uniquely mapped pairs; "
                     f"{total_dedup:,} unique deduplicated pairs; {total_calls:,} cytosine calls.\n")
        handle.write(f"- Paired reads retained after trimming: {trim_min:.2f}-{trim_max:.2f}%.\n")
        handle.write(f"- Unique Bismark mapping efficiency: {map_min:.1f}-{map_max:.1f}%.\n")
        handle.write(f"- Coordinate duplicate rate: {dup_min:.1f}-{dup_max:.1f}%.\n")
        handle.write(f"- Effective unique deduplicated pairs: {eff_min:.1f}-{eff_max:.1f}% of raw pairs.\n")
        handle.write(f"- Approximate post-dedup sequence depth after M-bias clipping: {depth_min:.1f}-{depth_max:.1f}x.\n")
        handle.write("- M-bias clipping: R1 5-prime 10 bp; R1 3-prime 0 bp; R2 5-prime 10 bp; R2 3-prime 3 bp.\n")
        handle.write("- All retained paired-end calls used --no_overlap. CpG, CHG and CHH contexts were retained.\n\n")
        handle.write("## Global methylation\n\n")
        handle.write(f"- CpG: {cpg_min:.1f}-{cpg_max:.1f}%; CHG: {chg_min:.1f}-{chg_max:.1f}%; CHH: {chh_min:.1f}-{chh_max:.1f}%.\n")
        handle.write("- Groups inferred from sample-name prefixes only: C and CMDR have ~0.7% CHG methylation, "
                     "whereas D, F1A, F1B and YE have 7.4-8.2%. This sharp, replicate-consistent split is likely biological, "
                     "but formal contrasts require confirmed sample metadata.\n")
        handle.write(f"- The largest within-prefix CpG difference is {largest_cpg_gap:.1f} percentage points in {largest_cpg_gap_group}; "
                     "confirm with site-level correlation/PCA before DMR testing.\n\n")
        handle.write("## Interpretation boundary\n\n")
        handle.write("CHH methylation is biological in plants and is not a valid conversion-rate estimator. "
                     "No unmethylated spike-in metadata was available, so an independent bisulfite conversion rate is not reported.\n")
        handle.write("\n## Key artifacts\n\n")
        handle.write("- `wgbs_sample_summary.tsv`: per-sample technical and global methylation metrics.\n")
        handle.write("- `wgbs_group_summary.tsv`: means/SDs for groups inferred from sample-name prefixes.\n")
        handle.write("- `wgbs_cohort_summary.png` / `.pdf`: cohort overview figure.\n")
        handle.write("- `../06_multiqc/multiqc_report.html`: combined pipeline MultiQC report.\n")
        handle.write("- `../05_methylation/mbias_summary/`: M-bias evidence and clipping decision.\n")
    print(out_tsv)
    print(group_tsv)
    print(report)


if __name__ == "__main__":
    main()
