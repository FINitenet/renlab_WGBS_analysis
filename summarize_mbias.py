#!/usr/bin/env python3
"""Aggregate Bismark M-bias text reports across a WGBS cohort."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SECTION_RE = re.compile(r"^(CpG|CHG|CHH) context \((R[12])\)$")


def parse_report(path: Path):
    context = read = None
    sample = path.parent.name
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            match = SECTION_RE.match(line)
            if match:
                context, read = match.groups()
                continue
            if not context or not line or line.startswith(("=", "position")):
                continue
            fields = line.split("\t")
            if len(fields) != 5 or not fields[0].isdigit():
                continue
            yield {
                "sample": sample,
                "context": context,
                "read": read,
                "position": int(fields[0]),
                "methylated": int(fields[1]),
                "unmethylated": int(fields[2]),
                "coverage": int(fields[4]),
            }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    args = parser.parse_args()

    reports = sorted(args.input_dir.glob("*/*M-bias.txt"))
    if not reports:
        raise SystemExit(f"No M-bias reports found below {args.input_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    totals = defaultdict(lambda: [0, 0, 0])
    rows = []
    for report in reports:
        for row in parse_report(report):
            rows.append(row)
            key = (row["context"], row["read"], row["position"])
            totals[key][0] += row["methylated"]
            totals[key][1] += row["unmethylated"]
            totals[key][2] += row["coverage"]

    tsv = args.output_dir / "mbias_weighted_by_position.tsv"
    with tsv.open("w", newline="") as handle:
        fields = ["context", "read", "position", "methylated", "unmethylated", "coverage", "methylation_percent"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(totals, key=lambda x: (x[1], x[0], x[2])):
            methylated, unmethylated, coverage = totals[key]
            pct = 100.0 * methylated / (methylated + unmethylated) if methylated + unmethylated else 0.0
            writer.writerow(dict(zip(fields, [*key, methylated, unmethylated, coverage, f"{pct:.6f}"])))

    import matplotlib.pyplot as plt

    colors = {"CpG": "#276FBF", "CHG": "#F28E2B", "CHH": "#2CA02C"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    for axis, read in zip(axes, ("R1", "R2")):
        for context in ("CpG", "CHG", "CHH"):
            points = []
            for (ctx, rd, pos), values in totals.items():
                if ctx == context and rd == read:
                    m, u, _ = values
                    points.append((pos, 100.0 * m / (m + u)))
            points.sort()
            axis.plot([p for p, _ in points], [v for _, v in points], label=context, color=colors[context], linewidth=1.5)
        axis.set_title(f"Cohort-weighted M-bias ({read})")
        axis.set_xlabel("Read position (bp)")
        axis.set_ylabel("Methylation (%)")
        axis.set_xlim(1, 150)
        axis.grid(alpha=0.2)
        axis.legend()
    fig.tight_layout()
    fig.savefig(args.output_dir / "mbias_cohort_weighted.png", dpi=180)
    fig.savefig(args.output_dir / "mbias_cohort_weighted.pdf")
    print(f"reports={len(reports)}")
    print(tsv)


if __name__ == "__main__":
    main()
