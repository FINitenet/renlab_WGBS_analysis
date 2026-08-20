#!/usr/bin/env python3
"""Restartable end-to-end WGBS analysis and publication figure workflow.

This is the single user-facing entry point.  It reuses the small, tested stage
scripts in this repository while owning configuration, stage order, provenance
checks and downstream figure generation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_GFF = Path(
    "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/"
    "TAIR10_GFF3_genes_transposons.gff"
)
STAGES = ("qc", "validate", "align", "dedup", "mbias", "extract", "summary", "figures")


def log(message):
    print(f"[{datetime.now():%F %T}] [WGBS] {message}", flush=True)


def die(message):
    raise SystemExit(f"[ERROR] {message}")


def existing(path, label):
    path = Path(path).expanduser().resolve()
    if not path.exists():
        die(f"{label} does not exist: {path}")
    return path


def find_input(analysis_dir):
    for name in ("1_rawdata_merged", "1_rawdata"):
        candidate = analysis_dir / name
        if candidate.is_dir():
            return candidate
    die(f"No 1_rawdata_merged or 1_rawdata directory below {analysis_dir}; use --input-dir")


def find_genome_dir(analysis_dir):
    name = "TAIR10_plus_transgene_bismark_bt2"
    for root in (analysis_dir, *analysis_dir.parents):
        candidate = root / "reference" / name
        if candidate.is_dir():
            return candidate
    die("Combined Bismark index was not found in this project tree; use --genome-dir")


def fasta_in(genome_dir):
    fastas = sorted(genome_dir.glob("*.fa")) + sorted(genome_dir.glob("*.fasta"))
    if len(fastas) != 1:
        die(f"Expected exactly one FASTA in {genome_dir}, found {len(fastas)}")
    return fastas[0]


def parse_clips(value):
    try:
        clips = tuple(int(item) for item in value.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError("M-bias clips must be four comma-separated integers")
    if len(clips) != 4 or min(clips) < 0:
        raise argparse.ArgumentTypeError("M-bias clips must be R1_5,R1_3,R2_5,R2_3")
    return clips


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis_dir", type=Path,
                        help="analysis directory containing 1_rawdata[_merged]")
    parser.add_argument("--input-dir", type=Path, help="override FASTQ input directory")
    parser.add_argument("--outdir", type=Path,
                        help="result root (default: <analysis_dir>/wgbs_joint_results)")
    parser.add_argument("--qc-dir", type=Path,
                        help="raw QC root (default: <analysis_dir>/wgbs_qc)")
    parser.add_argument("--genome-dir", type=Path,
                        help="Bismark genome directory (auto-detected from project/reference)")
    parser.add_argument("--genome-fa", type=Path, help="combined FASTA (auto-detected)")
    parser.add_argument("--gff", type=Path,
                        help="gene/TE GFF (default: project standard TAIR10 GFF when available)")
    parser.add_argument("--feature-bed", type=Path,
                        help="locus features (default: <genome-dir>/transgene_features.bed)")
    parser.add_argument("--locus", help="locus figure region (default: entire feature contig)")
    parser.add_argument("--library-type", choices=("directional", "non_directional", "pbat"),
                        default="directional")
    parser.add_argument("--mbias-clips", type=parse_clips, default=(10, 0, 10, 3),
                        metavar="R1_5,R1_3,R2_5,R2_3",
                        help="reviewed clipping values (default: 10,0,10,3 for this cohort)")
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--tool-bin", type=Path, help="environment bin directory")
    parser.add_argument("--from-stage", choices=STAGES, default="qc")
    parser.add_argument("--to-stage", choices=STAGES, default="figures")
    parser.add_argument("--dry-run", action="store_true", help="validate and print commands only")
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fai(path):
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                rows.append((fields[0], int(fields[1])))
    if not rows:
        die(f"Empty FASTA index: {path}")
    return rows


def infer_locus(feature_bed, fai):
    if not feature_bed or not feature_bed.is_file():
        return None
    contigs = dict(read_fai(fai))
    names = []
    with feature_bed.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                names.append(line.split("\t", 1)[0])
    if not names or len(set(names)) != 1:
        die(f"Feature BED must contain one contig for automatic locus inference: {feature_bed}")
    chrom = names[0]
    if chrom not in contigs:
        die(f"Feature contig {chrom} is absent from {fai}")
    return f"{chrom}:1-{contigs[chrom]}"


def command_text(command):
    return " ".join(shlex.quote(str(item)) for item in command)


def run(command, env, dry_run=False):
    log(command_text(command))
    if not dry_run:
        subprocess.run([str(item) for item in command], check=True, env=env)


def write_state(path, state, dry_run):
    if dry_run:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
    os.replace(temp, path)


def validate_state(outdir, state, dry_run):
    path = outdir / ".oneclick_config.json"
    if path.exists():
        previous = json.loads(path.read_text())
        keys = ("genome_dir", "genome_sha256", "samples_tsv")
        changed = [key for key in keys if previous.get(key) != state.get(key)]
        if changed:
            die(f"Refusing to reuse {outdir}: analysis identity changed ({', '.join(changed)})")
    write_state(path, state, dry_run)


def main():
    args = parse_args()
    if args.threads < 1:
        die("--threads must be at least 1")
    start, stop = STAGES.index(args.from_stage), STAGES.index(args.to_stage)
    if start > stop:
        die("--from-stage must not come after --to-stage")

    analysis_dir = existing(args.analysis_dir, "analysis directory")
    input_dir = existing(args.input_dir, "FASTQ input") if args.input_dir else find_input(analysis_dir)
    genome_dir = existing(args.genome_dir, "Bismark genome") if args.genome_dir else find_genome_dir(analysis_dir)
    genome_fa = existing(args.genome_fa, "genome FASTA") if args.genome_fa else fasta_in(genome_dir)
    fai = genome_fa.with_suffix(genome_fa.suffix + ".fai")
    if not fai.is_file():
        die(f"Missing FASTA index: {fai}; run samtools faidx first")
    gff = args.gff.expanduser().resolve() if args.gff else DEFAULT_GFF
    if stop >= STAGES.index("figures") and not gff.is_file():
        die(f"Gene/TE GFF does not exist: {gff}; use --gff")
    feature_bed = args.feature_bed.expanduser().resolve() if args.feature_bed else genome_dir / "transgene_features.bed"
    if not feature_bed.is_file():
        feature_bed = None
    locus = args.locus or infer_locus(feature_bed, fai)

    outdir = (args.outdir or analysis_dir / "wgbs_joint_results").expanduser().resolve()
    qc_dir = (args.qc_dir or analysis_dir / "wgbs_qc").expanduser().resolve()
    samples_tsv = qc_dir / "samples.tsv"
    figures = outdir / "08_reference_figures"
    tool_bin = args.tool_bin.expanduser().resolve() if args.tool_bin else None
    env = os.environ.copy()
    if tool_bin:
        env["PATH"] = str(tool_bin) + os.pathsep + env.get("PATH", "")
    env.update({
        "PROJECT_DIR": str(analysis_dir), "SAMPLES_TSV": str(samples_tsv),
        "OUTDIR": str(outdir), "GENOME_FA": str(genome_fa),
        "RAW_QC_DIR": str(qc_dir),
        "BISMARK_GENOME_DIR": str(genome_dir), "LIBRARY_TYPE": args.library_type,
        "THREADS": str(args.threads), "FASTQC_THREADS": str(min(args.threads, 12)),
        "BISMARK_PARALLEL": str(max(1, min(args.threads // 4, 4))),
        "SAMTOOLS_THREADS": str(min(args.threads, 8)),
        "EXTRACT_PARALLEL": str(max(1, min(args.threads // 6, 4))),
    })
    if tool_bin:
        env["TOOL_BIN"] = str(tool_bin)

    selected = set(STAGES[start:stop + 1])
    if "qc" in selected:
        command = [sys.executable, HERE / "wgbs_qc.py", input_dir, "-o", qc_dir,
                   "--threads", min(args.threads, 12)]
        run(command, env, args.dry_run)
    if not args.dry_run and not samples_tsv.is_file():
        die(f"Sample sheet not found: {samples_tsv}; include qc stage or set --qc-dir")

    state = {
        "analysis_dir": str(analysis_dir), "input_dir": str(input_dir),
        "outdir": str(outdir), "samples_tsv": str(samples_tsv),
        "genome_dir": str(genome_dir), "genome_fasta": str(genome_fa),
        "genome_sha256": sha256(genome_fa), "gff": str(gff),
        "library_type": args.library_type, "mbias_clips": list(args.mbias_clips),
    }
    validate_state(outdir, state, args.dry_run)

    stage_commands = {
        "validate": [HERE / "00_validate.sh"],
        "align": [HERE / "02_trim_decontam_map.sh"],
        "dedup": [HERE / "03_dedup.sh"],
        "mbias": [HERE / "04_meth_extractor.sh", "mbias"],
        "extract": [HERE / "04_meth_extractor.sh", "extract"],
    }
    for stage in ("validate", "align", "dedup", "mbias"):
        if stage in selected:
            run(stage_commands[stage], env, args.dry_run)
            if stage == "mbias":
                run([sys.executable, HERE / "summarize_mbias.py", outdir / "05_methylation",
                     "-o", outdir / "05_methylation" / "mbias_summary"], env, args.dry_run)

    if "extract" in selected:
        r1_5, r1_3, r2_5, r2_3 = args.mbias_clips
        env.update({"MBIAS_REVIEWED": "true", "IGNORE_R1_5": str(r1_5),
                    "IGNORE_R1_3": str(r1_3), "IGNORE_R2_5": str(r2_5),
                    "IGNORE_R2_3": str(r2_3)})
        run(stage_commands["extract"], env, args.dry_run)
        run([HERE / "05_multiqc.sh"], env, args.dry_run)

    if "summary" in selected:
        genome_size = sum(length for _, length in read_fai(fai))
        run([sys.executable, HERE / "summarize_results.py", outdir, samples_tsv,
             qc_dir / "02_qc" / "qc_summary.tsv", "-o", outdir / "07_summary",
             "--genome-size", genome_size,
             "--effective-r1", 150 - args.mbias_clips[0] - args.mbias_clips[1],
             "--effective-r2", 150 - args.mbias_clips[2] - args.mbias_clips[3]], env, args.dry_run)

    if "figures" in selected:
        figure_cmd = [sys.executable, HERE / "wgbs_reference_figures.py", "prepare",
                      "--methylation-dir", outdir / "05_methylation", "--gff", gff,
                      "--chrom-sizes", fai, "--outdir", figures,
                      "--workers", 2, "--threads", min(args.threads, 8)]
        run(figure_cmd, env, args.dry_run)
        if locus:
            locus_cmd = [sys.executable, HERE / "wgbs_reference_figures.py", "locus",
                         "--outdir", figures, "--region", locus]
            if feature_bed:
                locus_cmd.extend(["--feature-bed", feature_bed])
            run(locus_cmd, env, args.dry_run)

    state["completed_through"] = STAGES[stop]
    state["completed_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
    write_state(outdir / ".oneclick_config.json", state, args.dry_run)
    log(f"complete through {STAGES[stop]}: {outdir}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        die(f"stage failed with exit code {error.returncode}: {command_text(error.cmd)}")
    except KeyboardInterrupt:
        die("interrupted")
