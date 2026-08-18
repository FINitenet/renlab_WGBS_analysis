#!/usr/bin/env python3
"""From an input FASTQ folder to paired, sample-level WGBS QC reports.

The workflow discovers R1/R2 files, joins lanes in a stable order, validates
mate names, runs FastQC and MultiQC, and writes a compact TSV plus PNG/PDF
overview. Concatenating gzip members is standards-compliant and does not
decompress/recompress already gzipped inputs.
"""

import argparse
import csv
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path


FASTQ_RE = re.compile(
    r"^(?P<stem>.+?)(?:[._])(?:R)?(?P<read>[12])(?:_001)?"
    r"\.(?:fastq|fq)(?P<gz>\.gz)?$",
    re.IGNORECASE,
)
LANE_RE = re.compile(r"([._-])L\d{1,3}(?=[._-]|$)", re.IGNORECASE)
SAFE_SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Discover paired FASTQ lanes and generate FastQC/MultiQC figures."
    )
    parser.add_argument("input_dir", type=Path, help="Folder containing FASTQ files")
    parser.add_argument(
        "-o", "--outdir", type=Path,
        help="QC output root (default: <input parent>/wgbs_qc)",
    )
    parser.add_argument(
        "--tool-bin",
        type=Path,
        help="Optional directory containing FastQC and MultiQC (otherwise use PATH)",
    )
    parser.add_argument("--threads", type=int, default=12, help="FastQC threads (default: 12)")
    parser.add_argument(
        "--pair-check-reads",
        type=int,
        default=10000,
        help="Number of leading read pairs whose names are checked (default: 10000)",
    )
    parser.add_argument(
        "--no-recursive", action="store_true", help="Do not search input subdirectories"
    )
    parser.add_argument("--force", action="store_true", help="Rebuild staged files and QC")
    parser.add_argument("--dry-run", action="store_true", help="Discover and print the plan only")
    return parser.parse_args()


def log(message):
    print(f"[WGBS-QC] {message}", flush=True)


def die(message):
    raise SystemExit(f"[ERROR] {message}")


def strip_fastq_suffix(name):
    return re.sub(r"\.(?:fastq|fq)(?:\.gz)?$", "", name, flags=re.IGNORECASE)


def normalized_sample(stem):
    sample = LANE_RE.sub("", stem)
    sample = re.sub(r"[._-]{2,}", "_", sample).strip("._-")
    if not sample or not SAFE_SAMPLE_RE.fullmatch(sample):
        die(f"Unsafe or empty sample name inferred from '{stem}': '{sample}'")
    return sample


def discover_fastqs(input_dir, recursive=True):
    if not input_dir.is_dir():
        die(f"Input directory does not exist: {input_dir}")
    iterator = input_dir.rglob("*") if recursive else input_dir.glob("*")
    grouped = defaultdict(lambda: defaultdict(dict))
    unrecognized = []
    for path in sorted(iterator):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if not lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            continue
        match = FASTQ_RE.match(path.name)
        if not match:
            unrecognized.append(path)
            continue
        stem = match.group("stem")
        read = f"R{match.group('read')}"
        sample = normalized_sample(stem)
        if read in grouped[sample][stem]:
            die(f"Duplicate {read} for sample/lane {sample}/{stem}")
        grouped[sample][stem][read] = path.resolve()

    if unrecognized:
        names = "\n  ".join(str(p) for p in unrecognized[:10])
        die(f"FASTQ-like files with unsupported names:\n  {names}")
    if not grouped:
        die(f"No paired FASTQ files found under {input_dir}")

    samples = {}
    for sample, lanes in sorted(grouped.items()):
        r1_files, r2_files = [], []
        for lane, reads in sorted(lanes.items()):
            if set(reads) != {"R1", "R2"}:
                die(f"Unpaired lane {sample}/{lane}: found {sorted(reads)}")
            r1_files.append(reads["R1"])
            r2_files.append(reads["R2"])
        samples[sample] = {"R1": r1_files, "R2": r2_files}
    return samples


def file_signature(path):
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def gzip_file(path):
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def stage_one_read(sample, read, sources, staged_dir, state_dir, force=False):
    output = staged_dir / f"{sample}.{read}.fastq.gz"
    marker = state_dir / f"{sample}.{read}.merge.json"
    expected = {"inputs": [file_signature(p) for p in sources]}
    if not force and output.exists() and output.stat().st_size > 0 and marker.exists():
        try:
            if json.loads(marker.read_text()) == expected:
                log(f"stage skip: {sample} {read}")
                return output
        except (OSError, json.JSONDecodeError):
            pass

    staged_dir.mkdir(parents=True, exist_ok=True)
    state_dir.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{sample}.{read}.", suffix=".tmp", dir=staged_dir)
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        if len(sources) == 1 and gzip_file(sources[0]):
            # Avoid duplicating large single-lane inputs. The state signature
            # still detects a replaced or modified source on the next run.
            temp_path.unlink()
            os.symlink(sources[0], temp_path)
        elif all(gzip_file(p) for p in sources):
            with temp_path.open("wb") as target:
                for source in sources:
                    with source.open("rb") as source_handle:
                        shutil.copyfileobj(source_handle, target, length=16 * 1024 * 1024)
        else:
            with gzip.open(temp_path, "wb", compresslevel=6) as target:
                for source in sources:
                    opener = gzip.open if gzip_file(source) else open
                    with opener(source, "rb") as source_handle:
                        shutil.copyfileobj(source_handle, target, length=16 * 1024 * 1024)
        if temp_path.stat().st_size == 0:
            die(f"Staged output is empty: {output}")
        os.replace(temp_path, output)
        marker.write_text(json.dumps(expected, indent=2) + "\n")
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return output


def mate_key(header):
    token = header[1:].split(None, 1)[0]
    return re.sub(r"/[12]$", "", token)


def read_fastq_record(handle, path):
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if any(not line for line in lines[1:]):
        die(f"Truncated FASTQ record in {path}")
    header, sequence, plus, quality = [line.rstrip("\r\n") for line in lines]
    if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
        die(f"Invalid FASTQ record in {path}: {header[:80]}")
    return header


def check_pair_prefix(r1, r2, max_reads):
    checked = 0
    with gzip.open(r1, "rt", encoding="latin-1") as h1, gzip.open(
        r2, "rt", encoding="latin-1"
    ) as h2:
        while checked < max_reads:
            header1 = read_fastq_record(h1, r1)
            header2 = read_fastq_record(h2, r2)
            if header1 is None or header2 is None:
                if header1 != header2:
                    die(f"R1/R2 have unequal record counts near pair {checked + 1}: {r1}, {r2}")
                break
            if mate_key(header1) != mate_key(header2):
                die(f"Mate-name mismatch at pair {checked + 1}: {header1} != {header2}")
            checked += 1
    return checked


def write_manifest(samples, staged, manifest):
    manifest.parent.mkdir(parents=True, exist_ok=True)
    temp = manifest.with_suffix(manifest.suffix + ".tmp")
    with temp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "R1", "R2"])
        for sample in samples:
            writer.writerow(
                [sample, staged[sample]["R1"], staged[sample]["R2"]]
            )
    os.replace(temp, manifest)


def command_path(name, tool_bin):
    if tool_bin is not None:
        candidate = tool_bin / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    found = shutil.which(name)
    if found:
        return found
    checked = f"{tool_bin} and PATH" if tool_bin is not None else "PATH"
    die(f"Required command not found: {name} (checked {checked})")


def run_command(command, log_path, env):
    log("run: " + " ".join(shlex.quote(str(x)) for x in command))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=env)
    if completed.returncode:
        die(f"Command failed ({completed.returncode}); see {log_path}")


def fastqc_zip_for(fastq, fastqc_dir):
    return fastqc_dir / f"{strip_fastq_suffix(fastq.name)}_fastqc.zip"


def run_fastqc(staged, fastqc_dir, logs_dir, tool_bin, threads, force, env):
    fastqc = command_path("fastqc", tool_bin)
    fastqc_dir.mkdir(parents=True, exist_ok=True)
    fastqs = [staged[s][r] for s in staged for r in ("R1", "R2")]
    pending = [p for p in fastqs if force or not fastqc_zip_for(p, fastqc_dir).is_file()]
    if pending:
        command = [fastqc, "--threads", str(max(1, threads)), "--outdir", str(fastqc_dir)]
        command.extend(str(p) for p in pending)
        run_command(command, logs_dir / "fastqc.log", env)
    else:
        log("FastQC skip: all reports are complete")
    missing = [p for p in fastqs if not fastqc_zip_for(p, fastqc_dir).is_file()]
    if missing:
        die("Missing FastQC reports for: " + ", ".join(str(p) for p in missing))


def parse_fastqc_zip(zip_path):
    with zipfile.ZipFile(zip_path) as archive:
        members = [name for name in archive.namelist() if name.endswith("/fastqc_data.txt")]
        if len(members) != 1:
            die(f"Expected one fastqc_data.txt in {zip_path}, found {len(members)}")
        text = archive.read(members[0]).decode("utf-8", errors="replace")

    basic, statuses = {}, {}
    base_quality = []
    deduplicated_percent = None
    current = None
    for line in text.splitlines():
        if line.startswith(">>") and not line.startswith(">>END_MODULE"):
            fields = line[2:].split("\t")
            current = fields[0]
            statuses[current] = fields[1] if len(fields) > 1 else "unknown"
            continue
        if line.startswith(">>END_MODULE"):
            current = None
            continue
        if current == "Basic Statistics" and line and not line.startswith("#"):
            fields = line.split("\t", 1)
            if len(fields) == 2:
                basic[fields[0]] = fields[1]
        elif current == "Per base sequence quality" and line and not line.startswith("#"):
            fields = line.split("\t")
            try:
                position = fields[0]
                span = 1
                if "-" in position:
                    start, end = map(int, position.split("-", 1))
                    span = end - start + 1
                base_quality.append((float(fields[1]), span))
            except (ValueError, IndexError):
                pass
        elif current == "Sequence Duplication Levels" and line.startswith("#Total Deduplicated Percentage"):
            try:
                deduplicated_percent = float(line.split("\t")[-1])
            except ValueError:
                pass

    weighted_q = None
    min_q = None
    if base_quality:
        weighted_q = sum(q * span for q, span in base_quality) / sum(span for _, span in base_quality)
        min_q = min(q for q, _ in base_quality)
    return {
        "filename": basic.get("Filename", ""),
        "total_sequences": int(basic.get("Total Sequences", "0").replace(",", "")),
        "sequence_length": basic.get("Sequence length", ""),
        "gc_percent": float(basic.get("%GC", "nan")),
        "mean_base_quality": weighted_q,
        "min_base_quality": min_q,
        "deduplicated_percent": deduplicated_percent,
        "pass_modules": sum(value == "pass" for value in statuses.values()),
        "warn_modules": sum(value == "warn" for value in statuses.values()),
        "fail_modules": sum(value == "fail" for value in statuses.values()),
    }


def build_summary(samples, staged, fastqc_dir, summary_path):
    rows = []
    for sample in samples:
        pair_rows = []
        for read in ("R1", "R2"):
            metrics = parse_fastqc_zip(fastqc_zip_for(staged[sample][read], fastqc_dir))
            metrics.update({"sample": sample, "read": read, "lanes": len(samples[sample][read])})
            rows.append(metrics)
            pair_rows.append(metrics)
        if pair_rows[0]["total_sequences"] != pair_rows[1]["total_sequences"]:
            die(
                f"FastQC found unequal R1/R2 counts for {sample}: "
                f"{pair_rows[0]['total_sequences']} vs {pair_rows[1]['total_sequences']}"
            )

    fields = [
        "sample", "read", "lanes", "total_sequences", "sequence_length", "gc_percent",
        "mean_base_quality", "min_base_quality", "deduplicated_percent",
        "pass_modules", "warn_modules", "fail_modules", "filename",
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def make_plots(rows, output_prefix):
    os.environ.setdefault("MPLCONFIGDIR", str(output_prefix.parent / ".mplconfig"))
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        log(f"WARNING: matplotlib unavailable; static plots skipped: {exc}")
        return

    by_sample = defaultdict(dict)
    for row in rows:
        by_sample[row["sample"]][row["read"]] = row
    names = list(by_sample)
    x = list(range(len(names)))
    width = 0.38
    fig_width = max(11, len(names) * 0.75)
    fig, axes = plt.subplots(2, 2, figsize=(fig_width, 9), constrained_layout=True)

    axes[0, 0].bar(x, [by_sample[s]["R1"]["total_sequences"] / 1e6 for s in names], color="#4472C4")
    axes[0, 0].set_title("Read pairs per sample")
    axes[0, 0].set_ylabel("Million pairs")

    for offset, read, color in [(-width / 2, "R1", "#4472C4"), (width / 2, "R2", "#ED7D31")]:
        axes[0, 1].bar(
            [value + offset for value in x], [by_sample[s][read]["gc_percent"] for s in names],
            width=width, label=read, color=color,
        )
        axes[1, 0].bar(
            [value + offset for value in x],
            [by_sample[s][read]["mean_base_quality"] or 0 for s in names],
            width=width, label=read, color=color,
        )
    axes[0, 1].set_title("GC content")
    axes[0, 1].set_ylabel("GC (%)")
    axes[0, 1].legend(frameon=False)
    axes[1, 0].axhline(30, color="#C00000", linestyle="--", linewidth=1, label="Q30")
    axes[1, 0].set_title("Mean per-base quality")
    axes[1, 0].set_ylabel("Phred score")
    axes[1, 0].legend(frameon=False)

    warn = [sum(by_sample[s][r]["warn_modules"] for r in ("R1", "R2")) for s in names]
    fail = [sum(by_sample[s][r]["fail_modules"] for r in ("R1", "R2")) for s in names]
    axes[1, 1].bar(x, warn, color="#FFC000", label="WARN")
    axes[1, 1].bar(x, fail, bottom=warn, color="#C00000", label="FAIL")
    axes[1, 1].set_title("FastQC module flags (R1 + R2)")
    axes[1, 1].set_ylabel("Number of modules")
    axes[1, 1].legend(frameon=False)

    for axis in axes.flat:
        axis.set_xticks(x)
        axis.set_xticklabels(names, rotation=55, ha="right", fontsize=8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.2)
    fig.suptitle("WGBS raw-read QC overview", fontsize=16)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=200)
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)


def run_multiqc(fastqc_dir, multiqc_dir, logs_dir, tool_bin, env):
    multiqc = command_path("multiqc", tool_bin)
    multiqc_dir.mkdir(parents=True, exist_ok=True)
    command = [
        multiqc, "--force", "--module", "fastqc", "--data-format", "json",
        "--filename", "wgbs_qc_multiqc.html", "--outdir", str(multiqc_dir), str(fastqc_dir),
    ]
    run_command(command, logs_dir / "multiqc.log", env)
    report = multiqc_dir / "wgbs_qc_multiqc.html"
    if not report.is_file():
        die(f"MultiQC report was not created: {report}")
    return report


def main():
    args = parse_args()
    if args.threads < 1 or args.pair_check_reads < 0:
        die("--threads must be >=1 and --pair-check-reads must be >=0")
    samples = discover_fastqs(args.input_dir.resolve(), recursive=not args.no_recursive)
    lane_count = sum(len(value["R1"]) for value in samples.values())
    log(f"discovered {len(samples)} paired samples across {lane_count} lanes")
    for sample, reads in samples.items():
        log(f"plan: {sample}: {len(reads['R1'])} lane(s)")
    if args.dry_run:
        return 0

    outdir = (args.outdir or (args.input_dir.resolve().parent / "wgbs_qc")).resolve()
    staged_dir = outdir / "01_merged_fastq"
    qc_dir = outdir / "02_qc"
    fastqc_dir = qc_dir / "fastqc"
    multiqc_dir = qc_dir / "multiqc"
    logs_dir = outdir / "logs"
    state_dir = outdir / ".state"
    for directory in (staged_dir, qc_dir, logs_dir, state_dir):
        directory.mkdir(parents=True, exist_ok=True)

    staged = {}
    for sample, reads in samples.items():
        staged[sample] = {}
        for read in ("R1", "R2"):
            staged[sample][read] = stage_one_read(
                sample, read, reads[read], staged_dir, state_dir, force=args.force
            )
        checked = check_pair_prefix(
            staged[sample]["R1"], staged[sample]["R2"], args.pair_check_reads
        )
        log(f"pair-name check: {sample}: {checked} pairs")

    manifest = outdir / "samples.tsv"
    write_manifest(samples, staged, manifest)
    env = os.environ.copy()
    if args.tool_bin is not None:
        env["PATH"] = str(args.tool_bin) + os.pathsep + env.get("PATH", "")
    env.setdefault("MPLCONFIGDIR", str(outdir / ".mplconfig"))
    run_fastqc(staged, fastqc_dir, logs_dir, args.tool_bin, args.threads, args.force, env)
    rows = build_summary(samples, staged, fastqc_dir, qc_dir / "qc_summary.tsv")
    make_plots(rows, qc_dir / "fastqc_overview")
    report = run_multiqc(fastqc_dir, multiqc_dir, logs_dir, args.tool_bin, env)
    log(f"complete: {report}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        die("Interrupted")
