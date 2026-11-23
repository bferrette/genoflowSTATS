#!/usr/bin/env python3
"""genoflowSTATS command-line interface.

Reference-free genome assembly metrics, direct samtools stats parsing, robust
coverage estimates, BAM QC, and multi-assembly comparison.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import shlex
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Iterator, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter

try:
    import pysam
except Exception:  # optional for commands that do not read BAM/CRAM
    pysam = None

try:
    from rich.console import Console
    console = Console(stderr=True)
except Exception:
    class _Console:
        def log(self, *args: Any, **kwargs: Any) -> None:
            print(*args, file=sys.stderr)
    console = _Console()

from . import __version__

HEAT = LinearSegmentedColormap.from_list(
    "genoflow_heat", ["#26002f", "#6a176e", "#b73779", "#ed5f5f", "#fa9f42", "#ffe66d"]
)
HEAT_R = LinearSegmentedColormap.from_list(
    "genoflow_heat_r", ["#ffe66d", "#fa9f42", "#ed5f5f", "#b73779", "#6a176e", "#26002f"]
)


PALETTES = {
    "heat": ("#6a176e", "#ed5f5f", "#fa9f42", "magma"),
    "ocean": ("#175676", "#2a9d8f", "#e9c46a", "viridis"),
    "colorblind": ("#0072B2", "#D55E00", "#009E73", "cividis"),
}
COLORS = PALETTES["colorblind"]


def open_text(path: str | Path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "rt", encoding="utf-8")


def finite(v: Any) -> bool:
    try:
        return math.isfinite(float(v))
    except Exception:
        return False


def parse_number(text: str) -> int | float | str:
    s = text.strip()
    try:
        if re.fullmatch(r"[-+]?\d+", s):
            return int(s)
        return float(s)
    except Exception:
        return s


def write_kv_tsv(path: Path, values: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        for k, v in values.items():
            if v is None or (isinstance(v, float) and not math.isfinite(v)):
                v = "NA"
            w.writerow([k, v])


def write_json(path: Path, obj: Any) -> None:
    def clean(x: Any) -> Any:
        if isinstance(x, dict):
            return {k: clean(v) for k, v in x.items()}
        if isinstance(x, list):
            return [clean(v) for v in x]
        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, (np.floating,)):
            x = float(x)
        if isinstance(x, float) and not math.isfinite(x):
            return None
        return x
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(obj), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def style_ax(ax: plt.Axes) -> None:
    ax.set_facecolor("#f7f9fb")
    ax.figure.set_facecolor("white")
    ax.grid(True, color="#dce3e8", linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)


def save_svg(fig: plt.Figure, path: Path, dpi: int = 600) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, format="svg", dpi=max(600, int(dpi)), bbox_inches="tight")
    plt.close(fig)


# ------------------------------- FASTA metrics -------------------------------

def fasta_records(path: str | Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line_number, raw in enumerate(fh, 1):
            if raw.startswith(">"):
                if name is not None:
                    if not chunks:
                        raise ValueError(f"Empty FASTA sequence for {name}")
                    yield name, "".join(chunks)
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                name = fields[0]
                chunks = []
            else:
                sequence = raw.strip()
                if not sequence:
                    continue
                if name is None:
                    raise ValueError(f"Sequence before FASTA header at line {line_number}")
                if any(c.isspace() for c in sequence):
                    raise ValueError(f"Whitespace within FASTA sequence at line {line_number}")
                chunks.append(sequence)
        if name is not None:
            if not chunks:
                raise ValueError(f"Empty FASTA sequence for {name}")
            yield name, "".join(chunks)


def nxx(lengths: Sequence[int], denominator: float, pct: float) -> tuple[int, int]:
    if denominator <= 0 or not lengths:
        return 0, 0
    target = denominator * pct / 100.0
    total = 0
    for i, length in enumerate(sorted(lengths, reverse=True), start=1):
        total += length
        if total >= target:
            return int(length), i
    return 0, 0


def assembly_metrics(path: str | Path, genome_size: int | None = None, gap_min: int = 10) -> tuple[dict[str, Any], list[dict[str, Any]], np.ndarray, np.ndarray]:
    lengths: list[int] = []
    contig_lengths: list[int] = []
    scaffold_rows: list[dict[str, Any]] = []
    gap_lengths: list[int] = []
    total_gc = total_atgc = total_n = total_ambig = 0
    gap_re = re.compile(rf"N{{{int(gap_min)},}}", re.I)

    for name, seq in fasta_records(path):
        s = seq.upper()
        L = len(s)
        lengths.append(L)
        n = s.count("N")
        gc = s.count("G") + s.count("C")
        atgc = s.count("A") + s.count("T") + gc
        ambig = max(0, L - atgc - n)
        total_gc += gc; total_atgc += atgc; total_n += n; total_ambig += ambig
        gaps = list(gap_re.finditer(s))
        gl = [m.end() - m.start() for m in gaps]
        gap_lengths.extend(gl)
        pieces = [x for x in gap_re.split(s) if x]
        contig_lengths.extend(len(x) for x in pieces)
        scaffold_rows.append({
            "scaffold": name,
            "length": L,
            "gc_percent": (100.0 * gc / atgc if atgc else math.nan),
            "n_bases": n,
            "gap_count": len(gl),
            "gap_bp": sum(gl),
        })

    if not lengths:
        raise ValueError(f"No FASTA records found: {path}")
    lens = np.asarray(lengths, dtype=np.int64)
    contigs = np.asarray(contig_lengths, dtype=np.int64)
    gaps = np.asarray(gap_lengths, dtype=np.int64)
    total = int(lens.sum())
    g = int(genome_size) if genome_size and int(genome_size) > 0 else None
    n50, l50 = nxx(lengths, total, 50); n90, l90 = nxx(lengths, total, 90)
    ng50, lg50 = nxx(lengths, g, 50) if g else (0, 0)
    ng90, lg90 = nxx(lengths, g, 90) if g else (0, 0)
    cn50, cl50 = nxx(contig_lengths, int(contigs.sum()), 50) if contigs.size else (0, 0)
    cn90, cl90 = nxx(contig_lengths, int(contigs.sum()), 90) if contigs.size else (0, 0)
    sdesc = np.sort(lens)[::-1]
    sumsq = float(np.dot(lens.astype(float), lens.astype(float)))
    gap_bp = int(gaps.sum()) if gaps.size else 0

    def topfrac(n: int) -> float:
        return float(sdesc[:n].sum() / total) if total else math.nan

    vals: dict[str, Any] = {
        "assembly": Path(path).name,
        "scaffold_count": int(lens.size),
        "total_length": total,
        "canonical_bases": total_atgc,
        "n_bases": total_n,
        "ambiguous_non_N_bases": total_ambig,
        "gc_percent": 100.0 * total_gc / total_atgc if total_atgc else math.nan,
        "largest_scaffold": int(lens.max()),
        "smallest_scaffold": int(lens.min()),
        "mean_scaffold_length": float(lens.mean()),
        "median_scaffold_length": float(np.median(lens)),
        "scaffold_length_mad": float(np.median(np.abs(lens - np.median(lens)))),
        "n50": n50, "l50": l50, "n90": n90, "l90": l90,
        "ng50": ng50 or None, "lg50": lg50 or None,
        "ng90": ng90 or None, "lg90": lg90 or None,
        "auN": sumsq / total if total else math.nan,
        "auNG": sumsq / g if g else None,
        "genome_size_estimate": g,
        "assembly_to_genome_ratio": total / g if g else None,
        "assembly_vs_genome_percent": 100.0 * (total - g) / g if g else None,
        "contig_count_gap_split": int(contigs.size),
        "contig_n50": cn50, "contig_l50": cl50, "contig_n90": cn90, "contig_l90": cl90,
        "gap_min_run": int(gap_min), "gap_count": int(gaps.size), "gap_bp": gap_bp,
        "gap_fraction": gap_bp / total if total else math.nan,
        "gaps_per_100mb": gaps.size * 1e8 / total if total else math.nan,
        "ns_per_100kbp": total_n * 1e5 / total if total else math.nan,
        "median_gap": float(np.median(gaps)) if gaps.size else 0.0,
        "mean_gap": float(gaps.mean()) if gaps.size else 0.0,
        "p95_gap": float(np.percentile(gaps, 95)) if gaps.size else 0.0,
        "largest_gap": int(gaps.max()) if gaps.size else 0,
        "top1_fraction": topfrac(1), "top5_fraction": topfrac(5), "top10_fraction": topfrac(10),
        "top20_fraction": topfrac(20), "top50_fraction": topfrac(50),
    }
    return vals, scaffold_rows, lens, gaps


def nx_profile(lengths: np.ndarray, denominator: int) -> np.ndarray:
    if denominator <= 0 or not lengths.size:
        return np.zeros(100, dtype=float)
    ordered = np.sort(lengths)[::-1]
    indexes = np.searchsorted(np.cumsum(ordered, dtype=np.int64), denominator * np.arange(1, 101) / 100, side="left")
    return np.asarray([ordered[i] if i < len(ordered) else np.nan for i in indexes], dtype=float)


def plot_fasta_metrics(vals: dict[str, Any], lens: np.ndarray, outdir: Path, prefix: str, dpi: int) -> None:
    s = np.sort(lens)[::-1]
    cum = np.cumsum(s) / s.sum()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(np.arange(1, s.size + 1), cum, color=COLORS[0], linewidth=2)
    ax.fill_between(np.arange(1, s.size + 1), cum, color=COLORS[2], alpha=0.28)
    ax.set_xscale("log"); ax.set_ylim(0, 1.01)
    ax.set_xlabel("Scaffold rank (log scale)"); ax.set_ylabel("Cumulative assembly fraction")
    ax.set_title(f"{prefix}: cumulative scaffold length\n{len(s):,} scaffolds · {int(s.sum()):,} bp", loc="left")
    style_ax(ax); save_svg(fig, outdir / f"{prefix}.scaffold_cumulative.svg", dpi)

    x = np.arange(1, 101)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, nx_profile(lens, int(lens.sum())), color=COLORS[0], linewidth=2, label="Nx")
    if vals.get("genome_size_estimate"):
        ax.plot(x, nx_profile(lens, int(vals["genome_size_estimate"])), color=COLORS[2], linewidth=2, label="NGx")
    ax.axvline(50, color="#8996a3", linestyle=":", linewidth=1)
    ax.set_xlabel("Assembly or estimated genome covered (%)"); ax.set_ylabel("Scaffold length (bp)"); ax.set_title(f"{prefix}: Nx/NGx profile", loc="left")
    ax.legend(frameon=False); style_ax(ax); save_svg(fig, outdir / f"{prefix}.nx_ngx.svg", dpi)


# -------------------------- histogram utility metrics -------------------------

def inverse_ecdf_quantile(hist: np.ndarray, q: float) -> float:
    total = int(hist.sum())
    if total <= 0:
        return math.nan
    rank = max(1, int(math.ceil(q * total)))
    cs = np.cumsum(hist, dtype=np.int64)
    return float(np.searchsorted(cs, rank, side="left"))


def exact_mad_hist(hist: np.ndarray, median: int) -> float:
    if hist.sum() <= 0:
        return math.nan
    maxdev = max(median, len(hist) - 1 - median)
    dev = np.zeros(maxdev + 1, dtype=np.int64)
    for d in np.flatnonzero(hist):
        dev[abs(int(d) - median)] += int(hist[d])
    return inverse_ecdf_quantile(dev, 0.5)


def depth_metrics_from_hist(hist: np.ndarray, scale_mad: bool = True, mad_factor: float = 3.0) -> dict[str, Any]:
    n = int(hist.sum())
    if n <= 0:
        raise ValueError("Empty coverage histogram")
    med = int(inverse_ecdf_quantile(hist, 0.5))
    positive = hist.copy()
    if positive.size:
        positive[0] = 0
    npos = int(positive.sum())
    med_pos = int(inverse_ecdf_quantile(positive, 0.5)) if npos else 0
    mode = int(np.argmax(hist))
    mode_pos = int(np.argmax(positive)) if npos else 0
    mad_raw = exact_mad_hist(hist, med)
    mad = mad_raw * 1.4826 if scale_mad else mad_raw
    q = {p: inverse_ecdf_quantile(hist, p / 100.0) for p in (1, 5, 25, 50, 75, 95, 99)}
    depths = np.arange(hist.size, dtype=np.float64)
    weighted = float(np.dot(depths, hist.astype(float)))
    mean = weighted / n
    rev = np.cumsum(hist[::-1], dtype=np.int64)[::-1]
    vals: dict[str, Any] = {
        "coverage_sites_total": n,
        "coverage_zero_sites": int(hist[0]) if hist.size else 0,
        "coverage_zero_fraction": (int(hist[0]) / n if hist.size else 0.0),
        "coverage_covered_sites": npos,
        "coverage_fraction_ge_1x": npos / n,
        "coverage_mean": mean,
        "coverage_median": med,
        "coverage_median_positive": med_pos,
        "coverage_mode": mode,
        "coverage_mode_positive": mode_pos,
        "coverage_mad_raw": mad_raw,
        "coverage_mad_scaled": mad,
        "coverage_iqr": q[75] - q[25],
        "coverage_robust_cv": (mad / med if med > 0 else None),
        "coverage_p01": q[1], "coverage_p05": q[5], "coverage_p25": q[25],
        "coverage_p50": q[50], "coverage_p75": q[75], "coverage_p95": q[95], "coverage_p99": q[99],
        "coverage_mad_min": max(0.0, med - mad_factor * mad),
        "coverage_mad_max": med + mad_factor * mad,
    }
    for t in (1, 5, 10, 20, 30, 50, 100):
        vals[f"coverage_fraction_ge_{t}x"] = float(rev[t] / n) if t < rev.size else 0.0
    baseline = med_pos or med
    for label, mult in (("025", .25), ("050", .5), ("150", 1.5), ("200", 2.0)):
        threshold = int(math.ceil(baseline * mult))
        vals[f"coverage_threshold_{label}_baseline"] = threshold
        vals[f"coverage_fraction_ge_{label}_baseline"] = float(rev[threshold] / n) if threshold < rev.size else 0.0
    return vals


def plot_coverage(hist: np.ndarray, vals: dict[str, Any], outdir: Path, prefix: str, dpi: int) -> None:
    n = int(hist.sum()); x = np.arange(hist.size)
    p99 = vals.get("coverage_p99") or 0
    xmax = max(5, min(hist.size - 1, int(max(p99 * 1.15, (vals.get("coverage_mad_max") or 0) * 1.1))))
    mask = x <= xmax
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x[mask], hist[mask] / n, color=plt.get_cmap(COLORS[3])(np.linspace(.25, .85, int(mask.sum()))), width=1.0, linewidth=0)
    ax.axvline(vals["coverage_median_positive"], color="#26002f", linewidth=1.8, label=f"Median covered={vals['coverage_median_positive']}×")
    ax.set_xlabel("Depth (×)"); ax.set_ylabel("Site fraction"); ax.yaxis.set_major_formatter(PercentFormatter(1))
    scope = "reference/target sites" if vals.get("coverage_histogram_denominator") is not None else "covered sites only; zero depth unknown"
    ax.set_title(f"{prefix}: coverage distribution\n{scope}; median covered {vals['coverage_median_positive']}×", loc="left")
    ax.legend(frameon=False); style_ax(ax); save_svg(fig, outdir / f"{prefix}.coverage_histogram.svg", dpi)

    cdf = np.cumsum(hist, dtype=np.int64) / n
    surv = np.cumsum(hist[::-1], dtype=np.int64)[::-1] / n
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x[:xmax+1], cdf[:xmax+1], color=COLORS[0], linewidth=2, label="CDF")
    ax.plot(x[:xmax+1], surv[:xmax+1], color=COLORS[1], linewidth=2, label="Fraction ≥ depth")
    ax.set_ylim(0, 1.01); ax.set_xlabel("Depth"); ax.set_ylabel("Fraction")
    ax.set_title(f"{prefix}: cumulative coverage\n{scope}", loc="left"); ax.yaxis.set_major_formatter(PercentFormatter(1)); ax.legend(frameon=False)
    style_ax(ax); save_svg(fig, outdir / f"{prefix}.coverage_cumulative.svg", dpi)


# ------------------------------- samtools stats ------------------------------

def _sn_key(key: str) -> str:
    k = key.rstrip(":").strip().lower()
    k = re.sub(r"\s*\(%\)\s*", " percent ", k)
    k = re.sub(r"[^a-z0-9]+", "_", k).strip("_")
    return k


def parse_samtools_stats(path: str | Path, collect_rfs: bool = False) -> dict[str, Any]:
    sn: dict[str, Any] = {}
    cov_rows: list[tuple[str, int, int]] = []
    mapq: Counter[int] = Counter()
    indels: list[tuple[int, int, int]] = []
    gcd: list[dict[str, float]] = []
    rfs_summary: dict[str, Any] | None = None
    rfs_rows: list[dict[str, Any]] = []
    rfs_detail_count = 0
    rl: Counter[int] = Counter()
    unknown_tags: Counter[str] = Counter()
    command_line: str | None = None

    with open_text(path) as fh:
        for raw in fh:
            if raw.startswith("# The command line was:"):
                command_line = raw.split(":", 1)[1].strip()
            if not raw or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if not p:
                continue
            tag = p[0]
            try:
                if tag == "SN" and len(p) >= 3:
                    sn[_sn_key(p[1])] = parse_number(p[2])
                elif tag == "COV" and len(p) >= 4:
                    cov_rows.append((p[1], int(float(p[2])), int(float(p[3]))))
                elif tag == "MAPQ" and len(p) >= 3:
                    mapq[int(float(p[1]))] += int(float(p[2]))
                elif tag == "ID" and len(p) >= 4:
                    indels.append((int(float(p[1])), int(float(p[2])), int(float(p[3]))))
                elif tag == "GCD" and len(p) >= 8:
                    gcd.append({
                        "gc_percent": float(p[1]), "unique_sequence_percent": float(p[2]),
                        "depth_p10": float(p[3]), "depth_p25": float(p[4]), "depth_p50": float(p[5]),
                        "depth_p75": float(p[6]), "depth_p90": float(p[7]),
                    })
                elif tag == "RFS" and len(p) >= 5:
                    # Summary row: total count, output count, average GC, min, max, avg, total length
                    if len(p) >= 8 and re.fullmatch(r"\d+", p[1]) and re.fullmatch(r"\d+", p[2]):
                        rfs_summary = {
                            "reference_total_count": int(p[1]), "reference_output_count": int(p[2]),
                            "reference_average_gc": float(p[3]), "reference_min_length": int(p[4]),
                            "reference_max_length": int(p[5]), "reference_average_length": float(p[6]),
                            "reference_total_length": int(p[7]),
                        }
                    else:
                        rfs_detail_count += 1
                        if collect_rfs:
                            rfs_rows.append({
                                "sequence": p[1], "length": int(float(p[2])), "gc_fraction": float(p[3]),
                                "unknown_bases": int(float(p[4])),
                            })
                elif tag == "RL" and len(p) >= 3:
                    rl[int(float(p[1]))] += int(float(p[2]))
                else:
                    unknown_tags[tag] += 1
            except (ValueError, IndexError):
                continue

    return {
        "sn": sn, "cov_rows": cov_rows, "mapq": mapq, "indels": indels,
        "gcd": gcd, "rfs_summary": rfs_summary, "rfs_rows": rfs_rows,
        "read_length": rl, "unknown_tags": unknown_tags, "rfs_detail_count": rfs_detail_count,
        "command_line": command_line,
    }


def stats_reference_command(command: str | None) -> tuple[bool, str | None]:
    """Recognize unrestricted samtools stats commands without executing them."""
    if not command:
        return False, None
    try:
        args = shlex.split(command)
    except ValueError:
        return False, None
    if args and Path(args[0]).name == "samtools":
        args = args[1:]
    if not args or args[0] != "stats":
        return False, None
    args = args[1:]
    value_options = {"-@", "--threads", "-c", "--coverage", "-r", "--ref-seq",
                     "--GC-depth", "--ref-stats-chunk"}
    switches = {"--ref-stats"}
    reference = None
    positional = []
    i = 0
    while i < len(args):
        arg = args[i]
        if arg in value_options:
            if i + 1 >= len(args):
                return False, None
            if arg in ("-r", "--ref-seq"):
                reference = args[i + 1]
            i += 2
        elif arg.startswith("--ref-seq="):
            reference = arg.split("=", 1)[1]
            i += 1
        elif arg.startswith("--coverage=") or arg.startswith("--threads=") or arg in switches:
            i += 1
        elif arg.startswith("-"):
            return False, None
        else:
            positional.append(arg)
            i += 1
    return len(positional) == 1, reference


def resolve_stats_reference(stats_path: str | Path, reference: str | None) -> Path | None:
    if not reference:
        return None
    candidate = Path(reference)
    if candidate.is_absolute():
        # A stats report read through GVFS may contain paths from the remote host.
        for parent in Path(stats_path).absolute().parents:
            if parent.name.startswith("sftp:"):
                remote_candidate = parent / reference.lstrip("/")
                if remote_candidate.is_file():
                    return remote_candidate
                return None
    elif not candidate.is_file():
        candidate = Path(stats_path).parent / candidate
    return candidate if candidate.is_file() else None


def fasta_total_length(path: str | Path) -> int:
    """Sum sequence-line lengths without retaining whole scaffolds."""
    total = 0
    name: str | None = None
    record_length = 0
    with open_text(path) as fh:
        for line_number, raw in enumerate(fh, 1):
            if raw.startswith(">"):
                if name is not None and record_length == 0:
                    raise ValueError(f"Empty FASTA sequence for {name}")
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                name = fields[0]
                record_length = 0
            else:
                sequence = raw.strip()
                if not sequence:
                    continue
                if name is None:
                    raise ValueError(f"Sequence before FASTA header at line {line_number}")
                if len(sequence.split()) != 1:
                    raise ValueError(f"Whitespace within FASTA sequence at line {line_number}")
                record_length += len(sequence)
                total += len(sequence)
    if name is None or record_length == 0:
        raise ValueError(f"No nonempty FASTA sequence at end of {path}")
    return total


def denominator_from_stats(data: dict[str, Any]) -> int | None:
    sn = data["sn"]
    target = sn.get("bases_inside_the_target")
    if finite(target) and float(target) > 0:
        return int(float(target))
    rfs = data.get("rfs_summary")
    if rfs and int(rfs.get("reference_total_length", 0)) > 0:
        return int(rfs["reference_total_length"])
    return None


def coverage_hist_from_stats(data: dict[str, Any], reference_length: int | None = None) -> tuple[np.ndarray | None, dict[str, Any]]:
    rows = data["cov_rows"]
    if not rows:
        return None, {}
    maxd = max(d for _, d, _ in rows)
    hist = np.zeros(maxd + 1, dtype=np.int64)
    overflow_sites = 0
    for label, d, n in rows:
        match = re.fullmatch(r"\[(\d+)-(\d+)\]", label)
        if match and (int(match.group(1)) != int(match.group(2)) or int(match.group(1)) != d):
            raise ValueError("COV requires unit-width bins; regenerate samtools stats with --coverage 1,1000,1")
        if not match and not re.fullmatch(r"\[\d+<\]", label):
            raise ValueError(f"Unsupported COV bin: {label}")
        hist[d] += n
        if "<" in label:
            overflow_sites += n
    covered = int(hist.sum())
    denominator = denominator_from_stats(data)
    if reference_length is not None:
        if reference_length <= 0:
            raise ValueError("reference length must be positive")
        if denominator is not None and denominator != reference_length:
            raise ValueError("reference length conflicts with the samtools target/reference denominator")
        denominator = reference_length
    if denominator is not None and denominator < covered:
        raise ValueError("coverage sites exceed the supplied reference/target length")
    zero = 0
    if denominator is not None:
        zero = max(0, denominator - covered)
        hist[0] += zero
    meta = {
        "coverage_histogram_covered_sites": covered,
        "coverage_histogram_denominator": denominator,
        "coverage_fraction_scope": "reference_or_target" if denominator is not None else "covered_sites_only",
        "coverage_histogram_reconstructed_zero_sites": zero if denominator is not None else None,
        "coverage_overflow_sites": overflow_sites,
        "coverage_overflow_depth_lower_bound": maxd if overflow_sites else None,
        "coverage_mean_is_lower_bound": bool(overflow_sites),
    }
    return hist, meta


def mapq_metrics(hist_counter: Counter[int]) -> tuple[dict[str, Any], np.ndarray]:
    if not hist_counter:
        return {}, np.zeros(0, dtype=np.int64)
    nmax = max(hist_counter)
    hist = np.zeros(nmax + 1, dtype=np.int64)
    for q, n in hist_counter.items():
        hist[q] = n
    n = int(hist.sum())
    x = np.arange(hist.size, dtype=float)
    rev = np.cumsum(hist[::-1], dtype=np.int64)[::-1]
    vals = {
        "mapq_observations": n,
        "mapq_mean": float(np.dot(x, hist) / n) if n else None,
        "mapq_median": inverse_ecdf_quantile(hist, .5) if n else None,
        "mapq_mode": int(np.argmax(hist)) if n else None,
        "mapq_p05": inverse_ecdf_quantile(hist, .05) if n else None,
        "mapq_p95": inverse_ecdf_quantile(hist, .95) if n else None,
        "mapq0_fraction": hist[0] / n if n else None,
        "mapq_lt10_fraction": hist[:10].sum() / n if n else None,
        "mapq_lt20_fraction": hist[:20].sum() / n if n else None,
        "mapq_ge30_fraction": (rev[30] / n if hist.size > 30 else 0.0) if n else None,
        "mapq_ge50_fraction": (rev[50] / n if hist.size > 50 else 0.0) if n else None,
    }
    return vals, hist


def indel_metrics(rows: list[tuple[int, int, int]], mapped_cigar: float | None) -> dict[str, Any]:
    ins_events = sum(i for _, i, _ in rows); del_events = sum(d for _, _, d in rows)
    ins_bases = sum(k * i for k, i, _ in rows); del_bases = sum(k * d for k, _, d in rows)
    vals: dict[str, Any] = {
        "insertion_events": ins_events, "deletion_events": del_events,
        "insertion_bases": ins_bases, "deletion_bases": del_bases,
        "indel_events_total": ins_events + del_events,
        "indel_bases_total": ins_bases + del_bases,
    }
    if mapped_cigar and mapped_cigar > 0:
        vals["indel_bases_per_aligned_base"] = (ins_bases + del_bases) / mapped_cigar
    return vals


def gc_depth_metrics(rows: list[dict[str, float]]) -> dict[str, Any]:
    if not rows:
        return {}
    med = np.asarray([r["depth_p50"] for r in rows], dtype=float)
    gc = np.asarray([r["gc_percent"] for r in rows], dtype=float)
    cumulative = np.asarray([r["unique_sequence_percent"] for r in rows], dtype=float)
    weights = np.diff(np.r_[0.0, cumulative])
    good = np.isfinite(med) & np.isfinite(gc)
    med = med[good]; gc = gc[good]; weights = weights[good]
    if not med.size:
        return {}
    center = float(np.median(med)); mad = float(np.median(np.abs(med - center)))
    vals = {
        "gc_depth_bins": int(med.size),
        "gc_depth_median_of_medians": center,
        "gc_depth_median_mad": mad,
        "gc_depth_robust_cv": (1.4826 * mad / center if center > 0 else None),
        "gc_depth_p50_min": float(med.min()), "gc_depth_p50_max": float(med.max()),
        "gc_depth_p50_range_ratio": ((np.percentile(med,95)-np.percentile(med,5))/center if center > 0 else None),
    }
    if med.size >= 3:
        coeff = np.polyfit(gc, med, 1)
        vals["gc_depth_linear_slope_per_gc_percent"] = float(coeff[0])
    if np.all(np.isfinite(weights)) and np.all(weights >= 0) and weights.sum() > 0:
        vals["gc_depth_weighted_mean_p50"] = float(np.average(med, weights=weights))
    return vals


def read_length_metrics(hist_counter: Counter[int]) -> dict[str, Any]:
    if not hist_counter:
        return {}
    maxl = max(hist_counter)
    hist = np.zeros(maxl + 1, dtype=np.int64)
    for l, n in hist_counter.items(): hist[l] = n
    n = int(hist.sum()); x = np.arange(hist.size, dtype=float)
    return {
        "read_length_observations": n,
        "read_length_mean": float(np.dot(x, hist)/n) if n else None,
        "read_length_median": inverse_ecdf_quantile(hist,.5) if n else None,
        "read_length_p10": inverse_ecdf_quantile(hist,.1) if n else None,
        "read_length_p90": inverse_ecdf_quantile(hist,.9) if n else None,
        "read_length_max": int(np.flatnonzero(hist)[-1]) if n else None,
    }


def samtools_stats_metrics(path: str | Path, hifi_bases: int | None = None, collect_rfs: bool = False, reference_length: int | None = None, reference_fasta: str | Path | None = None, known_fasta_length: int | None = None) -> tuple[dict[str, Any], dict[str, Any]]:
    data = parse_samtools_stats(path, collect_rfs=collect_rfs)
    vals: dict[str, Any] = {f"sn_{k}": v for k, v in data["sn"].items()}
    vals["rfs_detail_rows"] = data.get("rfs_detail_count", 0)
    if data["rfs_summary"]:
        vals.update(data["rfs_summary"])
    denominator = denominator_from_stats(data)
    unrestricted, recorded_fasta = stats_reference_command(data.get("command_line"))
    target = data["sn"].get("bases_inside_the_target")
    target_known = finite(target) and float(target) > 0
    # A known target denominator always takes precedence over whole FASTA length.
    # Unknown or restricted recorded commands are not evidence of a whole reference.
    length_source = "explicit" if reference_length is not None else ("samtools" if denominator else None)
    if reference_length is None and not target_known and denominator is None:
        candidate = Path(reference_fasta) if reference_fasta else (resolve_stats_reference(path, recorded_fasta) if unrestricted else None)
        if known_fasta_length is not None and unrestricted:
            reference_length = known_fasta_length
            length_source = "manifest_fasta"
        elif candidate is not None and (unrestricted or not data.get("command_line")):
            reference_length = fasta_total_length(candidate)
            length_source = "reference_fasta"
            vals["coverage_reference_fasta"] = str(candidate)
    vals["coverage_denominator_source"] = length_source

    hist, covmeta = coverage_hist_from_stats(data, reference_length)
    denominator = covmeta.get("coverage_histogram_denominator", denominator)
    vals.update(covmeta)
    if hist is not None and int(hist.sum()) > 0:
        cov = depth_metrics_from_hist(hist, scale_mad=True)
        if covmeta.get("coverage_mean_is_lower_bound"):
            cov["coverage_mean_lower_bound"] = cov.pop("coverage_mean")
            overflow_depth = covmeta["coverage_overflow_depth_lower_bound"]
            for t in (1, 5, 10, 20, 30, 50, 100):
                if t > overflow_depth:
                    cov[f"coverage_fraction_ge_{t}x"] = None
            for label in ("025", "050", "150", "200"):
                if cov[f"coverage_threshold_{label}_baseline"] > overflow_depth:
                    cov[f"coverage_fraction_ge_{label}_baseline"] = None
        if denominator is None:
            # These fractions are conditional on covered sites, not breadth.
            for key in list(cov):
                if key.startswith("coverage_fraction_ge_"):
                    cov[key + "_covered_sites"] = cov.pop(key)
            cov["coverage_zero_sites"] = None
            cov["coverage_zero_fraction"] = None
            cov["coverage_sites_total"] = None
        vals.update(cov)

    mq, mqhist = mapq_metrics(data["mapq"]); vals.update(mq)
    mapped_cigar = data["sn"].get("bases_mapped_cigar")
    vals.update(indel_metrics(data["indels"], float(mapped_cigar) if finite(mapped_cigar) else None))
    vals.update(gc_depth_metrics(data["gcd"]))
    vals.update(read_length_metrics(data["read_length"]))

    err = data["sn"].get("error_rate")
    if finite(err) and float(err) > 0:
        vals["samtools_nm_edit_qv"] = -10.0 * math.log10(float(err))
    if denominator and finite(mapped_cigar):
        vals["mapped_cigar_coverage"] = float(mapped_cigar) / denominator
    baseline = vals.get("coverage_median_positive")
    if finite(mapped_cigar) and finite(baseline) and float(baseline) > 0:
        vals["depth_implied_sequence_space"] = float(mapped_cigar) / float(baseline)
    if hifi_bases is not None and denominator and denominator > 0:
        vals["nominal_input_coverage"] = float(hifi_bases) / denominator

    detail = {"parsed": data, "coverage_hist": hist, "mapq_hist": mqhist}
    return vals, detail


def write_samtools_detail(detail: dict[str, Any], outdir: Path, prefix: str) -> None:
    data = detail["parsed"]
    hist = detail.get("coverage_hist")
    if hist is not None:
        with (outdir / f"{prefix}.coverage_histogram.tsv").open("w", encoding="utf-8") as fh:
            fh.write("depth\tsites\n")
            for d, n in enumerate(hist):
                if n: fh.write(f"{d}\t{int(n)}\n")
    with (outdir / f"{prefix}.mapq_histogram.tsv").open("w", encoding="utf-8") as fh:
        fh.write("mapq\tcount\n")
        for q, n in sorted(data["mapq"].items()): fh.write(f"{q}\t{n}\n")
    with (outdir / f"{prefix}.indel_distribution.tsv").open("w", encoding="utf-8") as fh:
        fh.write("length\tinsertions\tdeletions\n")
        for k, i, d in data["indels"]: fh.write(f"{k}\t{i}\t{d}\n")
    if data["gcd"]:
        with (outdir / f"{prefix}.gc_depth.tsv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data["gcd"][0]), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(data["gcd"])
    if data["rfs_rows"]:
        with (outdir / f"{prefix}.reference_sequences.tsv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data["rfs_rows"][0]), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(data["rfs_rows"])


def plot_stats(detail: dict[str, Any], vals: dict[str, Any], outdir: Path, prefix: str, dpi: int) -> None:
    hist = detail.get("coverage_hist")
    if hist is not None and int(hist.sum()) > 0:
        plot_coverage(hist, vals, outdir, prefix, dpi)
    mq = detail.get("mapq_hist")
    if isinstance(mq, np.ndarray) and mq.size and mq.sum():
        fig, ax = plt.subplots(figsize=(10, 6)); x = np.arange(mq.size)
        ax.bar(x, mq/mq.sum(), color=plt.get_cmap(COLORS[3])(np.linspace(.25,.85,mq.size)), width=.9)
        ax.set_xlabel("MAPQ"); ax.set_ylabel("Alignment fraction"); ax.set_title(f"{prefix}: mapping-quality distribution")
        style_ax(ax); save_svg(fig, outdir/f"{prefix}.mapq.svg", dpi)
    rows = detail["parsed"]["indels"]
    if rows:
        x = np.asarray([r[0] for r in rows]); ins=np.asarray([r[1] for r in rows]); dele=np.asarray([r[2] for r in rows])
        fig, ax = plt.subplots(figsize=(10,6)); ax.plot(x,ins,color=COLORS[1],label="Insertions"); ax.plot(x,dele,color=COLORS[0],label="Deletions")
        ax.set_yscale("log"); ax.set_xlabel("Indel length (bp)"); ax.set_ylabel("Events (log scale)"); ax.set_title(f"{prefix}: indel-size distribution"); ax.legend(frameon=False)
        style_ax(ax); save_svg(fig,outdir/f"{prefix}.indels.svg",dpi)
    gcd = detail["parsed"]["gcd"]
    if gcd:
        gc=np.asarray([r["gc_percent"] for r in gcd]); p10=np.asarray([r["depth_p10"] for r in gcd]); p25=np.asarray([r["depth_p25"] for r in gcd]); med=np.asarray([r["depth_p50"] for r in gcd]); p75=np.asarray([r["depth_p75"] for r in gcd]); p90=np.asarray([r["depth_p90"] for r in gcd])
        fig,ax=plt.subplots(figsize=(10,6)); ax.fill_between(gc,p10,p90,color=COLORS[2],alpha=.18,label="P10–P90"); ax.fill_between(gc,p25,p75,color=COLORS[1],alpha=.25,label="P25–P75"); ax.plot(gc,med,color=COLORS[0],linewidth=2,label="Median")
        ax.set_xlabel("GC (%)"); ax.set_ylabel("samtools GCD depth"); ax.set_title(f"{prefix}: GC-dependent coverage\nShaded bands: 10–90 and 25–75 percentiles", loc="left"); ax.legend(frameon=False)
        style_ax(ax); save_svg(fig,outdir/f"{prefix}.gc_depth.svg",dpi)


# -------------------------------- BAM metrics --------------------------------

def bam_metrics(path: str | Path, include_secondary: bool = False) -> dict[str, Any]:
    if pysam is None:
        raise RuntimeError("pysam is required for BAM/CRAM metrics")
    totals = Counter(); mapq = Counter(); aligned = nm = ins = dele = soft = hard = 0
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for r in bam.fetch(until_eof=True):
            totals["records"] += 1
            if r.is_unmapped:
                totals["unmapped"] += 1; continue
            totals["mapped"] += 1
            if r.is_secondary:
                totals["secondary"] += 1
                if not include_secondary: continue
            if r.is_supplementary: totals["supplementary"] += 1
            if r.is_duplicate: totals["duplicate"] += 1
            if r.is_qcfail: totals["qcfail"] += 1
            mapq[int(r.mapping_quality)] += 1
            if r.cigartuples:
                for op, ln in r.cigartuples:
                    if op in (0,7,8): aligned += ln
                    elif op == 1: ins += ln
                    elif op == 2: dele += ln
                    elif op == 4: soft += ln
                    elif op == 5: hard += ln
            try: nm += int(r.get_tag("NM"))
            except KeyError: totals["missing_nm"] += 1
    vals = {
        "alignment_records": totals["records"], "mapped_records": totals["mapped"], "unmapped_records": totals["unmapped"],
        "mapped_fraction": totals["mapped"]/max(1,totals["records"]), "secondary_records": totals["secondary"],
        "supplementary_records": totals["supplementary"], "duplicate_flag_records": totals["duplicate"], "qcfail_records": totals["qcfail"],
        "aligned_bases": aligned, "nm_total": nm, "inserted_bases": ins, "deleted_bases": dele,
        "softclipped_bases": soft, "hardclipped_bases": hard, "missing_nm_records": totals["missing_nm"],
        "nm_rate_per_aligned_base": nm/aligned if aligned else None,
        "alignment_nm_qv": (-10*math.log10(nm/aligned) if aligned and nm>0 else None),
    }
    mqvals,_ = mapq_metrics(mapq); vals.update(mqvals)
    return vals


# ------------------------------ comparison mode ------------------------------
METRIC_META: dict[str, tuple[str,str]] = {
    "assembly.total_length": ("bp", "descriptive"), "assembly.scaffold_count": ("count", "lower"),
    "assembly.largest_scaffold": ("bp", "higher"), "assembly.n50": ("bp", "higher"), "assembly.n90": ("bp", "higher"),
    "assembly.ng50": ("bp", "higher"), "assembly.ng90": ("bp", "higher"), "assembly.auN": ("bp", "higher"), "assembly.auNG": ("bp", "higher"),
    "assembly.gap_fraction": ("fraction", "lower"), "assembly.gap_count": ("count", "lower"),
    "assembly.assembly_to_genome_ratio": ("ratio", "target_1"), "assembly.gc_percent": ("%", "descriptive"),
    "stats.coverage_median_positive": ("x", "descriptive"), "stats.coverage_fraction_ge_1x": ("fraction", "higher"),
    "stats.coverage_zero_fraction": ("fraction", "lower"), "stats.coverage_robust_cv": ("ratio", "lower"),
    "stats.mapq0_fraction": ("fraction", "lower"), "stats.mapq_median": ("MAPQ", "higher"),
    "stats.samtools_nm_edit_qv": ("QV", "higher"), "stats.indel_bases_per_aligned_base": ("fraction", "lower"),
    "stats.gc_depth_robust_cv": ("ratio", "lower"), "stats.mapped_cigar_coverage": ("x", "descriptive"),
    "stats.nominal_input_coverage": ("x", "descriptive"),
}


def read_manifest(path: str | Path) -> list[dict[str,str]]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or "name" not in rows[0] or "fasta" not in rows[0]:
        raise ValueError("Manifest requires tab-delimited columns: name and fasta; optional samtools_stats and hifi_bases")
    names = [r["name"].strip() for r in rows]
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("Manifest assembly names must be nonempty and unique")
    return rows


def _one_compare_fasta(task: tuple[str,str,int|None,int]) -> tuple[str,dict[str,Any],np.ndarray]:
    name, fasta, genome_size, gap_min = task
    vals,_,lens,_ = assembly_metrics(fasta, genome_size, gap_min)
    return name, vals, lens


def comparison_metrics(manifest: str | Path, genome_size: int | None, gap_min: int, threads: int) -> tuple[list[dict[str,Any]], dict[str,np.ndarray]]:
    rows = read_manifest(manifest)
    tasks = [(r["name"], r["fasta"], genome_size, gap_min) for r in rows]
    fmetrics: dict[str,dict[str,Any]] = {}; lengths: dict[str,np.ndarray] = {}
    if threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(threads,len(tasks))) as ex:
            for name, vals, lens in ex.map(_one_compare_fasta, tasks):
                fmetrics[name]=vals; lengths[name]=lens
    else:
        for t in tasks:
            name,vals,lens=_one_compare_fasta(t); fmetrics[name]=vals; lengths[name]=lens
    out=[]
    for r in rows:
        name=r["name"]; combined: dict[str,Any]={"name":name,"fasta":r["fasta"]}
        for k,v in fmetrics[name].items(): combined[f"assembly.{k}"]=v
        st=r.get("samtools_stats","").strip()
        if st:
            hb=int(r["hifi_bases"]) if r.get("hifi_bases","").strip() else None
            sv,_=samtools_stats_metrics(st,hb, known_fasta_length=fmetrics[name]["total_length"])
            for k,v in sv.items(): combined[f"stats.{k}"]=v
            combined["samtools_stats"]=st
        out.append(combined)
    return out,lengths


def write_comparison(rows: list[dict[str,Any]], outdir: Path) -> None:
    names=[r["name"] for r in rows]
    metrics=sorted({k for r in rows for k in r if k not in {"name","fasta","samtools_stats"}})
    with (outdir/"comparison.tsv").open("w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["metric","unit","direction",*names])
        for m in metrics:
            unit,direction=METRIC_META.get(m,("","descriptive"))
            vals=[]
            for r in rows:
                v=r.get(m); vals.append("NA" if v is None or (isinstance(v,float) and not math.isfinite(v)) else v)
            w.writerow([m,unit,direction,*vals])
    with (outdir/"comparison.long.tsv").open("w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["assembly","metric","value","unit","direction"])
        for r in rows:
            for m in metrics:
                v=r.get(m)
                if v is None or (isinstance(v,float) and not math.isfinite(v)): continue
                unit,direction=METRIC_META.get(m,("","descriptive")); w.writerow([r["name"],m,v,unit,direction])


def robust_scaled(values: np.ndarray, direction: str) -> np.ndarray:
    x=values.astype(float); out=np.full_like(x, np.nan)
    mask=np.isfinite(x)
    if mask.sum()<2:
        out[mask] = 0
        return out
    med=np.median(x[mask]); mad=np.median(np.abs(x[mask]-med)); scale=1.4826*mad
    if scale<=0: scale=np.std(x[mask]) or 1.0
    z=(x-med)/scale
    if direction=="lower": z=-z
    elif direction=="target_1": z=-np.abs(x-1.0)/scale
    out[mask]=z[mask]; out[~mask]=np.nan
    return out


def plot_comparison(rows: list[dict[str,Any]], lengths: dict[str,np.ndarray], genome_size: int|None, outdir: Path, dpi:int) -> None:
    # Nx/NGx overlay
    fig,ax=plt.subplots(figsize=(11,7)); colors=plt.get_cmap("tab20")(np.linspace(0,1,max(2,len(rows))))
    for color,r in zip(colors,rows):
        name=r["name"]; lens=lengths[name]; ax.plot(range(1,101),nx_profile(lens,int(lens.sum())),color=color,linewidth=2,label=f"{name} Nx")
        if genome_size:
            ax.plot(range(1,101),nx_profile(lens,int(genome_size)),color=color,linewidth=1.2,linestyle="--",alpha=.8,label=f"{name} NGx")
    ax.set_xlabel("x (%)"); ax.set_ylabel("Length (bp)"); ax.set_title("Assembly contiguity comparison"); ax.legend(frameon=False,fontsize=8,ncol=2)
    style_ax(ax); save_svg(fig,outdir/"comparison.nx_ngx.svg",dpi)

    selected=[m for m in METRIC_META if any(finite(r.get(m)) for r in rows) and METRIC_META[m][1] != "descriptive"]
    if selected:
        matrix=[]
        for m in selected:
            vals=np.asarray([float(r[m]) if finite(r.get(m)) else np.nan for r in rows]); matrix.append(robust_scaled(vals,METRIC_META[m][1]))
        mat=np.asarray(matrix).T
        fig,ax=plt.subplots(figsize=(max(10,.55*len(selected)+4),max(5,.45*len(rows)+3)))
        im=ax.imshow(mat,aspect="auto",cmap="RdBu",vmin=-2.5,vmax=2.5)
        ax.set_xticks(range(len(selected)),[m.replace("assembly.","").replace("stats.","") for m in selected],rotation=60,ha="right")
        ax.set_yticks(range(len(rows)),[r["name"] for r in rows]); ax.set_title("Robust comparative metric heatmap")
        fig.colorbar(im,ax=ax,label="Direction-aware robust score"); save_svg(fig,outdir/"comparison.heatmap.svg",dpi)


# ---------------------------------- commands ---------------------------------

def cmd_fasta(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); prefix=a.prefix or Path(a.fasta).stem
    vals,scaffolds,lens,_=assembly_metrics(a.fasta,a.genome_size,a.gap_min)
    write_kv_tsv(out/f"{prefix}.assembly_metrics.tsv",vals); write_json(out/f"{prefix}.assembly_metrics.json",vals)
    with (out/f"{prefix}.scaffolds.tsv").open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(scaffolds[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(scaffolds)
    if not a.no_plots: plot_fasta_metrics(vals,lens,out,prefix,a.dpi)


def cmd_stats(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); prefix=a.prefix or Path(a.stats).stem
    vals,detail=samtools_stats_metrics(a.stats,a.hifi_bases,collect_rfs=a.write_reference_table, reference_length=a.reference_length, reference_fasta=a.reference_fasta)
    write_kv_tsv(out/f"{prefix}.samtools_metrics.tsv",vals); write_json(out/f"{prefix}.samtools_metrics.json",vals)
    write_samtools_detail(detail,out,prefix)
    if not a.no_plots: plot_stats(detail,vals,out,prefix,a.dpi)


def cmd_bam(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); prefix=a.prefix or Path(a.bam).stem
    vals=bam_metrics(a.bam,a.include_secondary); write_kv_tsv(out/f"{prefix}.bam_metrics.tsv",vals); write_json(out/f"{prefix}.bam_metrics.json",vals)


def cmd_compare(a: argparse.Namespace) -> None:
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    rows,lengths=comparison_metrics(a.manifest,a.genome_size,a.gap_min,a.threads)
    write_comparison(rows,out); write_json(out/"comparison.json",rows)
    if not a.no_plots: plot_comparison(rows,lengths,a.genome_size,out,a.dpi)


def build_parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(prog="genoflowSTATS",description="Genome assembly statistics and read-support QC")
    p.add_argument("--version",action="version",version=f"%(prog)s {__version__}")
    sub=p.add_subparsers(dest="cmd",required=True)

    f=sub.add_parser("fasta",help="Evaluate one assembly FASTA")
    f.add_argument("--fasta",required=True); f.add_argument("--genome-size",type=int); f.add_argument("--gap-min",type=int,default=10)
    f.add_argument("--outdir",required=True); f.add_argument("--prefix"); f.add_argument("--dpi",type=int,default=600); f.add_argument("--no-plots",action="store_true"); f.set_defaults(func=cmd_fasta)

    s=sub.add_parser("stats",help="Parse samtools stats and estimate coverage/alignment QC")
    s.add_argument("--stats",required=True); s.add_argument("--hifi-bases",type=int,help="Optional total input HiFi bases for nominal coverage")
    s.add_argument("--reference-fasta", help="Matching full-reference FASTA; calculate length from sequence lines")
    s.add_argument("--reference-length", type=int, help="Reference or target length matching the COV domain, to reconstruct zero-depth sites")
    s.add_argument("--write-reference-table",action="store_true",help="Retain/write detailed RFS rows (can be very large for target-region stats)")
    s.add_argument("--outdir",required=True); s.add_argument("--prefix"); s.add_argument("--dpi",type=int,default=600); s.add_argument("--no-plots",action="store_true"); s.set_defaults(func=cmd_stats)

    b=sub.add_parser("bam",help="Direct one-pass BAM QC via pysam")
    b.add_argument("--bam",required=True); b.add_argument("--include-secondary",action="store_true"); b.add_argument("--outdir",required=True); b.add_argument("--prefix"); b.set_defaults(func=cmd_bam)

    c=sub.add_parser("compare",help="Evaluate multiple assemblies in one run")
    c.add_argument("--manifest",required=True,help="TSV: name, fasta, optional samtools_stats, hifi_bases")
    c.add_argument("--genome-size",type=int); c.add_argument("--gap-min",type=int,default=10); c.add_argument("--threads",type=int,default=max(1,min(8,os.cpu_count() or 1)))
    c.add_argument("--outdir",required=True); c.add_argument("--dpi",type=int,default=600); c.add_argument("--no-plots",action="store_true"); c.set_defaults(func=cmd_compare)
    for command in (f, s, c):
        command.add_argument("--palette", choices=sorted(PALETTES), default="colorblind", help="Plot color palette")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    global COLORS
    args=build_parser().parse_args(argv)
    COLORS = PALETTES.get(getattr(args, "palette", "colorblind"), COLORS)
    try:
        args.func(args); return 0
    except KeyboardInterrupt:
        console.log("Interrupted"); return 130
    except Exception as exc:
        console.log(f"[ERROR] {exc}"); return 1

if __name__ == "__main__":
    raise SystemExit(main())
