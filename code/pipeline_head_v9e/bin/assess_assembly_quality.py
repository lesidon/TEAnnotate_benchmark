#!/usr/bin/env python3
"""
Parse seqkit stats -a -T output and emit a JSON file with assembly statistics
and an auto-derived quality tier (draft | chromosome).

Quality tiers (gate satellite/HOR Stage B and the karyoplot):
  chromosome  — assembly is chromosome-scale by EITHER rule below
  draft       — neither rule met

The tier is decided in a genome-size-aware way, so compact genomes (e.g. fungi)
are not penalised by an absolute cutoff calibrated for large plant/animal
genomes:

  1. Absolute rule:  contig N50 >= --n50-floor-mb           (default 10 Mb)
                     Catches large genomes with large chromosomes.
  2. Relative rule:  contig N50 >= --n50-fraction * genome  AND  N50 >= 1 Mb
                     A chromosome-scale assembly of a k-chromosome genome has
                     N50 ~ genome/k, i.e. N50/genome ~ 1/k, so the default
                     fraction 0.05 captures organisms with up to ~20 chromosomes.
                     Example: Verticillium alfalfae (~32 Mb / 8 chromosomes) has
                     contig N50 ~4 Mb (~0.13 of the genome) -> chromosome.

The --override flag bypasses heuristic derivation and forces a specific tier.
Users can also set params.assembly_quality in nextflow.config to override.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

# Defaults for the auto-tier heuristic (overridable via CLI / nextflow params).
DEFAULT_N50_FLOOR_MB = 10.0   # Absolute contig N50 (Mb) -> chromosome-scale.
DEFAULT_N50_FRACTION = 0.05   # Relative: N50 >= fraction * genome -> chromosome-scale.
# Absolute minimum N50 for the relative rule, so a tiny or highly fragmented
# assembly is never promoted to chromosome-scale on the ratio alone.
RELATIVE_RULE_MIN_N50_BP = 1_000_000


def parse_seqkit_stats(stats_file: str) -> dict:
    """Read seqkit stats -a -T output and return a flat dict of values."""
    with open(stats_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    if not rows:
        sys.exit(f"ERROR: seqkit stats file '{stats_file}' contains no data rows.")

    # seqkit stats on a single file gives exactly one data row
    row = rows[0]

    def safe_int(v):
        try:
            return int(str(v).replace(",", ""))
        except (ValueError, TypeError):
            return 0

    def safe_float(v):
        try:
            return float(str(v).replace(",", "").rstrip("%"))
        except (ValueError, TypeError):
            return 0.0

    # seqkit -a column names (case may vary across versions; handle both)
    col = {k.strip().lower(): v for k, v in row.items()}

    return {
        "num_seqs" : safe_int(col.get("num_seqs", 0)),
        "sum_len"  : safe_int(col.get("sum_len",  0)),
        "n50"      : safe_int(col.get("n50",       0)),
        "min_len"  : safe_int(col.get("min_len",   0)),
        "max_len"  : safe_int(col.get("max_len",   0)),
        "avg_len"  : safe_float(col.get("avg_len", 0.0)),
        "gc"       : safe_float(col.get("gc(%)",   col.get("gc", 0.0))),
    }


def derive_quality(stats: dict, override: str, n50_floor_bp: int, n50_fraction: float):
    """Return (quality_tier, rule_applied), respecting the override if provided."""
    override = (override or "").strip()
    if override and override.lower() != "auto":
        valid = {"draft", "chromosome", "t2t"}
        if override.lower() not in valid:
            sys.exit(
                f"ERROR: --override '{override}' is not a valid quality tier. "
                f"Choose from: {', '.join(sorted(valid))}."
            )
        return override.lower().replace("t2t", "T2T"), "override"

    n50    = stats["n50"]
    genome = stats["sum_len"]

    # Rule 1 — absolute: large chromosomes (plant/animal genomes).
    if n50 >= n50_floor_bp:
        return "chromosome", "absolute_n50"
    # Rule 2 — relative: chromosome-scale compact genomes (e.g. fungi) whose
    # chromosomes sit well below the absolute floor.
    if genome > 0 and n50 >= RELATIVE_RULE_MIN_N50_BP and (n50 / genome) >= n50_fraction:
        return "chromosome", "relative_n50"
    return "draft", "below_thresholds"


def main():
    parser = argparse.ArgumentParser(
        description="Derive assembly quality tier from seqkit stats output."
    )
    parser.add_argument(
        "--stats",
        required=True,
        help="Path to seqkit stats -a -T TSV file.",
    )
    parser.add_argument(
        "--override",
        default="",
        help=(
            "Force quality tier (draft | chromosome | T2T). "
            "Empty or 'auto' = use the N50 heuristic."
        ),
    )
    parser.add_argument(
        "--n50-floor-mb",
        type=float,
        default=DEFAULT_N50_FLOOR_MB,
        help=(
            "Absolute contig N50 (Mb) at/above which an assembly is "
            "chromosome-scale. Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--n50-fraction",
        type=float,
        default=DEFAULT_N50_FRACTION,
        help=(
            "Relative rule: contig N50 >= fraction * genome size (and >= 1 Mb) "
            "also qualifies as chromosome-scale. Rescues compact genomes such as "
            "fungi. Default: %(default)s (~1/fraction chromosomes)."
        ),
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output JSON file path.",
    )
    args = parser.parse_args()

    stats = parse_seqkit_stats(args.stats)
    n50_floor_bp = int(args.n50_floor_mb * 1_000_000)
    quality, rule = derive_quality(stats, args.override, n50_floor_bp, args.n50_fraction)

    ratio = round(stats["n50"] / stats["sum_len"], 4) if stats["sum_len"] else 0.0

    result = {
        "assembly_quality"    : quality,
        "quality_rule"        : rule,
        "genome_size_mb"      : round(stats["sum_len"] / 1_000_000, 2),
        "num_seqs"            : stats["num_seqs"],
        "n50_mb"              : round(stats["n50"]     / 1_000_000, 3),
        "n50_to_genome_ratio" : ratio,
        "max_seq_mb"          : round(stats["max_len"] / 1_000_000, 3),
        "gc_percent"          : round(stats["gc"], 2),
    }

    Path(args.out).write_text(json.dumps(result, indent=2))

    print(
        f"[assess_assembly_quality] {Path(args.stats).stem}: "
        f"quality={quality} (rule={rule}), N50={result['n50_mb']} Mb "
        f"({ratio:.3f} of genome), size={result['genome_size_mb']} Mb, "
        f"seqs={result['num_seqs']}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
