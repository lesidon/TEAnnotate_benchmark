#!/usr/bin/env python3
"""
leakage_filter.py - build a leakage-controlled TE test set.

Several classifiers in the ensemble were trained on Repbase and/or Dfam
(TEclass2 and CREATE were trained on Dfam directly; Terrier/DeepTE/NeuralTE on
Repbase, which overlaps Dfam at the family level). Evaluating on the whole Dfam
partition therefore leaks training data into the test set and inflates the
metrics.

This script removes that overlap: it searches the candidate Dfam families
(--query, from dfam_partitioner.py) against one or more training corpora
(--training) with MMseqs2, and drops any candidate that matches a training
sequence above an identity + coverage threshold. What survives is a
held-out set the models have (approximately) not seen.

Outputs (to --outdir, stem = --prefix):
  <prefix>.heldout.fasta        candidates with NO qualifying training hit
  <prefix>.heldout.ids          their ids - feed to evaluate_classification.py --keep-list
  <prefix>.leaked.fasta         removed candidates            (--write-leaked)
  <prefix>.leaked.ids           their ids - evaluate_classification.py --exclude-list
  <prefix>.leakage_report.tsv   per candidate: status + best training hit
  <prefix>.leakage_summary.tsv  per source removal counts + totals

Requires MMseqs2 on PATH (or --mmseqs /path/to/mmseqs). Inputs must be plain
(uncompressed) nucleotide FASTA - MMseqs2 cannot read .gz.

Example
-------
  python leakage_filter.py \
      --query dfam_partitions/dfam_fungi.fasta \
      --training repbase=RepBase29.10.fasta dfam38=Dfam_3.8.fasta pgsb=mipsREdat.fasta \
      --identity 0.7 --coverage 0.5 --threads 16 \
      --outdir leakage_out --prefix dfam_fungi --write-leaked
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from typing import Dict, List, Optional, Tuple

import te_label_utils as tu

LOG = tu.setup_logger("leakage_filter")

# MMseqs2 easy-search custom output columns (see --format-output).
FMT_COLS = ["query", "target", "fident", "qcov", "tcov", "evalue", "bits"]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def parse_sources(items: List[str]) -> List[Tuple[str, str]]:
    """Parse NAME=PATH / PATH training entries into (name, path) pairs."""
    out: List[Tuple[str, str]] = []
    for it in items:
        if "=" in it:
            name, path = it.split("=", 1)
        else:
            name = os.path.splitext(os.path.basename(it))[0]
            path = it
        if not os.path.isfile(path):
            raise SystemExit(f"ERROR: training FASTA not found: {path}")
        out.append((name, path))
    return out


def iter_fasta_ids(path: str):
    """Stream (id, record_length_hint) - id is the first header token."""
    from Bio import SeqIO
    with tu.open_text(path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            yield rec.id.split("#", 1)[0]


def run_mmseqs(mmseqs: str, query: str, target: str, out_m8: str, tmpdir: str,
               identity: float, coverage: float, sensitivity: float,
               threads: int) -> None:
    cmd = [
        mmseqs, "easy-search", query, target, out_m8, tmpdir,
        "--search-type", "3",              # nucleotide-nucleotide
        "-s", str(sensitivity),
        "--min-seq-id", str(identity),
        "-c", str(coverage),
        "--cov-mode", "2",                 # coverage of the query (test) sequence
        "--threads", str(threads),
        "--format-output", ",".join(FMT_COLS),
        "-v", "1",
    ]
    LOG.info("mmseqs: %s", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"ERROR: mmseqs easy-search failed (rc={proc.returncode})")


def parse_hits(m8: str, identity: float, coverage: float
               ) -> Dict[str, Tuple[str, float, float, float, float]]:
    """Best qualifying hit per query id from an m8 file (streaming).

    Returns {query_id: (target, fident, qcov, tcov, bits)} keeping the max-bits
    hit that meets identity + coverage (re-checked exactly, since MMseqs's
    --min-seq-id prefilter is approximate).
    """
    best: Dict[str, Tuple[str, float, float, float, float]] = {}
    if not os.path.isfile(m8):
        return best
    with open(m8, "r", encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(FMT_COLS):
                continue
            q = f[0].split("#", 1)[0]
            target = f[1]
            try:
                fident = float(f[2]); qcov = float(f[3])
                tcov = float(f[4]); bits = float(f[6])
            except ValueError:
                continue
            if fident < identity or max(qcov, tcov) < coverage:
                continue
            prev = best.get(q)
            if prev is None or bits > prev[4]:
                best[q] = (target, fident, qcov, tcov, bits)
    return best


def write_fasta_subset(query: str, keep_out: str, leaked_out: Optional[str],
                       leaked_ids: set) -> Tuple[int, int]:
    """Stream the query FASTA, splitting records into held-out vs leaked."""
    from Bio import SeqIO
    n_keep = n_leak = 0
    lk_fh = open(leaked_out, "w", encoding="utf-8") if leaked_out else None
    try:
        with tu.open_text(query) as ih, open(keep_out, "w", encoding="utf-8") as kh:
            for rec in SeqIO.parse(ih, "fasta"):
                base = rec.id.split("#", 1)[0]
                block = f">{rec.description}\n{str(rec.seq)}\n"
                if base in leaked_ids:
                    n_leak += 1
                    if lk_fh:
                        lk_fh.write(block)
                else:
                    n_keep += 1
                    kh.write(block)
    finally:
        if lk_fh:
            lk_fh.close()
    return n_keep, n_leak


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Remove training-set overlap from a Dfam TE test partition "
                    "using MMseqs2 (leakage control).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--query", required=True,
                   help="Candidate test FASTA (e.g. dfam_partitioner output).")
    p.add_argument("--training", required=True, nargs="+", metavar="NAME=FILE",
                   help="Training corpora to screen against (Repbase, Dfam, PGSB...).")
    p.add_argument("--identity", type=float, default=0.7,
                   help="Min fraction identity for a hit to count as leakage.")
    p.add_argument("--coverage", type=float, default=0.5,
                   help="Min alignment coverage (max of query/target) for leakage.")
    p.add_argument("--sensitivity", type=float, default=7.5,
                   help="MMseqs2 -s sensitivity.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--mmseqs", default="mmseqs", help="mmseqs executable.")
    p.add_argument("--tmpdir", default=None, help="Parent dir for MMseqs tmp.")
    p.add_argument("--keep-tmp", action="store_true", help="Keep MMseqs tmp dirs.")
    p.add_argument("--write-leaked", action="store_true",
                   help="Also write the removed (leaked) sequences to a FASTA.")
    p.add_argument("--outdir", default="leakage_out")
    p.add_argument("--prefix", default=None,
                   help="Output stem (default: query file stem).")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")

    mmseqs = shutil.which(args.mmseqs) or (
        args.mmseqs if os.path.isfile(args.mmseqs) else None)
    if not mmseqs:
        raise SystemExit(f"ERROR: mmseqs not found ('{args.mmseqs}'). Install "
                         "MMseqs2 or pass --mmseqs /path/to/mmseqs.")
    if not os.path.isfile(args.query):
        raise SystemExit(f"ERROR: query FASTA not found: {args.query}")

    sources = parse_sources(args.training)
    tu.ensure_dir(args.outdir)
    prefix = args.prefix or os.path.splitext(os.path.basename(args.query))[0]
    tmp_root = tempfile.mkdtemp(prefix="leakage_", dir=args.tmpdir)

    # All candidate ids (complete audit, even those with no hit).
    all_ids = list(iter_fasta_ids(args.query))
    LOG.info("Loaded %d candidate families from %s", len(all_ids), args.query)

    # Search each training source; collect best qualifying hit per source.
    per_source_best: Dict[str, Dict[str, tuple]] = {}
    source_hit_counts: "Counter[str]" = Counter()
    try:
        for name, path in sources:
            m8 = os.path.join(tmp_root, f"{name}.m8")
            stmp = os.path.join(tmp_root, f"tmp_{name}")
            run_mmseqs(mmseqs, args.query, path, m8, stmp,
                       args.identity, args.coverage, args.sensitivity, args.threads)
            best = parse_hits(m8, args.identity, args.coverage)
            per_source_best[name] = best
            source_hit_counts[name] = len(best)
            LOG.info("  %-12s: %d candidate families matched", name, len(best))
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmp_root, ignore_errors=True)

    # Aggregate: a candidate leaks if it matches ANY source. Keep the globally
    # best hit (max bits) for the report.
    leaked_best: Dict[str, Tuple[str, str, float, float, float, float]] = {}
    for name, best in per_source_best.items():
        for q, (target, fident, qcov, tcov, bits) in best.items():
            prev = leaked_best.get(q)
            if prev is None or bits > prev[5]:
                leaked_best[q] = (name, target, fident, qcov, tcov, bits)
    leaked_ids = set(leaked_best)

    # --- Report (every candidate) ------------------------------------------
    report = os.path.join(args.outdir, f"{prefix}.leakage_report.tsv")
    with open(report, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "status", "best_source", "best_target",
                    "fident", "qcov", "tcov", "bits", "n_sources_hit"])
        for fid in all_ids:
            if fid in leaked_best:
                src, tgt, fi, qc, tc, bt = leaked_best[fid]
                n_src = sum(1 for b in per_source_best.values() if fid in b)
                w.writerow([fid, "removed", src, tgt, f"{fi:.4f}", f"{qc:.4f}",
                            f"{tc:.4f}", f"{bt:.1f}", n_src])
            else:
                w.writerow([fid, "kept", "", "", "", "", "", "", 0])

    # --- FASTA + id-list outputs -------------------------------------------
    keep_fa = os.path.join(args.outdir, f"{prefix}.heldout.fasta")
    leak_fa = os.path.join(args.outdir, f"{prefix}.leaked.fasta") \
        if args.write_leaked else None
    n_keep, n_leak = write_fasta_subset(args.query, keep_fa, leak_fa, leaked_ids)

    with open(os.path.join(args.outdir, f"{prefix}.heldout.ids"), "w",
              encoding="utf-8") as fh:
        for fid in all_ids:
            if fid not in leaked_ids:
                fh.write(fid + "\n")
    with open(os.path.join(args.outdir, f"{prefix}.leaked.ids"), "w",
              encoding="utf-8") as fh:
        for fid in all_ids:
            if fid in leaked_ids:
                fh.write(fid + "\n")

    # ── Summary ────────────────────────────────────────────────────────────
    total = len(all_ids)
    summary = os.path.join(args.outdir, f"{prefix}.leakage_summary.tsv")
    with open(summary, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["query", args.query])
        w.writerow(["identity_threshold", args.identity])
        w.writerow(["coverage_threshold", args.coverage])
        w.writerow(["n_candidates", total])
        w.writerow(["n_heldout", n_keep])
        w.writerow(["n_leaked", n_leak])
        w.writerow(["pct_leaked", f"{(100.0 * n_leak / total) if total else 0:.2f}"])
        for name, _ in sources:
            w.writerow([f"matched_by_{name}", source_hit_counts[name]])

    LOG.info("Held-out: %d | leaked: %d | %.1f%% removed",
             n_keep, n_leak, (100.0 * n_leak / total) if total else 0.0)
    LOG.info("Held-out FASTA: %s", keep_fa)
    LOG.info("Feed %s.heldout.ids to evaluate_classification.py --keep-list", prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
