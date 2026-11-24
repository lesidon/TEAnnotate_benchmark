#!/usr/bin/env python3
"""Proxy TE-annotation metrics for genomes without a curated gold-standard library.

Given each method's RepeatMasker ``.out`` for one genome, compute:

* coverage %      = masked bp / genome bp                 (per method)
* unique bp       = bp masked by a method and no other    (per method)
* pairwise Jaccard concordance = |A∩B| / |A∪B| of masked bp (method x method)
* LAI             = whole-genome LTR Assembly Index         (per genome; from 05)

Interval arithmetic is pure-python (merge + two-pointer intersection) — no bedtools.

Outputs (into --outdir):
    <species>.coverage.tsv     kingdom species method masked_bp genome_bp coverage_pct unique_bp
    <species>.concordance.tsv  species method_a method_b intersect_bp union_bp jaccard
    <species>.lai.tsv          species lai
"""
import argparse
import gzip
import os
import sys


# ── FASTA length (streaming, gzip-aware; no Biopython dep) ──────────────────
def genome_length(path):
    op = gzip.open if path.endswith(".gz") else open
    total = 0
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == ">":
                continue
            total += len(line.strip())
    return total


# ── RepeatMasker .out -> {seq: [(begin, end), ...]} ─────────────────────────
def read_rm_out(path):
    ivs = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 7:
                continue
            try:
                int(f[0])                      # SW score => a data line
                begin = int(f[5]); end = int(f[6])
            except ValueError:
                continue
            seq = f[4]
            if end < begin:
                begin, end = end, begin
            ivs.setdefault(seq, []).append((begin, end))
    return ivs


# ── interval ops (inclusive coordinates) ────────────────────────────────────
def merge(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= out[-1][1] + 1:                # overlapping or contiguous
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


def total_bp(merged):
    return sum(e - s + 1 for s, e in merged)


def intersect_bp(a, b):
    """bp overlap between two merged interval lists (two-pointer)."""
    i = j = t = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if e >= s:
            t += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t


def merged_by_seq(ivs):
    return {seq: merge(v) for seq, v in ivs.items()}


def coverage_bp(mbyseq):
    return sum(total_bp(m) for m in mbyseq.values())


def pair_intersect_bp(m1, m2):
    return sum(intersect_bp(m1[s], m2[s]) for s in m1.keys() & m2.keys())


def parse_out_args(pairs):
    methods = {}
    for p in pairs:
        if "=" not in p:
            sys.exit(f"[proxy_metrics] --out expects method=path, got '{p}'")
        name, path = p.split("=", 1)
        if not os.path.isfile(path):
            sys.stderr.write(f"[proxy_metrics] WARNING: missing {path} (method {name}) — skipped\n")
            continue
        methods[name] = merged_by_seq(read_rm_out(path))
    return methods


def read_lai(path):
    if not path or not os.path.isfile(path):
        return "NA"
    with open(path, encoding="utf-8", errors="replace") as fh:
        val = fh.read().strip().split()
    return val[-1] if val else "NA"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="genome FASTA (.fa/.fna[.gz])")
    ap.add_argument("--species", required=True)
    ap.add_argument("--kingdom", required=True)
    ap.add_argument("--out", nargs="+", required=True,
                    help="method=RepeatMasker.out  (repeatable)")
    ap.add_argument("--lai-file", default=None, help="text file with the whole-genome LAI")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    gbp = genome_length(args.genome)
    if gbp == 0:
        sys.exit(f"[proxy_metrics] genome length 0 for {args.genome}")
    methods = parse_out_args(args.out)
    if not methods:
        sys.exit("[proxy_metrics] no readable .out files")

    names = sorted(methods)
    cov = {m: coverage_bp(methods[m]) for m in names}

    # unique bp: coverage minus overlap with the union of all other methods.
    def union_of_others(target):
        acc = {}
        for m in names:
            if m == target:
                continue
            for seq, mv in methods[m].items():
                acc.setdefault(seq, []).extend((s, e) for s, e in mv)
        return {seq: merge(v) for seq, v in acc.items()}

    unique = {}
    for m in names:
        others = union_of_others(m)
        overlap = sum(intersect_bp(methods[m][s], others.get(s, [])) for s in methods[m])
        unique[m] = cov[m] - overlap

    cov_path = os.path.join(args.outdir, f"{args.species}.coverage.tsv")
    with open(cov_path, "w", encoding="utf-8") as fh:
        fh.write("kingdom\tspecies\tmethod\tmasked_bp\tgenome_bp\tcoverage_pct\tunique_bp\n")
        for m in names:
            fh.write(f"{args.kingdom}\t{args.species}\t{m}\t{cov[m]}\t{gbp}\t"
                     f"{100.0 * cov[m] / gbp:.4f}\t{unique[m]}\n")

    con_path = os.path.join(args.outdir, f"{args.species}.concordance.tsv")
    with open(con_path, "w", encoding="utf-8") as fh:
        fh.write("species\tmethod_a\tmethod_b\tintersect_bp\tunion_bp\tjaccard\n")
        for i, a in enumerate(names):
            for b in names[i:]:
                inter = cov[a] if a == b else pair_intersect_bp(methods[a], methods[b])
                union = cov[a] + cov[b] - inter
                jac = (inter / union) if union else 0.0
                fh.write(f"{args.species}\t{a}\t{b}\t{inter}\t{union}\t{jac:.4f}\n")

    lai_path = os.path.join(args.outdir, f"{args.species}.lai.tsv")
    with open(lai_path, "w", encoding="utf-8") as fh:
        fh.write("species\tlai\n")
        fh.write(f"{args.species}\t{read_lai(args.lai_file)}\n")

    sys.stderr.write(f"[proxy_metrics] {args.species}: {len(names)} methods, "
                     f"genome {gbp:,} bp -> {cov_path}\n")


if __name__ == "__main__":
    main()
