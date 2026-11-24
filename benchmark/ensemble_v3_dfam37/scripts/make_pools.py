#!/usr/bin/env python3
"""make_pools.py - per-kingdom leakage-controlled pools for ensemble_v2.

plant / animal: PanTEon v1.5.1 held-out sets (already leakage-filtered against RepBase
29.10 + Dfam 3.8 curated at 70% identity / 50% coverage).
fungi: the PanTEon fungi held-out set + MycoMobilome v1.1 protein-evidence consensi that
passed the same filter.

Writes pools/<kingdom>.fasta (plain ids; the classifiers' input) and
pools/<kingdom>.truth.tsv (family_id, species, kingdom, source, label_raw, order,
superfamily, heldout_species). Truth labels without a TE order are dropped.
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ens_labels import NON_TE, is_heldout, load_label_overrides, order_sf  # noqa: E402


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0].split("#")[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panteon-leakage-dir", required=True)
    ap.add_argument("--panteon-tsv-dir", required=True)
    ap.add_argument("--mm-fasta", required=True)
    ap.add_argument("--mm-tsv", required=True)
    ap.add_argument("--mm-heldout-ids", required=True)
    ap.add_argument("--label-map", default=None)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    load_label_overrides(a.label_map)
    os.makedirs(a.outdir, exist_ok=True)

    for kingdom in ("plant", "animal", "fungi"):
        rows = {}
        for r in csv.DictReader(open(os.path.join(a.panteon_tsv_dir, "panteon_%s.tsv" % kingdom)),
                                delimiter="\t"):
            rows[r["family_id"]] = ("PanTEon_v1.5.1", r["species"], r["label"].split("#", 1)[-1])
        seqs = [(fid, s) for fid, s in read_fasta(
            os.path.join(a.panteon_leakage_dir, "panteon_%s.heldout.fasta" % kingdom))]
        if kingdom == "fungi":
            keep = set(l.strip() for l in open(a.mm_heldout_ids) if l.strip())
            for r in csv.DictReader(open(a.mm_tsv), delimiter="\t"):
                if r["family_id"] in keep:
                    rows[r["family_id"]] = ("MycoMobilome_v1.1", r["species"], r["label"])
            seqs += [(fid, s) for fid, s in read_fasta(a.mm_fasta) if fid in keep]
        stats, orders, by_src = Counter(), Counter(), defaultdict(Counter)
        with open(os.path.join(a.outdir, kingdom + ".fasta"), "w") as fa, \
             open(os.path.join(a.outdir, kingdom + ".truth.tsv"), "w") as tv:
            tv.write("family_id\tspecies\tkingdom\tsource\tlabel_raw\torder\tsuperfamily\theldout_species\n")
            for fid, s in seqs:
                if fid not in rows:
                    stats["no_truth"] += 1
                    continue
                src, sp, raw = rows[fid]
                o, sf = order_sf(raw)
                if o in NON_TE:
                    stats["non_te"] += 1
                    continue
                ho = is_heldout(sp)
                stats["heldout_species"] += ho
                stats["kept"] += 1
                orders[o] += 1
                by_src[src][o] += 1
                fa.write(">%s\n" % fid)
                for i in range(0, len(s), 80):
                    fa.write(s[i:i + 80] + "\n")
                tv.write("\t".join([fid, sp, kingdom, src, raw, o, sf, "1" if ho else "0"]) + "\n")
        print("%s: %s" % (kingdom, dict(stats)))
        for src, c in by_src.items():
            print("   %s: %s" % (src, c.most_common()))


if __name__ == "__main__":
    main()

