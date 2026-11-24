#!/usr/bin/env python3
"""split.py - cluster-aware 50/50 train/test split per kingdom for ensemble_v2.

Sequences are grouped by MMseqs2 clusters (50% identity, 50% bidirectional coverage) so
near-identical families never sit on both sides. Clusters holding any sequence from a
held-out genome of the whole-pipeline benchmark go to TEST. The remaining clusters are
split with StratifiedGroupKFold (2 folds, stratified by order).

usage: split.py --truth pools/fungi.truth.tsv --clusters split/fungi_cluster.tsv
                --out split/fungi [--seed 7]
"""
import argparse
import csv
from collections import Counter

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--clusters", required=True, help="mmseqs <prefix>_cluster.tsv (rep, member)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    truth = list(csv.DictReader(open(a.truth), delimiter="\t"))
    clu = {}
    for line in open(a.clusters):
        rep, mem = line.split()[:2]
        clu[mem] = rep
    ids = [r["family_id"] for r in truth]
    groups = np.array([clu.get(i, i) for i in ids])
    y = np.array([r["order"] for r in truth])
    forced = {groups[i] for i, r in enumerate(truth) if r["heldout_species"] == "1"}
    is_forced = np.array([g in forced for g in groups])

    split = np.array(["train"] * len(ids), dtype=object)
    split[is_forced] = "test"
    free = np.where(~is_forced)[0]
    sgkf = StratifiedGroupKFold(n_splits=2, shuffle=True, random_state=a.seed)
    _, test_idx = next(sgkf.split(np.zeros(len(free)), y[free], groups[free]))
    split[free[test_idx]] = "test"

    with open(a.out + ".split.tsv", "w") as fh:
        fh.write("family_id\tcluster\tsplit\torder\tsource\theldout_species\n")
        for i, r in enumerate(truth):
            fh.write("\t".join([ids[i], groups[i], split[i], y[i], r["source"], r["heldout_species"]]) + "\n")
    for part in ("train", "test"):
        with open("%s.%s.ids" % (a.out, part), "w") as fh:
            fh.write("\n".join(ids[i] for i in range(len(ids)) if split[i] == part) + "\n")

    n_clu = len(set(groups))
    print("%s: %d sequences, %d clusters (%.1f seq/cluster); forced-to-test clusters %d (%d seqs)"
          % (a.truth, len(ids), n_clu, len(ids) / n_clu, len(forced), int(is_forced.sum())))
    for part in ("train", "test"):
        sel = split == part
        c = Counter(y[sel])
        src = Counter(r["source"] for i, r in enumerate(truth) if sel[i])
        print("  %-5s %6d  %s  %s" % (part, int(sel.sum()), dict(sorted(c.items())), dict(src)))
    # leakage check: no cluster on both sides
    tr = set(groups[split == "train"]); te = set(groups[split == "test"])
    print("  clusters in both splits: %d" % len(tr & te))


if __name__ == "__main__":
    main()

