#!/usr/bin/env python3
"""agree_check.py <sp>... - for each consensus the v10 vote relabelled from Unknown: how many
models of the kingdom's two-step Order set called the winning Order, and the share of its
masked bases on reference regions of the same class (vs any reference TE / unannotated)."""
import os, sys
from collections import defaultdict
import numpy as np
from vote_v10_common import ARMS, E, V, SUBSET, paths

sys.path.insert(0, E); sys.path.insert(0, f"{V}/bin")
import libtest_harmonised as LH  # noqa: E402
import ensemble_vote as PV       # noqa: E402
TOP = {"LTR": "LTR", "LINE": "nonLTR", "SINE": "nonLTR", "TIR": "TIR", "MITE": "TIR", "RC": "Helitron"}


def top(lab):
    c = LH.canon(lab)
    return c if c in ("Unknown", "Unspecified") else TOP.get(c.split("/")[0], "Unknown")


tot = defaultdict(lambda: [0, 0, 0])   # (class, n_agree) -> [n, bp, same-class bp]
for sp in sys.argv[1:]:
    _, size = LH.sizes(sp)
    ref = {}
    for i, l in enumerate(open(f"{E}/{sp}/reference/{sp}.fna.out")):
        f = l.split()
        if i >= 3 and len(f) >= 11 and f[4] in size:
            ref.setdefault(top(f[10]), {}).setdefault(f[4], np.zeros(size[f[4]], bool))[int(f[5]) - 1:int(f[6])] = True
    for arm in ARMS:
        d = f"{V}/apply/{sp}/{arm}"
        if not os.path.exists(f"{d}/relabel.tsv"):
            continue
        org = paths(sp, arm)[1]
        ow = PV.load_two_step(org)["order"]["weights"]
        ch = {}
        for i, l in enumerate(open(f"{d}/relabel.tsv")):
            f = l.rstrip("\n").split("\t")
            if i and f[2] and f[2] != f[1] and f[1].lower().startswith("unknown"):
                ch[f[0]] = f[2]
        hdr, votes = None, {}
        for l in open(f"{d}/v10_voting_table.tsv"):
            f = l.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [h.lower() for h in f]; continue
            votes[f[0]] = dict(zip(hdr, f))
        iv = defaultdict(lambda: defaultdict(list))
        for i, l in enumerate(open(paths(sp, arm)[4])):
            f = l.split()
            if i >= 3 and len(f) >= 11 and f[9] in ch and f[4] in size:
                iv[f[9]][f[4]].append((int(f[5]) - 1, int(f[6])))
        for name, new in ch.items():
            o = PV.vote_order(new)
            n_agree = sum(1 for m in ow if PV.vote_order(votes.get(name, {}).get(m, "NA")) == o
                          and votes[name][m] not in ("NA", ""))
            cls = top(new)
            bp = same = 0
            for c, ivs in iv[name].items():
                m = np.zeros(size[c], bool)
                for s, e in ivs:
                    m[s:e] = True
                bp += int(m.sum())
                if c in ref.get(cls, {}):
                    same += int((m & ref[cls][c]).sum())
            tot[(cls, min(n_agree, 3))][0] += 1
            tot[(cls, min(n_agree, 3))][1] += bp
            tot[(cls, min(n_agree, 3))][2] += same
            print(f"{sp:24s} {arm} {name[:40]:40s} -> {new:18s} agree={n_agree}  masked {bp/1e3:8.1f} kb  "
                  f"same class {100*same/max(1,bp):5.1f}%")
print("\nall genomes and arms: class, models agreeing on the Order -> n, masked Mb, same-class share")
for (cls, n), (k, bp, same) in sorted(tot.items()):
    print(f"  {cls:9s} agree {n if n < 3 else '3+'}: n={k:3d}  {bp/1e6:6.3f} Mb  same class {100*same/max(1,bp):5.1f}%")

