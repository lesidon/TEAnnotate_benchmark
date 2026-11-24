#!/usr/bin/env python3
"""diag_relabel.py <sp>... - where do the bases masked by consensi whose label the v10 vote
changed fall in the curated reference annotation (harmonised classes)? Per changed consensus
group (v9 label -> v10 label): bp on reference regions of the SAME class, of ANOTHER TE class,
on reference Unknown/Unspecified regions (the reference cannot judge), and outside any
reference annotation."""
import os, sys
from collections import defaultdict
import numpy as np
from vote_v10_common import ARMS, E, V, paths

sys.path.insert(0, E)
import libtest_harmonised as LH  # noqa: E402

TOP = {"LTR": "LTR", "LINE": "nonLTR", "SINE": "nonLTR", "TIR": "TIR", "MITE": "TIR", "RC": "Helitron"}


def top(lab):
    c = LH.canon(lab)
    if c in ("Unknown", "Unspecified"):
        return c
    return TOP.get(c.split("/")[0], "Unknown")


for sp in sys.argv[1:]:
    _, size = LH.sizes(sp)
    ref = {}
    for i, l in enumerate(open(f"{E}/{sp}/reference/{sp}.fna.out")):
        f = l.split()
        if i < 3 or len(f) < 11 or f[4] not in size:
            continue
        t = top(f[10])
        m = ref.setdefault(t, {})
        a = m.setdefault(f[4], np.zeros(size[f[4]], bool))
        a[int(f[5]) - 1:int(f[6])] = True
    for arm in ARMS:
        rel = f"{V}/apply/{sp}/{arm}/relabel.tsv"
        out = paths(sp, arm)[4]
        if not os.path.exists(rel):
            continue
        ch = {}
        for i, l in enumerate(open(rel)):
            f = l.rstrip("\n").split("\t")
            if i and f[2] and f[2] != f[1]:
                ch[f[0]] = (top(f[1]) if f[1] else "Unknown", top(f[2]), f[1], f[2])
        if not ch:
            continue
        # per changed consensus: union of its hit intervals
        iv = defaultdict(lambda: defaultdict(list))
        for i, l in enumerate(open(out)):
            f = l.split()
            if i >= 3 and len(f) >= 11 and f[9] in ch and f[4] in size:
                iv[f[9]][f[4]].append((int(f[5]) - 1, int(f[6])))
        agg = defaultdict(lambda: defaultdict(int))
        for name, chroms in iv.items():
            old, new, _, _ = ch[name]
            key = f"{old}->{new}"
            for c, ivs in chroms.items():
                m = np.zeros(size[c], bool)
                for s, e in ivs:
                    m[s:e] = True
                tot = int(m.sum())
                same = int((m & ref.get(new, {}).get(c, np.zeros(size[c], bool))).sum()) if new not in ("Unknown",) else 0
                unk = np.zeros(size[c], bool)
                for t in ("Unknown", "Unspecified"):
                    if c in ref.get(t, {}):
                        unk |= ref[t][c]
                anyte = np.zeros(size[c], bool)
                for t, mm in ref.items():
                    if t not in ("Unknown", "Unspecified") and c in mm:
                        anyte |= mm[c]
                agg[key]["n"] += 0
                agg[key]["bp"] += tot
                agg[key]["same_class"] += same
                agg[key]["other_te_class"] += int((m & anyte).sum()) - same
                agg[key]["ref_unknown"] += int((m & unk & ~anyte).sum())
                agg[key]["unannotated"] += int((m & ~anyte & ~unk).sum())
        for name in iv:
            agg[f"{ch[name][0]}->{ch[name][1]}"]["n"] += 1
        print(f"\n{sp} ({'with' if arm == 'w' else 'w/o'} TEtrimmer): changed consensi by class move "
              f"(Mb of their masked bases vs the reference)")
        for key, d in sorted(agg.items(), key=lambda kv: -kv[1]["bp"]):
            print(f"  {key:22s} n={d['n']:3d}  masked {d['bp']/1e6:6.3f}  same class {d['same_class']/1e6:6.3f}  "
                  f"other TE class {d['other_te_class']/1e6:6.3f}  ref Unknown {d['ref_unknown']/1e6:6.3f}  "
                  f"unannotated {d['unannotated']/1e6:6.3f}")

