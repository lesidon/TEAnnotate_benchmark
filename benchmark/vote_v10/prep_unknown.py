#!/usr/bin/env python3
"""prep_unknown.py <sp>... - the consensi each v9 run handed to the ensemble (Unknown at the
vote: the v9 voting-table rows plus every library entry still labelled Unknown), pooled per
organism_type with a '<code><arm>_' id prefix -> vote_v10/unknown/<org>.fasta."""
import os, sys
from vote_v10_common import ARMS, GENOMES, V, paths, read_fasta

pools = {}
for sp in sys.argv[1:]:
    for arm in ARMS:
        pre, org, lib, vt, _ = paths(sp, arm)
        if not (os.path.exists(lib) and os.path.exists(vt)):
            print(f"SKIP {sp}/{arm}: missing {lib if not os.path.exists(lib) else vt}")
            continue
        voted = {l.split("\t")[0] for i, l in enumerate(open(vt)) if i}
        n = 0
        for name, seq in read_fasta(lib):
            base, _, lab = name.split()[0].partition("#")
            if base in voted or lab.lower().startswith("unknown"):
                pools.setdefault(org, []).append((f"{pre}_{base}#Unknown", seq))
                n += 1
        print(f"{sp}/{arm} ({org}): {n} consensi Unknown at the vote")
TAG = os.environ.get("TAG", "")      # e.g. _ho: keeps a second batch apart
os.makedirs(f"{V}/unknown", exist_ok=True)
for org, recs in pools.items():
    with open(f"{V}/unknown/{org}{TAG}.fasta", "w") as fh:
        for h, seq in recs:
            fh.write(">" + h + "\n" + "\n".join(seq) + "\n")
    print(f"unknown/{org}{TAG}.fasta: {len(recs)} sequences")

