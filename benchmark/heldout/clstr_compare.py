#!/usr/bin/env python3
"""Compare two CD-HIT-format .clstr files over the same sequences: cluster counts and
pairwise co-clustering agreement (pairs together in A, in B, in both)."""
import re
import sys
from itertools import combinations

NAME = re.compile(r'>(.+?)\.\.\. ')


def read(path):
    cl, cur = {}, -1
    for line in open(path):
        if line.startswith('>Cluster'):
            cur += 1
            continue
        m = NAME.search(line)
        if m:
            cl[m.group(1).split()[0]] = cur
    return cl


def pairs(cl):
    by = {}
    for n, c in cl.items():
        by.setdefault(c, []).append(n)
    return {tuple(sorted(p)) for ms in by.values() for p in combinations(ms, 2)}


a, b = read(sys.argv[1]), read(sys.argv[2])
if set(a) != set(b):
    print(f'sequence sets differ: {len(set(a) ^ set(b))}')
pa, pb = pairs(a), pairs(b)
both = len(pa & pb)
print(f'clusters A(cd-hit) {len(set(a.values()))}  B(mmseqs) {len(set(b.values()))} | '
      f'co-clustered pairs A {len(pa)}  B {len(pb)}  both {both} | '
      f'B recovers {both / max(len(pa), 1):.3f} of A, A recovers {both / max(len(pb), 1):.3f} of B')
