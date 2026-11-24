#!/usr/bin/env python3
"""Family-level library evaluation after Flynn, Hubley et al. 2020 (RepeatModeler2, PNAS).

Input: for each method, the RepeatMasker .out of its library masked with the curated
reference library (RepeatMasker -lib reference.fa -nolow library.fa), so each hit maps a
library consensus (query) onto a reference family (repeat). Re-implements the paper's
get_family_summary_paper.sh:

  perfect    one library sequence covers >95% of the reference family, the reference
             covers >95% of that sequence (not over-extended), divergence <5%
  good       several library sequences with divergence <5% together cover >95% of the
             reference family (and it is not perfect)
  present    library sequences with divergence <20% together cover >80% of the
             reference family (includes perfect and good)
  not found  everything else

Coverage merges intervals closer than 10 bp (bedtools merge -d 10, as in the script).
Adds a per-class breakdown (reference family class group) and library-side counts: how
many library sequences hit any reference family at <20% divergence.

Usage: family_level.py <species> <species_dir> <out.csv>
  <species_dir>/reference/library.fa         curated reference (for the family list)
  <species_dir>/<method>/famlevel/library.fa.out
"""
import os
import re
import sys
from collections import defaultdict

GROUP = {
    'ltr': 'LTR', 'dirs': 'LTR', 'erv': 'LTR',
    'line': 'LINE', 'sine': 'SINE', 'ple': 'PLE', 'penelope': 'PLE',
    'dna': 'DNA', 'dnaauto': 'DNA', 'dnanona': 'DNA', 'mite': 'DNA', 'tir': 'DNA',
    'rc': 'RC', 'helitron': 'RC',
}


def group_of(label):
    top = label.split('/', 1)[0].strip().lower()
    if not top or top.startswith(('unknown', 'unspecified', 'unclassified')):
        return 'Unknown'
    return GROUP.get(top, 'Other')


def reference_families(path):
    fams = {}
    for line in open(path):
        if line.startswith('>'):
            tok = line[1:].split()[0]
            name, _, label = tok.partition('#')
            fams[name] = group_of(label)
    return fams


def library_names(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith('>')]


def read_out(path):
    """Hits as dicts; coordinates 1-based; query = library seq, rep = reference family."""
    hits = []
    for line in open(path):
        f = line.split()
        if len(f) < 14 or not f[0].isdigit():
            continue
        f = [x.strip('()') for x in f]
        div = float(f[1])
        q, qb, qe, qleft, strand, rep = f[4], int(f[5]), int(f[6]), int(f[7]), f[8], f[9]
        if strand == '+':
            rb, re_, rleft = int(f[11]), int(f[12]), int(f[13])
        else:                               # 'C': (left) end begin
            rleft, re_, rb = int(f[11]), int(f[12]), int(f[13])
        hits.append({
            'q': q, 'qlen': qe + qleft, 'qaln': qe - qb,
            'rep': rep, 'rlen': re_ + rleft, 'rb': rb, 're': re_, 'raln': re_ - rb,
            'div': div,
        })
    return hits


def merged_bp(spans, gap=10):
    total, cur = 0, None
    for s, e in sorted(spans):
        if cur and s <= cur[1] + gap:
            cur[1] = max(cur[1], e)
        else:
            if cur:
                total += cur[1] - cur[0]
            cur = [s, e]
    if cur:
        total += cur[1] - cur[0]
    return total


def classify(hits):
    perfect = set()
    for h in hits:
        if h['qlen'] <= 0 or h['rlen'] <= 0 or h['qaln'] / h['qlen'] <= 0.1:
            continue
        if h['raln'] / h['rlen'] > 0.95 and h['qaln'] / h['qlen'] > 0.95 and h['div'] < 5.0:
            perfect.add(h['rep'])
    by_rep = defaultdict(list)
    for h in hits:
        by_rep[h['rep']].append(h)
    present, good_or_perfect = set(), set()
    for rep, hs in by_rep.items():
        rlen = hs[0]['rlen']
        if rlen <= 0:
            continue
        if merged_bp([(h['rb'], h['re']) for h in hs if h['div'] < 20.0]) / rlen > 0.8:
            present.add(rep)
        if merged_bp([(h['rb'], h['re']) for h in hs if h['div'] < 5.0]) / rlen > 0.95:
            good_or_perfect.add(rep)
    good = good_or_perfect - perfect
    present |= perfect | good          # a perfect/good family is also present
    return perfect, good, present


def main():
    species, sdir, out_csv = sys.argv[1:4]
    ref = reference_families(os.path.join(sdir, 'reference', 'library.fa'))
    rows = []
    for method in sorted(os.listdir(sdir)):
        o = os.path.join(sdir, method, 'famlevel', 'library.fa.out')
        lib = os.path.join(sdir, method, 'library.fa')
        if method == 'reference' or not os.path.exists(o) or not os.path.exists(lib):
            continue
        hits = read_out(o)
        perfect, good, present = classify(hits)
        names = library_names(lib)
        matched = {h['q'] for h in hits if h['div'] < 20.0}
        classes = sorted(set(ref.values())) + ['all']
        for cls in classes:
            fams = {f for f, c in ref.items() if cls == 'all' or c == cls}
            p, g, pr = len(perfect & fams), len(good & fams), len(present & fams)
            rows.append([method, cls, len(fams), p, g, pr - p - g, len(fams) - pr,
                         len(names), len(matched)])
    with open(out_csv, 'w') as fh:
        fh.write('Method,class,n_reference,perfect,good,present_only,not_found,'
                 'n_library,n_library_matched\n')
        for r in rows:
            fh.write(','.join(str(x) for x in r) + '\n')
    print(out_csv, len(rows), 'rows')


if __name__ == '__main__':
    main()
