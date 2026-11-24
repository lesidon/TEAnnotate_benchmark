#!/usr/bin/env python3
"""Merge the per-chunk RepeatMasker .out files of one method into one genome-coordinate .out
(benchmark only).

    merge_out.py <method>      -> merged/<method>/zea_mays.fna.out

Query names go back to the genome sequence names, begin/end are shifted by the piece offset,
(left) is recomputed against the full sequence length, and hit IDs are offset per chunk so they
stay unique. Lines are sorted by genome sequence order and position. Fields are joined by single
spaces (lib-test and the harmonised scorer split on whitespace).
"""
import csv
import os
import sys

M = os.path.dirname(os.path.abspath(__file__))
m = sys.argv[1]
pieces = {r['piece']: r for r in csv.DictReader(open(os.path.join(M, 'pieces.tsv')), delimiter='\t')}
chunk_ids = sorted({r['chunk'] for r in pieces.values()})
header, rows, id_off, n_lines = None, [], 0, 0
for c in chunk_ids:
    f = os.path.join(M, 'runs', m, c, c + '.fa.out')
    if not os.path.exists(f):
        sys.exit('missing ' + f)
    lines = open(f).read().splitlines()
    data = [l for l in lines if l.split() and l.split()[0].isdigit()]
    if not data:
        if not any('no repetitive sequences' in l for l in lines):
            sys.exit('no hits and no "no repetitive sequences" note in ' + f)
        continue
    if header is None:
        header = lines[:3]
    max_id = 0
    for l in data:
        f_ = l.split()
        pc = pieces[f_[4]]
        assert pc['chunk'] == c, (f_[4], c)
        off, sl = int(pc['offset']), int(pc['seq_length'])
        b, e = int(f_[5]) + off, int(f_[6]) + off
        assert 1 <= b <= e <= sl, (l, pc)
        hid = int(f_[14])
        max_id = max(max_id, hid)
        f_[4], f_[5], f_[6], f_[7] = pc['seq'], str(b), str(e), '(%d)' % (sl - e)
        f_[14] = str(hid + id_off)
        rows.append((int(pc['seq_index']), b, e, ' '.join(f_)))
    id_off += max_id
    n_lines += len(data)
rows.sort()
os.makedirs(os.path.join(M, 'merged', m), exist_ok=True)
dst = os.path.join(M, 'merged', m, 'zea_mays.fna.out')
with open(dst + '.tmp', 'w') as fh:
    fh.write('\n'.join(header or ['', '', '']) + '\n')
    for *_, l in rows:
        fh.write(l + '\n')
os.replace(dst + '.tmp', dst)
print('%s: %d hits from %d chunks -> %s' % (m, n_lines, len(chunk_ids), dst))
