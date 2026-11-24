#!/usr/bin/env python3
"""Split a genome into ~TARGET-bp chunks for parallel RepeatMasker runs (benchmark only).

    split_genome.py <genome.fa> <outdir> [target_bp=50e6]

Sequences longer than 1.2 x target are cut into equal parts; each cut is moved to the middle of
the nearest run of >=100 Ns within 2 Mb, so that no repeat hit can span a cut where the assembly
has a gap. Shorter sequences are packed whole into chunks of up to the target size.
Writes <outdir>/chunks/cNNN.fa (records named pNNNNN), <outdir>/pieces.tsv
(piece, chunk, seq, seq_index, offset, length, seq_length; offset 0-based) and <outdir>/chunks.tsv.
"""
import os
import re
import sys

genome, out = sys.argv[1], sys.argv[2]
target = int(float(sys.argv[3])) if len(sys.argv) > 3 else 50_000_000
NRUN = re.compile(r'[Nn]{100,}')
WIN = 2_000_000
os.makedirs(os.path.join(out, 'chunks'), exist_ok=True)


def records(path):
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts)


def write_chunk(chunk, pieces):
    with open(os.path.join(out, 'chunks', chunk + '.fa'), 'w') as fh:
        for pid, s in pieces:
            fh.write('>' + pid + '\n')
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + '\n')


piece_rows, chunk_rows, small = [], [], []
n_piece = n_chunk = n_snapped = n_cuts = 0
total = 0


def new_piece():
    global n_piece
    n_piece += 1
    return 'p%05d' % n_piece


def new_chunk():
    global n_chunk
    n_chunk += 1
    return 'c%03d' % n_chunk


for idx, (name, seq) in enumerate(records(genome)):
    L = len(seq)
    total += L
    if L <= 1.2 * target:
        small.append((idx, name, seq))
        continue
    n = max(2, round(L / target))
    runs = [((m.start() + m.end()) // 2) for m in NRUN.finditer(seq)]
    cuts = []
    for i in range(1, n):
        t = round(i * L / n)
        near = [r for r in runs if abs(r - t) <= WIN and (not cuts or r > cuts[-1])]
        if near:
            cuts.append(min(near, key=lambda r: abs(r - t)))
            n_snapped += 1
        else:
            cuts.append(t)
        n_cuts += 1
    bounds = [0] + cuts + [L]
    for a, b in zip(bounds, bounds[1:]):
        pid, c = new_piece(), new_chunk()
        write_chunk(c, [(pid, seq[a:b])])
        piece_rows.append((pid, c, name, idx, a, b - a, L))
        chunk_rows.append((c, b - a, 1))
    del seq

small.sort(key=lambda x: -len(x[2]))
bins = []                      # [bp, [(idx, name, seq)]]
for rec in small:
    for b in bins:
        if b[0] + len(rec[2]) <= target:
            b[0] += len(rec[2])
            b[1].append(rec)
            break
    else:
        bins.append([len(rec[2]), [rec]])
for bp, recs in bins:
    c = new_chunk()
    pieces = []
    for idx, name, seq in recs:
        pid = new_piece()
        pieces.append((pid, seq))
        piece_rows.append((pid, c, name, idx, 0, len(seq), len(seq)))
    write_chunk(c, pieces)
    chunk_rows.append((c, bp, len(recs)))

# every base of every sequence must be in exactly one piece
cover = {}
for pid, c, name, idx, off, ln, sl in piece_rows:
    cover.setdefault(name, []).append((off, ln, sl))
for name, ps in cover.items():
    ps.sort()
    pos = 0
    for off, ln, sl in ps:
        assert off == pos, (name, off, pos)
        pos += ln
    assert pos == ps[0][2], (name, pos, ps[0][2])
assert sum(r[5] for r in piece_rows) == total

with open(os.path.join(out, 'pieces.tsv'), 'w') as fh:
    fh.write('piece\tchunk\tseq\tseq_index\toffset\tlength\tseq_length\n')
    for r in piece_rows:
        fh.write('\t'.join(map(str, r)) + '\n')
with open(os.path.join(out, 'chunks.tsv'), 'w') as fh:
    fh.write('chunk\tbp\tn_pieces\n')
    for r in chunk_rows:
        fh.write('\t'.join(map(str, r)) + '\n')
sizes = [r[1] for r in chunk_rows]
print(f'{len(cover)} sequences, {total:,} bp -> {len(chunk_rows)} chunks ({min(sizes):,}-{max(sizes):,} bp), '
      f'{len(piece_rows)} pieces; {n_snapped}/{n_cuts} cuts placed in N runs')
