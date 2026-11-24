#!/usr/bin/env python3
"""thr_variants.py <workdir> <outdir> <species>...

Copy-threshold re-check for the evidence gate on the v10 development runs, without
re-running the pipeline. The thresholds only decide consensi with no domain, sequence
support or homology (tiers 3 and 4 of evidence.tsv): such a consensus is kept when it
has >= F full-length genomic loci, or >= H half-length ones (H = 0: off). For every
(F, H) on the grid, the w/o-TEtrimmer arm's final library is rebuilt from the default
one (F=2, H=5): records that would now fail are dropped, and records the default gate
removed (gate_removed.fa) that would now pass are added back with their gate labels.
Variants with the same membership as the default are not written.

Writes <outdir>/<species>/thr_f<F>_h<H>.fa and prints a membership table.
"""
import os, sys

GRID_F = (1, 2, 3)
GRID_H = (0, 5, 10)
DEFAULT = (2, 5)


def fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                recs.append((name, ''.join(buf)))
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def passes(row, f, h):
    nf, nh = int(row['full_copies'] or 0), int(row['half_copies'] or 0)
    return nf >= f or (h > 0 and nh >= h)


work, out = sys.argv[1], sys.argv[2]
for sp in sys.argv[3:]:
    base = f'{work}/{sp}/out_wo/{sp}'
    ev_path = f'{base}/03_library/{sp}.evidence.tsv'
    rows = [l.rstrip('\n').split('\t') for l in open(ev_path)]
    head = rows[0]
    ev = {r[0]: dict(zip(head, r)) for r in rows[1:]}
    final = fasta(f'{base}/05_classification/{sp}.te_library.fa')
    removed = fasta(f'{base}/03_library/{sp}.gate_removed.fa')
    copy_tier = {n for n, r in ev.items() if r['tier'] in ('3', '4')}
    default_keep = {n for n in copy_tier if passes(ev[n], *DEFAULT)}
    os.makedirs(f'{out}/{sp}', exist_ok=True)
    print(f'\n{sp}: {len(final)} final, {len(removed)} gate-removed, {len(copy_tier)} decided by copies '
          f'(default keeps {len(default_keep)})')
    for f in GRID_F:
        for h in GRID_H:
            keep = {n for n in copy_tier if passes(ev[n], f, h)}
            added = keep - default_keep
            dropped = default_keep - keep
            tag = f'thr_f{f}_h{h}'
            if not added and not dropped:
                print(f'  {tag}: same as default' + (' (default)' if (f, h) == DEFAULT else ''))
                continue
            lib = [(n, s) for n, s in final if n.split('#')[0] not in dropped]
            lib += [(n, s) for n, s in removed if n.split('#')[0] in added]
            bp_add = sum(len(s) for n, s in removed if n.split('#')[0] in added)
            bp_drop = sum(len(s) for n, s in final if n.split('#')[0] in dropped)
            with open(f'{out}/{sp}/{tag}.fa', 'w') as fh:
                for n, s in lib:
                    fh.write(f'>{n}\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
            print(f'  {tag}: {len(lib)} seqs  (+{len(added)} consensi / {bp_add} bp, '
                  f'-{len(dropped)} consensi / {bp_drop} bp)')

