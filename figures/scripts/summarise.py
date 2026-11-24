"""Print the derived numbers quoted in the manuscript text (from ./data, written by extract_data.py)."""
import collections
import figstyle as S

acc = S.read('accuracy.csv')
by = {(r['species'], r['method']): r for r in acc}
gen = {r['species']: r for r in S.read('genomes.csv')}
comp = collections.defaultdict(dict)
for r in S.read('compute.csv'):
    comp[r['species']][r['stage']] = (S.num(r['cpu_h']), S.num(r['wall_h']))

print('== accuracy: best F1 / most perfect per genome (TEAnnotate rows: teannotate, teannotate_wo)')
for sp in S.ORDER:
    rs = [r for r in acc if r['species'] == sp]
    if not rs:
        print(f'{sp:28s} pending')
        continue
    bf = max(rs, key=lambda r: S.num(r['f1']))
    bp = max(rs, key=lambda r: int(r['perfect'] or 0))
    t, w = by[(sp, 'teannotate')], by[(sp, 'teannotate_wo')]
    print(f"{sp:28s} best F1 {bf['method']:14s} {S.num(bf['f1']):.3f} | most perfect {bp['method']:14s} {bp['perfect']:>4s}"
          f" | TEA {S.num(t['f1']):.3f}/{t['perfect']} (lib {t['n_library']}) wo {S.num(w['f1']):.3f}/{w['perfect']} (lib {w['n_library']})"
          f" | dF1 {S.num(t['f1']) - S.num(w['f1']):+.3f} dlib {100 * (1 - int(t['n_library']) / int(w['n_library'])):.0f}%")

print('\n== compute (measured CPU-h, wall h)')
for sp in sorted(comp, key=lambda s: int(gen[s]['genome_bp'])):
    st = comp[sp]
    if 'tetrimmer' not in st:
        continue
    tot = sum(v[0] for v in st.values())
    disc = sum(v[0] for k, v in st.items() if k not in ('tetrimmer', 'library_other'))
    tt, lo = st['tetrimmer'], st['library_other']
    n_in = int(by[(sp, 'teannotate_wo')]['n_library'])
    print(f"{sp:28s} {int(gen[sp]['genome_bp']) / 1e9:5.2f} Gb total {tot:7.0f} | discovery {disc:7.0f} ({100 * disc / tot:3.0f}%)"
          f" | TEtrimmer {tt[0]:6.0f} CPU-h {tt[1]:5.1f} h ({100 * tt[0] / tot:3.0f}%) {n_in} inputs {n_in / tt[1]:6.1f}/h"
          f" | other lib {lo[0]:5.1f} CPU-h {lo[1]:4.1f} h")
    for k in ('earlgrey', 'hite', 'edta', 'heliano'):
        if k in st:
            print(f"{'':30s}{k:9s} {st[k][0]:7.1f} CPU-h {st[k][1]:5.1f} h")

print('\n== gate and curation on held-out genomes (v10 library stages)')
for r in S.read('gate_curation.csv'):
    n = int(r['n_dedup'])
    rm = int(r['action_remove'] or 0)
    nonte = int(r['action_to_nonte'] or 0)
    tiers = {k: int(r[f'tier_{k}'] or 0) for k in '01234'}
    kept = n - rm - nonte
    print(f"{r['species']:28s} dedup {n:5d} removed {rm:4d} ({100 * rm / n:4.1f}%, {int(r['removed_bp']) / 1e3:7.0f} kb) non-TE {nonte:3d}"
          f" kept {kept:5d} | domain {tiers['1']} support/homology {tiers['2']} copies {tiers['3']}"
          f" | TEtrimmer keep {r.get('tt_keep')} reinjected {r.get('tt_reinject')} relabel {r.get('tt_relabel_from_input')}")

print('\n== robustness: largest F1 shortfall from the best method of each genome')
for role in ('dev', 'heldout'):
    gaps = collections.defaultdict(list)
    for sp in S.ORDER:
        rs = [r for r in acc if r['species'] == sp and r['role'] == role]
        if not rs:
            continue
        best = max(S.num(r['f1']) for r in rs)
        f = {r['method']: S.num(r['f1']) for r in rs}
        f['teannotate_best_arm'] = max(f['teannotate'], f['teannotate_wo'])
        for m, v in f.items():
            gaps[m].append((best - v, sp))
    for m, g in gaps.items():
        w = max(g)
        print(f'{role:8s} {m:20s} worst shortfall {w[0]:.3f} ({w[1]}), n={len(g)}')
