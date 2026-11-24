#!/usr/bin/env python3
"""E00 - where do the development runs lose reference TE families and bases? (diagnostic)

Optimistic, reference-informed diagnosis of the frozen v10 library stages on the six
development genomes (route combination of the discovery-tool ablation). Nothing here is a
deployable rule: the oracle arms use the curated reference to pick candidates and labels.

Stages of one run (ablation/reclass/<sp>/<route>/out_{wo,with}):
  tool_<t>     raw library of each discovery tool         ablation/libs/<sp>/<t>.fa
  raw_union    all tool libraries together (= the dedup input)
  dedup        cluster representatives                    03_library/<sp>.dedup_library.fa
  refined      after DEDUP_REFINE                         03_library/<sp>.refined_library.fa
  gate_removed candidates the evidence gate removed       03_library/<sp>.gate_removed.fa
  gated        after the evidence gate                    03_library/<sp>.gated_library.fa
  final_wo     classified library w/o TEtrimmer           (gated sequences, new labels)
  curated      TEtrimmer output                           04_curation/<sp>.curated_library.fa
  final_with   classified library with TEtrimmer          (curated + rescued, re-deduplicated)

Family level uses the Table 2 rules (family_level.py: each stage library masked with the
curated reference). Base level decomposes every reference order's bases into TP / covered
with another label / missed, and every test order's bases into TP / other reference order
/ no reference TE (unsupported) / reference Unknown (don't care), for
  final_<arm>   the existing Table 2 genome maskings
  raw_union, refined          all candidates before / after dedup (tool labels)
  oracle_<arm>  final library + the earlier-stage candidates that recover reference families
                the final library lost (labelled with the reference family's class)
each with its own labels and with oracle labels (dominant reference family of the consensus).

usage: e00.py collect <sp> | families <sp> | baselevel <sp> | report
"""
import csv
import functools
import os
import statistics
import sys
from collections import defaultdict

import numpy as np

B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E, A, X = f"{B}/edta_manual_bench", f"{B}/ablation", f"{B}/e00"
sys.path.insert(0, E)
from family_level import classify, merged_bp    # noqa: E402  same rules as Table 2
from libtest_harmonised import canon, sizes    # noqa: E402  same labels as Table 2

ROUTE = {'arabidopsis_thaliana': 'edta+heliano+hite', 'oryza_sativa': 'edta+heliano+hite',
         'drosophila_melanogaster': 'earlgrey+hite', 'caenorhabditis_elegans': 'earlgrey+hite',
         'zymoseptoria_tritici': 'earlgrey+hite', 'septoria_linicola': 'earlgrey+hite'}
PATHS = {'wo': ['raw_union', 'dedup', 'refined', 'gated', 'final_wo'],
         'with': ['raw_union', 'dedup', 'refined', 'gated', 'curated', 'final_with']}
STEP = {'dedup': 'dedup', 'refined': 'refine', 'gated': 'gate', 'final_wo': 'classify',
        'curated': 'TEtrimmer', 'final_with': 'rescue+redup'}
ORDERS = ['LTR', 'LINE', 'SINE', 'DNA', 'RC']
LEVELS = [(3, 'perfect'), (2, 'good+'), (1, 'present+')]
GENOME_LIBS = ['final_wo', 'final_with', 'raw_union', 'refined', 'oracle_wo', 'oracle_with']


@functools.lru_cache(maxsize=None)
def order(label):
    """Harmonised order: LTR, LINE, SINE, DNA (TIR + MITE), RC, Unknown or nonTE."""
    c = canon(label or 'Unknown')
    if c == 'Unspecified':
        return 'nonTE'
    return {'LTR': 'LTR', 'LINE': 'LINE', 'SINE': 'SINE', 'TIR': 'DNA', 'MITE': 'DNA',
            'RC': 'RC'}.get(c.split('/', 1)[0], 'Unknown')


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(seq)
            name, seq = line[1:].rstrip('\n'), []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def split_header(h):
    ident, _, label = h.split()[0].partition('#')
    return ident, label or 'Unknown'


def stage_sources(sp):
    cb = ROUTE[sp]
    R = f"{A}/reclass/{sp}/{cb}"
    L = f"{R}/out_wo/{sp}/03_library/{sp}"
    s = {f"tool_{t}": f"{A}/libs/{sp}/{t}.fa" for t in cb.split('+')}
    s.update(dedup=f"{L}.dedup_library.fa", refined=f"{L}.refined_library.fa",
             gate_removed=f"{L}.gate_removed.fa", gated=f"{L}.gated_library.fa",
             final_wo=f"{A}/combos/{sp}/{cb}/library_wo.fa",
             curated=f"{R}/out_with/{sp}/04_curation/{sp}.curated_library.fa",
             final_with=f"{A}/combos/{sp}/{cb}/library_with.fa")
    return s


def collect(sp):
    """Stage FASTAs renamed q1..qN (RepeatMasker-safe) and stages/map.tsv."""
    cb = ROUTE[sp]
    d = f"{X}/{sp}/stages"
    os.makedirs(d, exist_ok=True)
    recs = {st: list(read_fasta(p)) for st, p in stage_sources(sp).items()}
    recs['raw_union'] = [r for t in cb.split('+') for r in recs[f"tool_{t}"]]
    with open(f"{d}/map.tsv", 'w') as m:
        m.write('stage\tqid\tid\tlabel\tlength\n')
        for st, rs in recs.items():
            if not rs:
                continue
            with open(f"{d}/{st}.fa", 'w') as fa:
                for i, (h, seq) in enumerate(rs, 1):
                    ident, label = split_header(h)
                    m.write(f"{st}\tq{i}\t{ident}\t{label}\t{len(seq)}\n")
                    fa.write(f">q{i}\n{seq}\n")
    # genome-masking libraries of the two pre-gate stages (tool labels)
    g = f"{X}/{sp}/genome_libs"
    os.makedirs(g, exist_ok=True)
    for st in ('raw_union', 'refined'):
        with open(f"{g}/{st}.fa", 'w') as fa:
            for i, (h, seq) in enumerate(recs[st], 1):
                fa.write(f">q{i}#{split_header(h)[1]}\n{seq}\n")
    print(sp, ' '.join(f"{st}={len(rs)}" for st, rs in recs.items()))


def read_out(path):
    """family_level.read_out plus query coordinates (query = library seq, rep = reference)."""
    hits = []
    for line in open(path):
        f = line.split()
        if len(f) < 14 or not f[0].isdigit():
            continue
        f = [x.strip('()') for x in f]
        q, qb, qe, qleft, strand, rep = f[4], int(f[5]), int(f[6]), int(f[7]), f[8], f[9]
        if strand == '+':
            rb, re_, rleft = int(f[11]), int(f[12]), int(f[13])
        else:
            rleft, re_, rb = int(f[11]), int(f[12]), int(f[13])
        hits.append({'q': q, 'qb': qb, 'qe': qe, 'qlen': qe + qleft, 'qaln': qe - qb,
                     'rep': rep, 'rlen': re_ + rleft, 'rb': rb, 're': re_, 'raln': re_ - rb,
                     'div': float(f[1])})
    return hits


def statuses(hits):
    perfect, good, present = classify(hits)
    st = dict.fromkeys(present, 1)
    st.update(dict.fromkeys(good, 2))
    st.update(dict.fromkeys(perfect, 3))
    return st


def pair_cover(hits, max_div=20.0):
    """(q, rep) -> (reference coverage, query coverage) of the hits below max_div."""
    pairs = defaultdict(list)
    for h in hits:
        if h['div'] < max_div:
            pairs[(h['q'], h['rep'])].append(h)
    return {k: (merged_bp([(h['rb'], h['re']) for h in hs]) / max(1, hs[0]['rlen']),
                merged_bp([(h['qb'], h['qe']) for h in hs]) / max(1, hs[0]['qlen']), hs)
            for k, hs in pairs.items()}


def best_candidates(cover):
    """Reference family -> (q, refcov, qcov) of the candidate covering most of it."""
    best = {}
    for (q, rep), (rc, qc, _) in cover.items():
        if rep not in best or (rc, qc) > best[rep][1:]:
            best[rep] = (q, rc, qc)
    return best


def dominant(cover):
    """Library seq -> reference family covering >= half of it (most bases)."""
    dom = {}
    for (q, rep), (_, qc, _) in cover.items():
        if qc >= 0.5 and (q not in dom or qc > dom[q][1]):
            dom[q] = (rep, qc)
    return {q: rep for q, (rep, _) in dom.items()}


def boundary(rc, qc):
    if rc > 0.95:
        return 'complete' if qc > 0.95 else 'overextended'
    if rc > 0.5:
        return 'truncated' if qc > 0.95 else 'partial_extra'
    return 'fragment'


def ref_info(sp):
    label, rlen = {}, {}
    for h, seq in read_fasta(f"{E}/{sp}/reference/library.fa"):
        n, lab = split_header(h)
        label[n], rlen[n] = lab, len(seq)
    return label, rlen


def ref_strata(sp, rlen):
    """Reference family -> (copies, full-length copies >= 80%, median divergence) in the genome."""
    spans, divs = defaultdict(list), defaultdict(list)
    for i, line in enumerate(open(f"{E}/{sp}/reference/{sp}.fna.out")):
        f = line.split()
        if i < 3 or len(f) < 15 or not f[0].isdigit():
            continue
        ff = [x.strip('()') for x in f]
        if f[8] == '+':
            rb, re_ = int(ff[11]), int(ff[12])
        else:
            re_, rb = int(ff[12]), int(ff[13])
        spans[(f[9], f[14])].append((rb, re_))
        divs[f[9]].append(float(f[1]))
    copies, full = defaultdict(int), defaultdict(int)
    for (rep, _), sp_ in spans.items():
        copies[rep] += 1
        if rlen.get(rep) and merged_bp(sp_) >= 0.8 * rlen[rep]:
            full[rep] += 1
    return {r: (copies[r], full[r], statistics.median(divs[r]) if divs[r] else None) for r in rlen}


def bins(c, fl, dv):
    return {'copies': '0' if c == 0 else '1-5' if c <= 5 else '6-20' if c <= 20 else '>20',
            'full_length': '0' if fl == 0 else '1-2' if fl <= 2 else '>=3',
            'divergence': 'na' if dv is None else '<5' if dv < 5 else '5-15' if dv < 15 else '>=15'}


def families(sp):
    d = f"{X}/{sp}"
    rlabel, rlen = ref_info(sp)
    rorder = {n: order(l) for n, l in rlabel.items()}
    fams = [n for n in rlabel if rorder[n] in ORDERS]
    qmap, qseq = defaultdict(dict), defaultdict(dict)
    for r in csv.DictReader(open(f"{d}/stages/map.tsv"), delimiter='\t'):
        qmap[r['stage']][r['qid']] = (r['id'], r['label'], int(r['length']))
    stages = [s for s in qmap if os.path.exists(f"{d}/famlevel/{s}/library.fa.out")]
    lost = sorted(set(qmap) - set(stages))
    if lost:
        sys.exit(f"{sp}: no family-level masking for {lost}")
    hits = {s: read_out(f"{d}/famlevel/{s}/library.fa.out") for s in stages}
    st = {s: statuses(hits[s]) for s in stages}
    cover = {s: pair_cover(hits[s]) for s in stages}
    best = {s: best_candidates(cover[s]) for s in stages}
    strata = ref_strata(sp, rlen)
    tools = sorted(s for s in stages if s.startswith('tool_'))
    order_st = tools + [s for s in ['raw_union', 'dedup', 'refined', 'gate_removed', 'gated',
                                    'final_wo', 'curated', 'final_with'] if s in stages]
    seen = {arm: [s for s in order_st if s in tools or s in PATHS[arm] or s == 'gate_removed']
            for arm in PATHS}
    S = lambda s, f: st[s].get(f, 0)

    def first_drop(f, path, lvl):
        if S(path[0], f) < lvl or S(path[-1], f) >= lvl:
            return ''
        return next(STEP[s] for s in path[1:] if S(s, f) < lvl)

    rows, summ, strat, lab, bnd = [], defaultdict(int), defaultdict(int), defaultdict(int), defaultdict(int)
    for f in fams:
        c, fl, dv = strata[f]
        row = {'family': f, 'order': rorder[f], 'label': rlabel[f], 'length': rlen[f],
               'copies': c, 'full_length': fl, 'median_div': '' if dv is None else f"{dv:.1f}"}
        row.update({f"st_{s}": S(s, f) for s in order_st})
        for arm, path in PATHS.items():
            fin = path[-1]
            orc = max(S(s, f) for s in seen[arm])
            row[f"{arm}_oracle"] = orc
            q, rc, qc = best[fin].get(f, ('', 0.0, 0.0))
            qid, qlab = (qmap[fin][q][:2] if q else ('', ''))
            row.update({f"{arm}_loss_perfect": first_drop(f, path, 3),
                        f"{arm}_loss_present": first_drop(f, path, 1),
                        f"{arm}_best": qid, f"{arm}_best_label": qlab,
                        f"{arm}_label_ok": '' if not q else int(order(qlab) == rorder[f]),
                        f"{arm}_refcov": f"{rc:.2f}", f"{arm}_qcov": f"{qc:.2f}",
                        f"{arm}_boundary": boundary(rc, qc) if q else 'none'})
            for o in (rorder[f], 'all'):
                for lvl, name in LEVELS:
                    k = (arm, o, name)
                    summ[k + ('n_ref',)] += 1
                    summ[k + ('raw',)] += S('raw_union', f) >= lvl
                    summ[k + ('final',)] += S(fin, f) >= lvl
                    summ[k + ('oracle',)] += orc >= lvl
                    summ[k + ('gained',)] += S('raw_union', f) < lvl <= S(fin, f)
                    step = first_drop(f, path, lvl)
                    if step:
                        summ[k + (f"lost_{step}",)] += 1
                if S(fin, f) >= 2:
                    lab[(arm, o, 'n')] += 1
                    lab[(arm, o, 'ok')] += order(qlab) == rorder[f]
                    lab[(arm, o, 'unknown')] += order(qlab) in ('Unknown', 'nonTE')
                if S(fin, f) >= 1:
                    bnd[(arm, o, boundary(rc, qc))] += 1
            for kind, b in bins(c, fl, dv).items():
                k = (arm, kind, b)
                strat[k + ('n',)] += 1
                strat[k + ('raw_perfect',)] += S('raw_union', f) == 3
                strat[k + ('final_perfect',)] += S(fin, f) == 3
                strat[k + ('raw_present',)] += S('raw_union', f) >= 1
                strat[k + ('final_present',)] += S(fin, f) >= 1
        rows.append(row)
    out = f"{d}/results"
    os.makedirs(out, exist_ok=True)
    with open(f"{out}/families.tsv", 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    steps = sorted({k[3] for k in summ if k[3].startswith('lost_')})
    with open(f"{out}/family_summary.tsv", 'w') as fh:
        cols = ['n_ref', 'raw', 'final', 'oracle', 'gained'] + steps
        fh.write('genome\tarm\torder\tlevel\t' + '\t'.join(cols) + '\n')
        for arm in PATHS:
            for o in ORDERS + ['all']:
                for _, name in LEVELS:
                    if summ[(arm, o, name, 'n_ref')]:
                        fh.write(f"{sp}\t{arm}\t{o}\t{name}\t" +
                                 '\t'.join(str(summ[(arm, o, name, c)]) for c in cols) + '\n')
    with open(f"{out}/label_boundary.tsv", 'w') as fh:
        kinds = ['complete', 'overextended', 'truncated', 'partial_extra', 'fragment']
        fh.write('genome\tarm\torder\tgood+_families\tlabel_ok\tlabel_unknown\t' + '\t'.join(kinds) + '\n')
        for arm in PATHS:
            for o in ORDERS + ['all']:
                fh.write(f"{sp}\t{arm}\t{o}\t{lab[(arm, o, 'n')]}\t{lab[(arm, o, 'ok')]}\t"
                         f"{lab[(arm, o, 'unknown')]}\t" + '\t'.join(str(bnd[(arm, o, k)]) for k in kinds) + '\n')
    with open(f"{out}/strata.tsv", 'w') as fh:
        cols = ['n', 'raw_perfect', 'final_perfect', 'raw_present', 'final_present']
        fh.write('genome\tarm\tstratum\tbin\t' + '\t'.join(cols) + '\n')
        for k in sorted({k[:3] for k in strat}):
            fh.write(f"{sp}\t" + '\t'.join(k) + '\t' + '\t'.join(str(strat[k + (c,)]) for c in cols) + '\n')

    # oracle libraries and oracle labels for the genome maskings
    for s in stages:
        if s in ('raw_union', 'refined', 'dedup', 'gated', 'curated', 'final_wo', 'final_with'):
            qseq[s] = {h.split()[0]: seq for h, seq in read_fasta(f"{d}/stages/{s}.fa")}
    g = f"{d}/genome_libs"
    # RepeatMasker names the internal / terminal part of a full-length LTR consensus
    # <family>-int / <family>-LTR; such hits never count for a family's status (as in
    # family_level.py) but give the label of their family here
    def ref_label(rep):
        for suf in ('-int', '-LTR', '_LTR', '_I', '-I'):
            if rep not in rlabel and rep.endswith(suf):
                rep = rep[:-len(suf)]
        return rlabel.get(rep)
    unnamed = {h['rep'] for s in stages for h in hits[s] if h['rep'] not in rlabel}
    if unnamed:
        print(f"{sp}: {len(unnamed)} hit names not in the reference, e.g. {sorted(unnamed)[:5]}; "
              f"{sum(ref_label(r) is None for r in unnamed)} without a family")
    dom = {s: {q: ref_label(rep) for q, rep in dominant(cover[s]).items() if ref_label(rep)} for s in stages}
    for s in ('raw_union', 'refined'):
        with open(f"{g}/{s}.oracle_labels.tsv", 'w') as fh:
            fh.writelines(f"{q}\t{l}\n" for q, l in dom[s].items())
    for arm, path in PATHS.items():
        fin = path[-1]
        base = {qmap[fin][q][0]: l for q, l in dom[fin].items()}
        with open(f"{g}/final_{arm}.oracle_labels.tsv", 'w') as fh:
            fh.writelines(f"{n}\t{l}\n" for n, l in base.items())
        used, n_add, fixed = set(), 0, defaultdict(int)
        with open(f"{g}/oracle_{arm}.fa", 'w') as fa, open(f"{g}/oracle_{arm}.map.tsv", 'w') as mp, \
                open(f"{g}/oracle_{arm}.oracle_labels.tsv", 'w') as ol:
            mp.write('name\tfamily\tfamily_label\tfrom_stage\tqid\tsource_id\tfrom_status\tfinal_status\n')
            for h, seq in read_fasta(stage_sources(sp)[fin]):
                fa.write(f">{h}\n{seq}\n")
            ol.writelines(f"{n}\t{l}\n" for n, l in base.items())
            for f in fams:
                pool = [(S(s, f), i, s) for i, s in enumerate(path[:-1])]
                bst, _, s = max(pool)
                if bst <= S(fin, f):
                    continue
                if bst == 3:        # the single perfect candidate (family_level.classify rule)
                    ok = [h for h in hits[s] if h['rep'] == f and h['qlen'] > 0 and h['rlen'] > 0
                          and h['qaln'] / h['qlen'] > 0.95 and h['raln'] / h['rlen'] > 0.95 and h['div'] < 5.0]
                    chosen = [max(ok, key=lambda h: h['qaln'])['q']]
                else:               # greedy: fewest candidates reaching the good / present cover
                    mdiv, target = (5.0, 0.95) if bst == 2 else (20.0, 0.8)
                    cand = {q: [(h['rb'], h['re']) for h in hs if h['div'] < mdiv]
                            for (q, rep), (_, _, hs) in cover[s].items() if rep == f}
                    cand = {q: v for q, v in cand.items() if v}
                    chosen, spans = [], []
                    while cand and merged_bp(spans) / rlen[f] <= target:
                        q = max(cand, key=lambda q: merged_bp(spans + cand[q]))
                        if merged_bp(spans + cand[q]) == merged_bp(spans):
                            break
                        spans += cand.pop(q)
                        chosen.append(q)
                for q in chosen:
                    if (s, q) in used:
                        continue
                    used.add((s, q))
                    n_add += 1
                    name = f"e00r{n_add}"
                    fa.write(f">{name}#{rlabel[f]}\n{qseq[s][q]}\n")
                    ol.write(f"{name}\t{rlabel[f]}\n")
                    mp.write(f"{name}\t{f}\t{rlabel[f]}\t{s}\t{q}\t{qmap[s][q][0]}\t{bst}\t{S(fin, f)}\n")
                fixed[rorder[f]] += 1
        print(f"{sp} oracle_{arm}: +{n_add} candidates for {sum(fixed.values())} families "
              f"({', '.join(f'{o} {n}' for o, n in sorted(fixed.items()))})")
    for arm in PATHS:
        print(f"{sp} {arm}: perfect raw {summ[(arm, 'all', 'perfect', 'raw')]} -> final "
              f"{summ[(arm, 'all', 'perfect', 'final')]} (oracle {summ[(arm, 'all', 'perfect', 'oracle')]}), "
              f"present+ raw {summ[(arm, 'all', 'present+', 'raw')]} -> final {summ[(arm, 'all', 'present+', 'final')]}")


def load_hits(path, relabel=None):
    by = defaultdict(list)
    for i, line in enumerate(open(path)):
        f = line.split()
        if i < 3 or len(f) < 11:
            continue
        lab = relabel.get(f[9], f[10]) if relabel else f[10]
        by[f[4]].append((int(f[5]) - 1, int(f[6]), order(lab)))
    return by


def decompose(ref, test, size):
    n = defaultdict(int)
    for c, L in size.items():
        R = {o: np.zeros(L, bool) for o in ORDERS}
        Rdc = np.zeros(L, bool)
        for s, e, o in ref.get(c, ()):
            (R[o] if o in R else Rdc)[s:e] = True
        T = {o: np.zeros(L, bool) for o in ORDERS + ['other']}
        for s, e, o in test.get(c, ()):
            T[o if o in T else 'other'][s:e] = True
        Rin = np.logical_or.reduce([R[o] for o in ORDERS])
        DC = Rdc & ~Rin
        Tany = np.logical_or.reduce(list(T.values()))
        n['any', 'tp'] += int((Tany & Rin).sum())
        n['any', 'fn'] += int((Rin & ~Tany).sum())
        n['any', 'fp'] += int((Tany & ~Rin & ~DC).sum())
        for o in ORDERS:
            To, Ro = T[o], R[o]
            Toth = np.logical_or.reduce([T[p] for p in ORDERS if p != o])
            n[o, 'ref'] += int(Ro.sum())
            n[o, 'test'] += int(To.sum())
            n[o, 'tp'] += int((To & Ro).sum())
            n[o, 'fn_other'] += int((Ro & ~To & Toth).sum())
            n[o, 'fn_unknown'] += int((Ro & ~To & ~Toth & T['other']).sum())
            n[o, 'fn_missed'] += int((Ro & ~Tany).sum())
            n[o, 'fp_cross'] += int((To & ~Ro & Rin).sum())
            n[o, 'fp_unsupported'] += int((To & ~Rin & ~DC).sum())
            n[o, 'dontcare'] += int((To & DC).sum())
    return n


def baselevel(sp):
    d, cb = f"{X}/{sp}", ROUTE[sp]
    outs = {'final_wo': f"{E}/{sp}/abl_wo_{cb}/{sp}.fna.out", 'final_with': f"{E}/{sp}/abl_with_{cb}/{sp}.fna.out"}
    outs.update({l: f"{d}/genome_rm/{l}/{sp}.fna.out" for l in GENOME_LIBS[2:]})
    _, size = sizes(sp)
    ref = load_hits(f"{E}/{sp}/reference/{sp}.fna.out")
    cols = ['ref', 'test', 'tp', 'fn_other', 'fn_unknown', 'fn_missed', 'fp_cross', 'fp_unsupported', 'dontcare']
    os.makedirs(f"{d}/results", exist_ok=True)
    with open(f"{d}/results/baselevel.tsv", 'w') as fh:
        fh.write('genome\tlibrary\tlabels\torder\t' + '\t'.join(cols) + '\tprecision\trecall\tF1\tdetection\n')
        for lib in GENOME_LIBS:
            if not os.path.exists(outs[lib]):
                print(f"{sp}/{lib}: no genome masking ({outs[lib]})")
                continue
            ol = f"{d}/genome_libs/{lib}.oracle_labels.tsv"
            relabel = dict(l.rstrip('\n').split('\t') for l in open(ol)) if os.path.exists(ol) else {}
            for mode, rl in (('own', None), ('oracle', relabel)):
                n = decompose(ref, load_hits(outs[lib], rl), size)
                tp, fn, fp = n['any', 'tp'], n['any', 'fn'], n['any', 'fp']
                fh.write(f"{sp}\t{lib}\t{mode}\tany_dc\t\t\t{tp}\t\t\t{fn}\t\t{fp}\t\t"
                         f"{tp / max(1, tp + fp):.4f}\t{tp / max(1, tp + fn):.4f}\t{2 * tp / max(1, 2 * tp + fp + fn):.4f}\t\n")
                for o in ORDERS:
                    v = [n[o, c] for c in cols]
                    tp, fp, ref_bp = n[o, 'tp'], n[o, 'fp_cross'] + n[o, 'fp_unsupported'], n[o, 'ref']
                    fn = ref_bp - tp
                    fh.write(f"{sp}\t{lib}\t{mode}\t{o}\t" + '\t'.join(map(str, v)) +
                             f"\t{tp / max(1, tp + fp):.4f}\t{tp / max(1, ref_bp):.4f}\t"
                             f"{2 * tp / max(1, 2 * tp + fp + fn):.4f}\t"
                             f"{(ref_bp - n[o, 'fn_missed']) / max(1, ref_bp):.4f}\n")
                print(f"{sp}/{lib}/{mode}: any_dc F1 {2 * n['any', 'tp'] / max(1, 2 * n['any', 'tp'] + n['any', 'fp'] + n['any', 'fn']):.4f}", flush=True)


def report():
    """One Markdown report over all genomes that have results."""
    tsv = lambda p: list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
    L = ['# E00 headroom diagnosis (development genomes, frozen v10 library stages)', '',
         'Optimistic diagnostic: oracle columns use the curated reference to choose candidates or labels.', '']
    for sp in ROUTE:
        r = f"{X}/{sp}/results"
        fs, lb, bl = tsv(f"{r}/family_summary.tsv"), tsv(f"{r}/label_boundary.tsv"), tsv(f"{r}/baselevel.tsv")
        if not fs:
            continue
        L += [f"## {sp} ({ROUTE[sp]})", '']
        steps = [c for c in fs[0] if c.startswith('lost_')]
        for arm in PATHS:
            L += [f"Family level, {arm} TEtrimmer: raw union -> final (oracle = best over the stages); "
                  "losses = families first dropping below the level at that step", '',
                  '| order | n | level | raw | final | oracle | gained | ' + ' | '.join(s[5:] for s in steps) + ' |',
                  '|' + '---|' * (6 + len(steps) + 1)]
            for x in fs:
                if x['arm'] == arm and x['level'] in ('perfect', 'present+'):
                    L.append(f"| {x['order']} | {x['n_ref']} | {x['level']} | {x['raw']} | {x['final']} | "
                             f"{x['oracle']} | {x['gained']} | " + ' | '.join(x[s] for s in steps) + ' |')
            L.append('')
            L += ['| order | good+ fam | label ok | label unknown | complete | overext | truncated | partial+extra | fragment |',
                  '|---|---|---|---|---|---|---|---|---|']
            for x in lb:
                if x['arm'] == arm:
                    L.append(f"| {x['order']} | {x['good+_families']} | {x['label_ok']} | {x['label_unknown']} | "
                             f"{x['complete']} | {x['overextended']} | {x['truncated']} | {x['partial_extra']} | {x['fragment']} |")
            L.append('')
        if bl:
            get = {(x['library'], x['labels'], x['order']): x for x in bl}
            L += ['Base level F1 per order (don\'t-care scoring); detection = reference bases covered by any label', '',
                  '| order | final_wo | +oracle labels | oracle_wo | +labels | final_with | +oracle labels | oracle_with | +labels | det final_with | det raw_union | det refined |',
                  '|---|---|---|---|---|---|---|---|---|---|---|---|']
            for o in ['any_dc'] + ORDERS:
                c = lambda lib, m, k='F1': get.get((lib, m, o), {}).get(k) or '-'
                L.append(f"| {o} | {c('final_wo', 'own')} | {c('final_wo', 'oracle')} | {c('oracle_wo', 'own')} | "
                         f"{c('oracle_wo', 'oracle')} | {c('final_with', 'own')} | {c('final_with', 'oracle')} | "
                         f"{c('oracle_with', 'own')} | {c('oracle_with', 'oracle')} | {c('final_with', 'own', 'detection')} | "
                         f"{c('raw_union', 'own', 'detection')} | {c('refined', 'own', 'detection')} |")
            L += ['', 'final_with base errors per order (Mb): FN covered by another order / Unknown label / missed; '
                  'FP on another reference order / on no reference TE', '',
                  '| order | ref Mb | FN other | FN unknown | FN missed | FP cross | FP unsupported |', '|---|---|---|---|---|---|---|']
            for o in ORDERS:
                x = get.get(('final_with', 'own', o))
                if x:
                    mb = lambda k: f"{int(x[k]) / 1e6:.2f}"
                    L.append(f"| {o} | {mb('ref')} | {mb('fn_other')} | {mb('fn_unknown')} | {mb('fn_missed')} | "
                             f"{mb('fp_cross')} | {mb('fp_unsupported')} |")
            L.append('')
    with open(f"{X}/results/E00_report.md", 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print(f"{X}/results/E00_report.md")


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'report':
        os.makedirs(f"{X}/results", exist_ok=True)
        report()
    else:
        {'collect': collect, 'families': families, 'baselevel': baselevel}[cmd](sys.argv[2])
