"""Collect every number used by the figures and tables of the TEAnnotate manuscript.

Reads the benchmark backup (read-only) and writes plain CSV tables to ./data/.
Standard library only.

    python extract_data.py [BENCHMARK_ROOT]

BENCHMARK_ROOT defaults to D:/teannotate_publication/whole_pipeline_benchmark
(the publication-core copy of $B on Zhores).

Sources
  accuracy, development genomes : ablation/results/ablation_{with,wo}{,_newdev}.tsv
                                  (+ sens/prec from edta_manual_bench/results/*_bench_harmonised_abl.csv)
  accuracy, held-out genomes    : edta_manual_bench/results/<Sp>_bench_harmonised.csv  (TE type any_dontcare)
                                  edta_manual_bench/results/<Sp>_family_level.csv      (class "all")
  classifier benchmark (Fig. 3) : ensemble_v3_dfam37/test/<kingdom>.wf1.two_step.tsv
                                  (ensemble_v2 consensi and split with the leakage filter extended to all
                                  Dfam 3.7 families; support-weighted F1 objective. The ensemble_v2 files
                                  are superseded and not used here)
  compute (Fig. 4)              : development - ablation/libs/<sp>/<tool>.trace.tsv (standalone discovery runs),
                                  ablation/combos/<sp>/<route>/trace_with.tsv (library stages)
                                  held-out    - work/heldout_v9{,e}/<sp>/out_with/pipeline_info (discovery),
                                  work/heldout_v10/<sp>/out_with/pipeline_info (library stages; heldout_v13 for maize)
  genome size                   : 01_genome_stats/<sp>.seqkit_stats.tsv of any pipeline run
  gate / curation counts        : work/heldout_v10/<sp>/out_with/<sp>/03_library/<sp>.evidence.tsv,
                                  .../04_curation/<sp>.tetrimmer_postprocess.tsv
"""
import csv, glob, os, re, sys, collections

ROOT = sys.argv[1] if len(sys.argv) > 1 else r'D:\teannotate_publication\whole_pipeline_benchmark'
RES = os.path.join(ROOT, 'edta_manual_bench', 'results')
ABL = os.path.join(ROOT, 'ablation')
WORK = os.path.join(ROOT, 'work')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# species id, CSV prefix, display name, role, route tools, kingdom (route)
SPECIES = [
    ('arabidopsis_thaliana', 'Arabidopsis_thaliana', 'A. thaliana', 'dev', 'edta+heliano+hite', 'plant'),
    ('oryza_sativa', 'Oryza_sativa', 'O. sativa', 'dev', 'edta+heliano+hite', 'plant'),
    ('drosophila_melanogaster', 'Drosophila_melanogaster', 'D. melanogaster', 'dev', 'earlgrey+hite', 'animal'),
    ('caenorhabditis_elegans', 'Caenorhabditis_elegans', 'C. elegans', 'dev', 'earlgrey+hite', 'animal'),
    ('zymoseptoria_tritici', 'Zymoseptoria_tritici', 'Z. tritici', 'dev', 'earlgrey+hite', 'fungus'),
    ('septoria_linicola', 'Septoria_linicola', 'S. linicola', 'dev', 'earlgrey+hite', 'fungus'),
    ('chlamydomonas_reinhardtii', 'Chlamydomonas_reinhardtii', 'C. reinhardtii', 'heldout', 'edta+heliano+hite', 'plant'),
    ('zea_mays', 'Zea_mays', 'Z. mays', 'heldout', 'edta+heliano+hite', 'plant'),
    ('caenorhabditis_briggsae', 'Caenorhabditis_briggsae', 'C. briggsae', 'heldout', 'earlgrey+hite', 'animal'),
    ('gallus_gallus', 'Gallus_gallus', 'G. gallus', 'heldout', 'earlgrey+hite', 'animal'),
    ('taeniopygia_guttata', 'Taeniopygia_guttata', 'T. guttata', 'heldout', 'earlgrey+hite', 'animal'),
    ('mus_musculus', 'Mus_musculus', 'M. musculus', 'heldout', 'earlgrey+hite', 'animal'),
    ('danio_rerio', 'Danio_rerio', 'D. rerio', 'heldout', 'earlgrey+hite', 'animal'),
    ('blumeria_hordei', 'Blumeria_hordei', 'B. hordei', 'heldout', 'earlgrey+hite', 'fungus'),
    ('ectocarpus_sp7', 'Ectocarpus_sp7', 'Ectocarpus sp. 7', 'heldout', 'earlgrey+hite', 'animal'),
]
TOOLS = ['earlgrey', 'edta', 'heliano', 'hite']
CLASSES = ['ltr', 'nonltr', 'tir', 'helitron']


def write(name, rows, fields):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'wrote data/{name} ({len(rows)} rows)')


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- Nextflow traces
_DUR = re.compile(r'([\d.]+)\s*(ms|d|h|m|s)')
_UNIT = {'d': 86400, 'h': 3600, 'm': 60, 's': 1, 'ms': 0.001}


def dur_s(text):
    """'1d 12h 59m 24s' / '28.4s' / '174ms' -> seconds (None for '-')."""
    if text in (None, '', '-'):
        return None
    return sum(float(v) * _UNIT[u] for v, u in _DUR.findall(text))


def proc_name(name):
    """'TEANNOTATE:CURATION:TETRIMMER (danio_rerio)' -> 'TETRIMMER'."""
    return name.split(' (')[0].split(':')[-1]


def read_trace(path):
    """Rows of a Nextflow trace with realtime_s and cpu_h, both trace styles (raw or human readable)."""
    rows = []
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rt, pc = r.get('realtime'), r.get('%cpu', '')
            if rt and re.fullmatch(r'\d+', rt):          # raw trace: ms and plain percent
                secs, pcpu = int(rt) / 1000, fnum(pc)
            else:
                secs, pcpu = dur_s(rt), fnum(pc.rstrip('%')) if pc not in ('', '-') else None
            r['proc'] = proc_name(r['name'])
            r['realtime_s'] = secs
            r['cpu_h'] = secs * pcpu / 100 / 3600 if secs is not None and pcpu is not None else None
            rows.append(r)
    return rows


def latest_tasks(trace_files):
    """Latest COMPLETED/CACHED row per process over several traces (sorted by their timestamped names)."""
    best = {}
    for f in sorted(trace_files):
        for r in read_trace(f):
            if r['status'] in ('COMPLETED', 'CACHED') and r['cpu_h'] is not None:
                best[r['proc']] = r
    return best


# ---------------------------------------------------------------- genome sizes
def genome_size(sp):
    pats = [os.path.join(WORK, '*', sp, 'out*', sp, '01_genome_stats', f'{sp}.seqkit_stats.tsv')]
    for p in pats:
        for f in sorted(glob.glob(p)):
            with open(f, newline='') as fh:
                r = next(csv.DictReader(fh, delimiter='\t'))
            return int(r['sum_len']), int(r['num_seqs'])
    return None, None


def count_fasta(path):
    if not os.path.exists(path):
        return None
    with open(path, errors='replace') as fh:
        return sum(1 for line in fh if line.startswith('>'))


# ---------------------------------------------------------------- accuracy
def harmonised(prefix, files):
    """{method: {te_type: row}} from one or more *_bench_harmonised*.csv files (later files win)."""
    out = collections.defaultdict(dict)
    for name in files:
        p = os.path.join(RES, f'{prefix}_{name}.csv')
        if not os.path.exists(p):
            continue
        with open(p, newline='') as fh:
            for r in csv.DictReader(fh):
                out[r['Method']][r['TE type']] = r
    return out


def family(prefix):
    p = os.path.join(RES, f'{prefix}_family_level.csv')
    out = {}
    if os.path.exists(p):
        with open(p, newline='') as fh:
            for r in csv.DictReader(fh):
                if r['class'] == 'all':
                    out[r['Method']] = r
    return out


def ablation_rows():
    rows = []
    for name in ('ablation_with.tsv', 'ablation_wo.tsv', 'ablation_with_newdev.tsv', 'ablation_wo_newdev.tsv'):
        p = os.path.join(ABL, 'results', name)
        with open(p, newline='') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                r['_file'] = name
                rows.append(r)
    # one row per (species, method); the _newdev files hold rice and C. elegans
    seen = {}
    for r in rows:
        seen[(r['species'], r['method'])] = r
    return seen


METHOD_LABEL = {'teannotate': 'TEAnnotate', 'teannotate_wo': 'TEAnnotate, no TEtrimmer', 'earlgrey': 'Earl Grey',
                'hite': 'HiTE', 'edta': 'EDTA', 'heliano': 'HELIANO'}


def accuracy():
    abl = ablation_rows()
    out = []
    for sp, prefix, disp, role, route, kingdom in SPECIES:
        fam = family(prefix)
        if role == 'dev':
            # Both ablation arms wrote to *_abl.csv, so the later arm overwrote the earlier one; on 12 Nov 2026
            # the lost arm was rescored into its own file: *_abl_with.csv for the first four dev genomes,
            # *_abl_wo.csv for C. elegans
            harm = harmonised(prefix, ['bench_harmonised_abl', 'bench_harmonised_abl_wo', 'bench_harmonised_abl_with'])
            src = {'teannotate': f'abl_with_{route}', 'teannotate_wo': f'abl_wo_{route}',
                   'earlgrey': 'fresh_earlgrey', 'hite': 'fresh_hite', 'edta': 'fresh_edta', 'heliano': 'fresh_heliano'}
            for m, meth in src.items():
                a = abl.get((sp, meth))
                if a is None:
                    continue
                h = harm.get(meth, {}).get('any_dontcare', {})
                fr = fam.get(meth, {})
                out.append(dict(species=sp, display=disp, role=role, kingdom=kingdom, method=m, label=METHOD_LABEL[m],
                                source=meth, f1=fnum(a['any_dc']), sens=fnum(h.get('sens')), prec=fnum(h.get('prec')),
                                f1_ltr=fnum(a['ltr']), f1_nonltr=fnum(a['nonltr']), f1_tir=fnum(a['tir']),
                                f1_helitron=fnum(a['helitron']), perfect=int(a['perfect']), good=int(a['good']),
                                n_library=int(a['n_lib']), n_reference=int(fr['n_reference']) if fr else None))
        else:
            harm = harmonised(prefix, ['bench_harmonised'])
            src = {'teannotate': ('teannotate_ref10', 'teannotate_ref10'),
                   'teannotate_wo': ('teannotate_wo_ref10', 'teannotate_wo_ref10'),
                   'earlgrey': ('earlgrey', 'earlgrey'), 'hite': ('hite', 'hite'), 'edta': ('edta', 'edta2'),
                   'heliano': ('heliano', 'heliano')}
            if sp == 'zea_mays':
                src['teannotate'] = ('teannotate_v13', 'teannotate_v13')
                src['teannotate_wo'] = ('teannotate_wo_v13', 'teannotate_wo_v13')
            for m, (hm, fm) in src.items():
                h = harm.get(hm)
                if not h:
                    continue
                fr = fam.get(fm, {})
                dc = h['any_dontcare']
                out.append(dict(species=sp, display=disp, role=role, kingdom=kingdom, method=m, label=METHOD_LABEL[m],
                                source=hm, f1=fnum(dc['F1']), sens=fnum(dc['sens']), prec=fnum(dc['prec']),
                                f1_ltr=fnum(h['ltr']['F1']), f1_nonltr=fnum(h['nonltr']['F1']),
                                f1_tir=fnum(h['tir']['F1']), f1_helitron=fnum(h['helitron']['F1']),
                                perfect=int(fr['perfect']) if fr else None, good=int(fr['good']) if fr else None,
                                n_library=int(fr['n_library']) if fr else None,
                                n_reference=int(fr['n_reference']) if fr else None))
    return out


ENS = os.path.join(ROOT, 'ensemble_v3_dfam37', 'test')
# rows of ensemble_v3_dfam37/test/<kingdom>.wf1.two_step.tsv (support-weighted objective, the headline run)
ENS_METHODS = [('TEclass2', 'TEclass2'), ('DeepTE', 'DeepTE'), ('CREATE', 'CREATE'), ('NeuralTE', 'NeuralTE'),
               ('Terrier', 'Terrier'), ('tuned, current voter', 'Weighted vote'),
               ('two-step tuned', 'Two-step vote')]


def classifiers():
    """Classifier-only benchmark on the leakage-filtered test halves (weighted-F1 objective files)."""
    out = []
    for k in ('plant', 'animal', 'fungi'):
        with open(os.path.join(ENS, f'{k}.wf1.two_step.tsv'), newline='') as fh:
            rows = {r['method']: r for r in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t')}
        for src, label in ENS_METHODS:
            r = rows[src]
            for level in ('order', 'sf'):
                lo, hi = (fnum(x) for x in r[f'{level}_wf1_ci'].split('-'))
                out.append(dict(kingdom=k, method=label, source=src, level=level, wf1=fnum(r[f'{level}_wf1']),
                                wf1_lo=lo, wf1_hi=hi, macro_f1=fnum(r[f'{level}_macro_f1']),
                                accuracy=fnum(r[f'{level}_accuracy']), coverage=fnum(r[f'{level}_coverage'])))
    return out


# ---------------------------------------------------------------- compute
def compute():
    out = []
    for sp, prefix, disp, role, route, kingdom in SPECIES:
        size, _ = genome_size(sp)
        gb = size / 1e9 if size else None
        tools = route.split('+')
        if role == 'dev':
            for t in tools:
                rows = read_trace(os.path.join(ABL, 'libs', sp, f'{t}.trace.tsv'))
                r = [x for x in rows if x['status'] == 'COMPLETED'][-1]
                out.append(dict(species=sp, display=disp, role=role, genome_gb=gb, stage=t, cpu_h=r['cpu_h'],
                                wall_h=r['realtime_s'] / 3600, source=f'ablation/libs/{sp}/{t}.trace.tsv'))
            lib = read_trace(os.path.join(ABL, 'combos', sp, route, 'trace_with.tsv'))
            src = f'ablation/combos/{sp}/{route}/trace_with.tsv'
        else:
            run = 'heldout_v9e' if kingdom == 'plant' else 'heldout_v9'
            disc = latest_tasks(glob.glob(os.path.join(WORK, run, sp, 'out_with', 'pipeline_info', 'execution_trace*.txt')))
            for t in tools:
                r = disc.get(t.upper())
                if r is None:
                    continue
                out.append(dict(species=sp, display=disp, role=role, genome_gb=gb, stage=t, cpu_h=r['cpu_h'],
                                wall_h=r['realtime_s'] / 3600, source=f'work/{run}/{sp}/out_with ({r["status"]})'))
            lib_run = 'heldout_v13' if sp == 'zea_mays' else 'heldout_v10'
            traces = sorted(glob.glob(os.path.join(WORK, lib_run, sp, 'out_with', 'pipeline_info', 'execution_trace*.txt')))
            if not traces:
                continue
            lib = list(latest_tasks(traces).values())
            src = f'work/{lib_run}/{sp}/out_with'
        tt = [r for r in lib if r['proc'] == 'TETRIMMER' and r['cpu_h'] is not None]
        other = [r for r in lib if r['proc'] != 'TETRIMMER' and r['cpu_h'] is not None]
        # The maize TETRIMMER trace records staging and a small retry, not the chunked run.
        if tt and sp != 'zea_mays':
            out.append(dict(species=sp, display=disp, role=role, genome_gb=gb, stage='tetrimmer',
                            cpu_h=tt[-1]['cpu_h'], wall_h=tt[-1]['realtime_s'] / 3600, source=src))
        out.append(dict(species=sp, display=disp, role=role, genome_gb=gb, stage='library_other',
                        cpu_h=sum(r['cpu_h'] for r in other), wall_h=sum(r['realtime_s'] for r in other) / 3600,
                        source=src))
    return out


# ---------------------------------------------------------------- genomes, gate, curation
def genomes(acc):
    ref = {(r['species']): r['n_reference'] for r in acc if r['method'] == 'teannotate' and r['n_reference']}
    out = []
    for sp, prefix, disp, role, route, kingdom in SPECIES:
        size, nseq = genome_size(sp)
        fam = family(prefix)
        nref = ref.get(sp) or next((int(v['n_reference']) for v in fam.values()), None)
        out.append(dict(species=sp, display=disp, role=role, route=kingdom, genome_bp=size, genome_mb=round(size / 1e6, 1) if size else None,
                        n_seqs=nseq, ref_families_famlevel=nref,
                        ref_records_fasta=count_fasta(os.path.join(ROOT, 'curated', sp, f'{sp}.curated.lib.fa'))))
    return out


def gate():
    out = []
    for sp, prefix, disp, role, route, kingdom in SPECIES:
        lib_run = 'heldout_v13' if sp == 'zea_mays' else 'heldout_v10'
        base = os.path.join(WORK, lib_run, sp, 'out_with', sp)
        ev = os.path.join(base, '03_library', f'{sp}.evidence.tsv')
        if not os.path.exists(ev):
            continue
        with open(ev, newline='') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        act = collections.Counter(r['action'] for r in rows)
        tier = collections.Counter(r['tier'] for r in rows)
        bp_rm = sum(int(r['length']) for r in rows if r['action'] == 'remove')
        pp = os.path.join(base, '04_curation', f'{sp}.tetrimmer_postprocess.tsv')
        ppa = collections.Counter()
        if os.path.exists(pp):
            with open(pp, newline='') as fh:
                ppa = collections.Counter(r['action'] for r in csv.DictReader(fh, delimiter='\t'))
        out.append(dict(species=sp, display=disp, n_dedup=len(rows),
                        **{f'action_{k}': v for k, v in sorted(act.items())},
                        **{f'tier_{k}': v for k, v in sorted(tier.items())}, removed_bp=bp_rm,
                        **{f'tt_{k}': v for k, v in sorted(ppa.items())}))
    return out


if __name__ == '__main__':
    acc = accuracy()
    write('accuracy.csv', acc, ['species', 'display', 'role', 'kingdom', 'method', 'label', 'source', 'f1', 'sens',
                                'prec', 'f1_ltr', 'f1_nonltr', 'f1_tir', 'f1_helitron', 'perfect', 'good',
                                'n_library', 'n_reference'])
    cls = classifiers()
    write('classifiers.csv', cls, ['kingdom', 'method', 'source', 'level', 'wf1', 'wf1_lo', 'wf1_hi', 'macro_f1',
                                   'accuracy', 'coverage'])
    comp = compute()
    write('compute.csv', comp, ['species', 'display', 'role', 'genome_gb', 'stage', 'cpu_h', 'wall_h', 'source'])
    gen = genomes(acc)
    write('genomes.csv', gen, ['species', 'display', 'role', 'route', 'genome_bp', 'genome_mb', 'n_seqs',
                               'ref_families_famlevel', 'ref_records_fasta'])
    g = gate()
    keys = sorted({k for r in g for k in r} - {'species', 'display'})
    write('gate_curation.csv', g, ['species', 'display'] + keys)
