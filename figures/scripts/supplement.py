"""Supplementary tables and figures of the TEAnnotate manuscript.

    python supplement.py [BENCHMARK_ROOT]

Writes LaTeX table bodies to ../../supp/ (input by supplementary.tex) and Figs. S1-S2 to ../ (figures/).
Reads data/accuracy.csv (extract_data.py), data/classifier_paired_bootstrap.tsv (paired_boot.py) and,
read-only, the benchmark backup:
  Table S3  edta_manual_bench/results/<Sp>_bench_harmonised_thr.csv, <Sp>_family_level.csv (thr_* rows)
  Table S4  ensemble_v3_dfam37/split/<kingdom>.split.tsv (ensemble_v2 split minus the Dfam 3.7 matches)
  Table S5  edta_manual_bench/pilot/pilot_results.csv, pilot/<sp>_family_level.csv and the
            'teannotate w/o tetrimmer (v4)' rows of results/<Sp>_bench_edta_manual.csv (no gate)
  Figs S1-2 per-order reference base pairs (TP+FN) from results/<Sp>_bench_harmonised*.csv
Table S2 is taken from bin/ensemble_vote.py TWO_STEP_DEFAULTS and nextflow.config of the pipeline.
"""
import csv, os, sys, collections
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import figstyle as S

ROOT = sys.argv[1] if len(sys.argv) > 1 else r'D:\teannotate_publication\whole_pipeline_benchmark'
RES = os.path.join(ROOT, 'edta_manual_bench', 'results')
PILOT = os.path.join(ROOT, 'edta_manual_bench', 'pilot')
SPLIT = os.path.join(ROOT, 'ensemble_v3_dfam37', 'split')
SUPP = os.path.join(os.path.dirname(S.FIGDIR), 'supp')
os.makedirs(SUPP, exist_ok=True)

PREFIX = {'arabidopsis_thaliana': 'Arabidopsis_thaliana', 'oryza_sativa': 'Oryza_sativa',
          'drosophila_melanogaster': 'Drosophila_melanogaster', 'caenorhabditis_elegans': 'Caenorhabditis_elegans',
          'zymoseptoria_tritici': 'Zymoseptoria_tritici', 'septoria_linicola': 'Septoria_linicola',
          'chlamydomonas_reinhardtii': 'Chlamydomonas_reinhardtii', 'zea_mays': 'Zea_mays',
          'caenorhabditis_briggsae': 'Caenorhabditis_briggsae', 'gallus_gallus': 'Gallus_gallus',
          'taeniopygia_guttata': 'Taeniopygia_guttata', 'mus_musculus': 'Mus_musculus',
          'danio_rerio': 'Danio_rerio', 'blumeria_hordei': 'Blumeria_hordei', 'ectocarpus_sp7': 'Ectocarpus_sp7'}
DEV4 = ['arabidopsis_thaliana', 'drosophila_melanogaster', 'zymoseptoria_tritici', 'septoria_linicola']
SHORT = {'arabidopsis_thaliana': r'\textit{A.~thaliana}', 'drosophila_melanogaster': r'\textit{D.~melanogaster}',
         'zymoseptoria_tritici': r'\textit{Z.~tritici}', 'septoria_linicola': r'\textit{S.~linicola}'}


def rows(path, delim=','):
    with open(path, newline='') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')), delimiter=delim))


def put(name, text):
    with open(os.path.join(SUPP, name), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(text)
    print('wrote supp/' + name)


def f3(x):
    return '%.3f' % float(x)


# ------------------------------------------------------------------ Table S2: vote settings
VOTE = {  # bin/ensemble_vote.py TWO_STEP_DEFAULTS (identical to ensemble_v2/tune/ and, after the extended
          # leakage filter, ensemble_v3_dfam37/tune/<k>.wf1.two_step.best.json)
    'Plants': ({'TEclass2': 1.5, 'CREATE': 1.0, 'NeuralTE': 2.0, 'Terrier': 2.0}, 0.3, 1,
               {'TEclass2': 1.0, 'CREATE': 1.0, 'NeuralTE': 3.0, 'Terrier': 1.0}, 0.5),
    'Animals': ({'TEclass2': 1.5, 'NeuralTE': 1.0, 'Terrier': 2.0}, 0.3, 1,
                {'NeuralTE': 1.0, 'Terrier': 2.0}, 0.5),
    'Fungi': ({'CREATE': 1.0, 'NeuralTE': 1.5, 'Terrier': 2.0}, 0.3, 1,
              {'NeuralTE': 2.0, 'Terrier': 1.0}, 0.0),
}
CLS = ['TEclass2', 'DeepTE', 'CREATE', 'NeuralTE', 'Terrier']
lines = []
for k, (ow, ofr, oag, sw, sfr) in VOTE.items():
    lines.append('%s, order & %s & %.1f & %d \\\\' % (k, ' & '.join('%g' % ow[c] if c in ow else '--' for c in CLS), ofr, oag))
    lines.append('%s, superfamily & %s & %s & -- \\\\' % (k, ' & '.join('%g' % sw[c] if c in sw else '--' for c in CLS),
                                                       'plurality' if sfr == 0 else '%.1f' % sfr))
put('tab_s2_votes.tex', '\n'.join(lines) + '\n')

# ------------------------------------------------------------------ Table S3: gate copy-number grid
VARIANTS = [('thr_f1_h0', 1, 0), ('thr_f1_h5', 1, 5), ('thr_f1_h10', 1, 10), ('thr_f2_h0', 2, 0),
            ('DEFAULT', 2, 5), ('thr_f2_h10', 2, 10), ('thr_f3_h0', 3, 0), ('thr_f3_h5', 3, 5),
            ('thr_f3_h10', 3, 10)]
# all six development genomes in the order of Table 2; the default (F2/H5) library of the first four is the
# v10 w/o-TEtrimmer run (thresholds/thr_bench.sh), of C. elegans and rice the ablation route library
# (thresholds/thr_bench_newdev.sh), i.e. the 'no TEtrimmer' library of Table 2
GRID_SP = ['arabidopsis_thaliana', 'oryza_sativa', 'drosophila_melanogaster', 'caenorhabditis_elegans',
           'zymoseptoria_tritici', 'septoria_linicola']
DEFAULT_M = {'caenorhabditis_elegans': 'abl_wo_earlgrey+hite', 'oryza_sativa': 'abl_wo_edta+heliano+hite'}
grid = {}
for sp in GRID_SP:
    p = PREFIX[sp]
    hp = os.path.join(RES, f'{p}_bench_harmonised_thr.csv')
    if not os.path.exists(hp):
        print('  S3: no grid results yet for', sp)
        continue
    h = {r['Method']: r for r in rows(hp) if r['TE type'] == 'any_dontcare'}
    fam = {r['Method']: r for r in rows(os.path.join(RES, f'{p}_family_level.csv')) if r['class'] == 'all'}
    for v, _, _ in VARIANTS:
        m = DEFAULT_M.get(sp, 'teannotate_wo_ref10') if v == 'DEFAULT' else v
        if m in h:
            grid[(sp, v)] = (float(h[m]['F1']), int(fam[m]['perfect']) if m in fam else None,
                             int(fam[m]['n_library']) if m in fam else None)
lines = []
for v, f, hcop in VARIANTS:
    cells = []
    for sp in GRID_SP:
        g = grid.get((sp, v))
        cells.append('--' if g is None else '%s / %s' % (f3(g[0]), g[1] if g[1] is not None else '--'))
    lab = '%d & %d%s' % (f, hcop, r' (default)' if v == 'DEFAULT' else '')
    lines.append(lab + ' & ' + ' & '.join(cells) + ' \\\\')
put('tab_s3_thresholds.tex', '\n'.join(lines) + '\n')
for (sp, m), v in sorted(grid.items()):
    print('  S3', sp, m, v)

# ------------------------------------------------------------------ Table S4: cluster split
lines = []
for k, disp in (('plant', 'Plants'), ('animal', 'Animals'), ('fungi', 'Fungi')):
    rs = rows(os.path.join(SPLIT, f'{k}.split.tsv'), '\t')
    by = collections.defaultdict(lambda: collections.Counter())
    cl = collections.defaultdict(set)
    for r in rs:
        by[r['source']][r['split']] += 1
        cl[(r['source'], r['split'])].add(r['cluster'])
    tot = collections.Counter(r['split'] for r in rs)
    tcl = {s: len({r['cluster'] for r in rs if r['split'] == s}) for s in ('train', 'test')}
    for src in sorted(by):
        lines.append('%s & %s & %s & %s & %s & %s \\\\' % (
            disp, src.replace('_v', ' v'), format(by[src]['train'], ','), format(by[src]['test'], ','),
            format(len(cl[(src, 'train')]), ','), format(len(cl[(src, 'test')]), ',')))
        disp = ''
    if len(by) > 1:
        lines.append(' & all & %s & %s & %s & %s \\\\' % (format(tot['train'], ','), format(tot['test'], ','),
                                                           format(tcl['train'], ','), format(tcl['test'], ',')))
    lines.append('\\addlinespace')
put('tab_s4_split.tex', '\n'.join(lines[:-1]).replace(',', '{,}') + '\n')

# ------------------------------------------------------------------ Table S5: evidence-gate pilot
pil = {(r['species'], r['variant'], r['TE type']): r for r in rows(os.path.join(PILOT, 'pilot_results.csv'))}
pfam = {}
for sp in DEV4:
    for r in rows(os.path.join(PILOT, f'{sp}_family_level.csv')):
        if r['class'] == 'all':
            pfam[(sp, r['Method'])] = r
base = {}
for sp in DEV4:
    for r in rows(os.path.join(RES, f'{PREFIX[sp]}_bench_edta_manual.csv')):
        if r['Method'] == 'teannotate w/o tetrimmer (v4)':
            base[(sp, r['TE type'])] = r
lines = []
for sp in DEV4:
    first = True
    for label, get, fam in (('No gate', lambda t: base.get((sp, t)), None),
                            ('Gate (adopted)', lambda t: pil.get((sp, 'evid_lite', t)), 'evid_lite'),
                            ('Strict gate', lambda t: pil.get((sp, 'evid_strict', t)), 'evid_strict')):
        a = get('any')
        cls = [get(t) for t in ('ltr', 'nonltr', 'tir', 'helitron')]
        fr = pfam.get((sp, fam)) if fam else None
        lines.append('%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\' % (
            SHORT[sp] if first else '', label, f3(a['sens']), f3(a['prec']), f3(a['F1']),
            *[f3(c['F1']) if c and int(c['TP']) + int(c['FN']) > 0 else '--' for c in cls],
            fr['perfect'] if fr else '--',
            fr['n_library'] if fr else '--'))
        first = False
    lines.append('\\addlinespace')
put('tab_s5_gate.tex', '\n'.join(lines[:-1]) + '\n')

# ------------------------------------------------------------------ Table S6: paired bootstrap
pb = rows(os.path.join(S.DATA, 'classifier_paired_bootstrap.tsv'), '\t')
lines = []
KD = {'plant': 'Plants', 'animal': 'Animals', 'fungi': 'Fungi'}
for r in pb:
    lo = float(r['ci_lo'])
    lo_s = ('%.4f' % abs(lo)) if abs(lo) < 0.001 else '%.3f' % abs(lo)
    lines.append('%s & %s & %s & %s & %s & %s & %s & %s & %s--%s \\\\' % (
        KD[r['kingdom']], r['level'], format(int(r['n_test']), ','), format(int(r['n_clusters']), ','),
        r['best_single'], f3(r['vote_wf1']), f3(r['single_wf1']), '%.3f' % float(r['delta']),
        ('$-$' if lo < 0 else '') + lo_s, '%.3f' % float(r['ci_hi'])))
put('tab_s6_paired.tex', '\n'.join(lines).replace(',', '{,}') + '\n')

# ------------------------------------------------------------------ Figs S1-S2: per-order F1
S.setup()
mpl.rcParams['hatch.linewidth'] = 0.55
acc = S.read('accuracy.csv')
by = {(r['species'], r['method']): r for r in acc}
disp = {r['species']: r['display'] for r in S.read('genomes.csv')}
TYPE = {'f1_ltr': 'ltr', 'f1_nonltr': 'nonltr', 'f1_tir': 'tir', 'f1_helitron': 'helitron'}


def ref_bp(sp):
    """{te type: reference bp (TP+FN)} from any harmonised scorecard of the genome (same truth for all methods)."""
    p = PREFIX[sp]
    for name in ('bench_harmonised', 'bench_harmonised_abl', 'bench_harmonised_abl_with'):
        f = os.path.join(RES, f'{p}_{name}.csv')
        if os.path.exists(f):
            out = {}
            for r in rows(f):
                if r['TE type'] in TYPE.values() and r['TE type'] not in out:
                    out[r['TE type']] = int(float(r['TP'])) + int(float(r['FN']))
            return out
    return {}


REF = {sp: ref_bp(sp) for sp in S.ORDER if sp not in S.PENDING}
for sp, d in REF.items():
    print('  ref bp', sp, d)
BAR_W, n_dev = 0.135, 6
x_of = {sp: i + (0.6 if i >= n_dev else 0) for i, sp in enumerate(S.ORDER)}


def style(m):
    col = S.METHODS[m][1]
    if m == 'teannotate_wo':
        return dict(facecolor='white', edgecolor=col, hatch='//////', linewidth=0.6)
    return dict(facecolor=col, edgecolor='white', linewidth=0.4)


def order_figure(stem, panels):
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 4.4), sharex=True, gridspec_kw={'hspace': 0.14})
    for ax, (col, ylabel) in zip(axes, panels):
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_ylabel(ylabel)
        ax.grid(axis='x', visible=False)
        for sp in S.ORDER:
            x = x_of[sp]
            if sp in S.PENDING or REF.get(sp, {}).get(TYPE[col], 0) < 1000:
                ax.axvspan(x - 0.44, x + 0.44, color=S.GRID, alpha=0.6, lw=0, zorder=0)
                ax.text(x, 0.5, 'not yet scored' if sp in S.PENDING else 'none in reference', rotation=90,
                        ha='center', va='center', fontsize=6, color=S.INK2)
                continue
            # HELIANO detects only Helitrons: its F1 for other orders says nothing about it
            present = [(m, S.num(by[(sp, m)][col])) for m in S.METHODS if (sp, m) in by
                       and (m != 'heliano' or col == 'f1_helitron')]
            present = [(m, v) for m, v in present if v is not None]
            x0 = x - BAR_W * (len(present) - 1) / 2
            for k, (m, v) in enumerate(present):
                ax.bar(x0 + k * BAR_W, v, width=BAR_W, align='center', zorder=3, **style(m))
                if v < 0.005:
                    ax.text(x0 + k * BAR_W, 0.012, '0', ha='center', va='bottom', fontsize=4.5, color=S.INK2, zorder=5)
        ax.axhline(0, color=S.AXIS, lw=0.6, zorder=4)
        ax.axvline((x_of[S.ORDER[n_dev - 1]] + x_of[S.ORDER[n_dev]]) / 2, color=S.AXIS, lw=0.6)
    top = axes[0]
    for xm, txt in (((x_of[S.ORDER[0]] + x_of[S.ORDER[n_dev - 1]]) / 2, 'Development genomes'),
                    ((x_of[S.ORDER[n_dev]] + x_of[S.ORDER[-1]]) / 2, 'Held-out genomes')):
        top.text(xm, 1.035, txt, transform=top.get_xaxis_transform(), ha='center', va='bottom', fontsize=7,
                 color=S.INK, fontweight='bold')
    for ax, letter in zip(axes, 'AB'):
        ax.text(-0.085, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
    axes[1].set_xticks([x_of[sp] for sp in S.ORDER])
    axes[1].set_xticklabels([S.italic(disp[sp]) for sp in S.ORDER], rotation=40, ha='right')
    axes[1].tick_params(axis='x', length=0)
    axes[1].set_xlim(-0.6, x_of[S.ORDER[-1]] + 0.6)
    shown = [m for m in S.METHODS if m != 'heliano' or any(c == 'f1_helitron' for c, _ in panels)]
    handles = [Patch(label=S.METHODS[m][0], **style(m)) for m in shown]
    top.legend(handles=handles, loc='lower center', ncol=len(shown), bbox_to_anchor=(0.5, 1.11), handlelength=1.4,
               handleheight=0.9, handletextpad=0.4, columnspacing=1.2, borderaxespad=0)
    S.save(fig, stem)


order_figure('figS1_orders_retro', [('f1_ltr', 'LTR retrotransposon F1'), ('f1_nonltr', 'Non-LTR retrotransposon F1')])
order_figure('figS2_orders_dna', [('f1_tir', 'TIR transposon F1'), ('f1_helitron', 'Helitron F1')])
