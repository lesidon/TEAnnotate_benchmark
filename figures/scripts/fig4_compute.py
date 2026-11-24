"""Fig. 4 - computational cost per stage and TEtrimmer throughput, genomes ordered by size.

(A) CPU-hours measured by Nextflow (realtime x %cpu) for each discovery tool of the route, TEtrimmer and
all other library stages (deduplication, evidence gate, post-processing, classification); the vertical
tick marks their sum. Masking is not included. (B) TEtrimmer throughput: gated consensi curated per
wall-clock hour of the TETRIMMER task.
Development genomes: standalone timed runs of each discovery tool and the route's library stages with TEtrimmer;
held-out genomes: the frozen discovery runs and the final library-stage runs.
Z. mays shows measured stages only. Its chunked TEtrimmer run is omitted from both panels.
Data: data/compute.csv, accuracy.csv.
"""
import collections
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator
import figstyle as S

S.setup()
comp = S.read('compute.csv')
acc = {(r['species'], r['method']): r for r in S.read('accuracy.csv')}
gen = {r['species']: r for r in S.read('genomes.csv')}

stages = collections.defaultdict(dict)
for r in comp:
    if r['species'] in S.PENDING:
        continue
    stages[r['species']][r['stage']] = (S.num(r['cpu_h']), S.num(r['wall_h']))

order = sorted(S.ORDER, key=lambda sp: int(gen[sp]['genome_bp']))
STYLE = {
    'earlgrey': S.METHODS['earlgrey'][1:4] + (22,), 'hite': S.METHODS['hite'][1:4] + (26,),
    'edta': S.METHODS['edta'][1:4] + (18,), 'heliano': S.METHODS['heliano'][1:4] + (22,),
    'tetrimmer': (S.MIDNIGHT, 'o', True, 30), 'library_other': (S.MIDNIGHT, 'o', False, 26),
}
LABEL = {'earlgrey': 'Earl Grey', 'hite': 'HiTE', 'edta': 'EDTA', 'heliano': 'HELIANO', 'tetrimmer': 'TEtrimmer',
         'library_other': 'Other library stages'}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.4, 3.6), sharey=True,
                             gridspec_kw={'width_ratios': [3.1, 1.0], 'wspace': 0.06})
ylab = []
for y, sp in enumerate(order):
    g = gen[sp]
    size = int(g['genome_bp']) / 1e9
    role = '' if g['role'] == 'heldout' else ' *'
    ylab.append(S.italic(g['display']) + f'  {size:.2f} Gb' + role)
    if sp in S.PENDING:
        for a in (ax, bx):
            a.axhspan(y - 0.42, y + 0.42, color=S.GRID, alpha=0.6, lw=0, zorder=0)
        ax.text(2.0, y, 'pending', va='center', fontsize=6, color=S.INK2)
        continue
    st = stages[sp]
    total = sum(v[0] for v in st.values())
    xs = [v[0] for v in st.values()]
    ax.hlines(y, min(xs), max(xs), color=S.GRID, lw=0.8, zorder=1)
    if 'tetrimmer' in st:
        ax.scatter(total, y, marker='|', s=60, color=S.INK, lw=1.0, zorder=5)
    for k, (cpu, wall) in st.items():
        col, mk, filled, size_ = STYLE[k]
        dy = 0.14 if k in ('tetrimmer', 'library_other') else -0.14   # tools above, pipeline stages below
        ax.scatter(cpu, y + dy, s=size_, marker=mk, facecolor=col if filled else 'white',
                   edgecolor=col if not filled else 'white', lw=0.9 if not filled else 0.4, zorder=4)
    if 'tetrimmer' in st:
        n_in = int(acc[(sp, 'teannotate_wo')]['n_library'])
        rate = n_in / st['tetrimmer'][1]
        bx.hlines(y, 1, rate, color=S.MIDNIGHT, lw=1.2, alpha=0.35, zorder=2)
        bx.scatter(rate, y, s=24, marker='o', color=S.MIDNIGHT, edgecolor='white', lw=0.4, zorder=3)
        bx.text(rate * 1.25, y, f'{rate:.0f}', va='center', fontsize=5.5, color=S.INK2)
    elif sp == 'zea_mays':
        bx.text(2, y, 'chunked run', va='center', fontsize=5.5, color=S.INK2)

for a in (ax, bx):
    a.set_xscale('log')
    a.grid(axis='y', visible=False)
    a.tick_params(axis='y', length=0)
ax.set_xlim(0.01, 5000)
ax.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10, 100, 1000]))
ax.set_xticklabels(['0.01', '0.1', '1', '10', '100', '1000'])
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlabel('CPU-hours (log scale)')
bx.set_xlim(1, 1000)
bx.xaxis.set_major_locator(FixedLocator([1, 10, 100, 1000]))
bx.set_xticklabels(['1', '10', '100', '1000'])
bx.xaxis.set_minor_locator(NullLocator())
bx.set_xlabel('TEtrimmer consensi per hour')
ax.set_yticks(range(len(order)))
ax.set_yticklabels(ylab, fontsize=6.3)
ax.set_ylim(len(order) - 0.5, -0.5)
for a, letter, dx in ((ax, 'A', -0.36), (bx, 'B', -0.08)):
    a.text(dx, 1.06, letter, transform=a.transAxes, fontsize=9, fontweight='bold', va='bottom')

handles = [Line2D([], [], ls='', marker=STYLE[k][1], ms=4.4, markerfacecolor=STYLE[k][0] if STYLE[k][2] else 'white',
                  markeredgecolor=STYLE[k][0] if not STYLE[k][2] else 'white',
                  markeredgewidth=0.9 if not STYLE[k][2] else 0.4, label=LABEL[k]) for k in STYLE]
handles.append(Line2D([], [], ls='', marker='|', ms=6, color=S.INK, markeredgewidth=1.0, label='Sum'))
fig.legend(handles=handles, loc='upper center', ncol=7, bbox_to_anchor=(0.5, 1.04), handletextpad=0.2,
           columnspacing=1.0, fontsize=6)
fig.text(0.0, -0.02, '* development genome', fontsize=5.5, color=S.MUTED)
S.save(fig, 'fig4_compute')
