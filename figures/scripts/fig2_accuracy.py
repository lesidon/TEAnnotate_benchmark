"""Fig. 2 (accuracy) - library accuracy on the 6 development and 9 held-out genomes,
as grouped bars from a zero baseline.

(A) harmonised all-repeat F1 (lib-test any_dontcare); (B) perfectly recovered reference families as a
fraction of the reference families (Flynn et al. 2020). Data: data/accuracy.csv (extract_data.py).
Within each genome the bars keep a fixed method order; methods not run on a genome are left out and the
remaining bars are centred. TEAnnotate without TEtrimmer is hatched (texture, not a second hue).
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import figstyle as S

S.setup()
mpl.rcParams['hatch.linewidth'] = 0.55
acc = S.read('accuracy.csv')
by = {(r['species'], r['method']): r for r in acc}
disp = {r['species']: r['display'] for r in S.read('genomes.csv')}

BAR_W = 0.135
fig, axes = plt.subplots(2, 1, figsize=(6.4, 4.4), sharex=True, gridspec_kw={'hspace': 0.14})
panels = [('f1', 'Harmonised all-repeat F1', (0, 1.0), 0.2),
          ('perfect_frac', 'Perfect families\n(fraction of reference)', (0, 0.6), 0.1)]
n_dev = 6
x_of = {sp: i + (0.6 if i >= n_dev else 0) for i, sp in enumerate(S.ORDER)}   # gap between the two groups


def style(m):
    label, col, _mk, _filled, _size = S.METHODS[m]
    if m == 'teannotate_wo':
        return dict(facecolor='white', edgecolor=col, hatch='//////', linewidth=0.6)
    return dict(facecolor=col, edgecolor='white', linewidth=0.4)


def value(r, metric):
    if metric == 'f1':
        return S.num(r['f1'])
    if r['perfect'] in ('', None):
        return None
    return int(r['perfect']) / int(r['n_reference'])


for ax, (metric, ylabel, ylim, step) in zip(axes, panels):
    ax.set_ylim(*ylim)
    ax.set_yticks([round(i * step, 2) for i in range(int(round(ylim[1] / step)) + 1)])
    ax.set_ylabel(ylabel)
    ax.grid(axis='x', visible=False)
    for sp in S.ORDER:
        x = x_of[sp]
        if sp in S.PENDING:
            ax.axvspan(x - 0.44, x + 0.44, color=S.GRID, alpha=0.6, lw=0, zorder=0)
            ax.text(x, ylim[0] + 0.5 * (ylim[1] - ylim[0]), 'not yet scored', rotation=90, ha='center',
                    va='center', fontsize=6, color=S.INK2)
            continue
        present = [(m, value(by[(sp, m)], metric)) for m in S.METHODS if (sp, m) in by]
        present = [(m, v) for m, v in present if v is not None]
        x0 = x - BAR_W * (len(present) - 1) / 2
        for k, (m, v) in enumerate(present):
            ax.bar(x0 + k * BAR_W, v, width=BAR_W, align='center', zorder=3, **style(m))
            if v == 0 or (metric == 'f1' and v < 0.005):   # label zero or rounded-zero F1
                ax.text(x0 + k * BAR_W, ylim[1] * 0.012, '0', ha='center', va='bottom', fontsize=4.5,
                        color=S.INK2, zorder=5)
    ax.axhline(0, color=S.AXIS, lw=0.6, zorder=4)
    ax.axvline((x_of[S.ORDER[n_dev - 1]] + x_of[S.ORDER[n_dev]]) / 2, color=S.AXIS, lw=0.6)

top = axes[0]
mid_dev = (x_of[S.ORDER[0]] + x_of[S.ORDER[n_dev - 1]]) / 2
mid_ho = (x_of[S.ORDER[n_dev]] + x_of[S.ORDER[-1]]) / 2
for xm, txt in ((mid_dev, 'Development genomes'), (mid_ho, 'Held-out genomes')):
    top.text(xm, 1.035, txt, transform=top.get_xaxis_transform(), ha='center', va='bottom', fontsize=7,
             color=S.INK, fontweight='bold')
for ax, letter in zip(axes, 'AB'):
    ax.text(-0.085, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

axes[1].set_xticks([x_of[sp] for sp in S.ORDER])
axes[1].set_xticklabels([S.italic(disp[sp]) for sp in S.ORDER], rotation=40, ha='right')
axes[1].tick_params(axis='x', length=0)
axes[1].set_xlim(-0.6, x_of[S.ORDER[-1]] + 0.6)

handles = [Patch(label=S.METHODS[m][0], **style(m)) for m in S.METHODS]
top.legend(handles=handles, loc='lower center', ncol=6, bbox_to_anchor=(0.5, 1.11), handlelength=1.4,
           handleheight=0.9, handletextpad=0.4, columnspacing=1.2, borderaxespad=0)
S.save(fig, 'fig2_accuracy')
