"""Fig. 3 - classifier-only benchmark of the deep-learning TE classifiers and the two ensemble votes.

Leakage-filtered test halves (PanTEon v1.5.1; fungi also MycoMobilome v1.1), per kingdom, at order
level (top) and superfamily level (bottom; only test sequences whose truth names a superfamily).
Filled marker: support-weighted F1 with its 95% cluster-bootstrap CI (the tuning objective);
small open marker: macro-F1 of the same configuration. The two-step vote is the pipeline default;
its first step is the weighted vote, so both votes coincide at order level.
Data: data/classifiers.csv (extract_data.py, from ensemble_v3_dfam37/test/*.wf1.two_step.tsv).
"""
import csv
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as S

S.setup()
rows = S.read('classifiers.csv')
METHODS = ['TEclass2', 'DeepTE', 'CREATE', 'NeuralTE', 'Terrier', 'Weighted vote', 'Two-step vote']
STYLE = {m: (S.SLATE, 'o', 22) for m in METHODS[:5]}
STYLE['Weighted vote'] = (S.GLACIER, 's', 26)
STYLE['Two-step vote'] = (S.MIDNIGHT, 'o', 36)
KING = [('plant', 'Plants'), ('animal', 'Animals'), ('fungi', 'Fungi')]
LEVEL = [('order', 'Order'), ('sf', 'Superfamily')]
by = {(r['kingdom'], r['level'], r['method']): r for r in rows}
# value labels of the two-step vote and the best single classifier: 6-decimal point estimates of
# paired_boot.py (the test files carry 4 decimals, which can round differently from Table S6)
PRECISE = {}
for r in csv.DictReader(open(os.path.join(S.DATA, 'classifier_paired_bootstrap.tsv'), newline=''), delimiter='\t'):
    lv = 'order' if r['level'] == 'order' else 'sf'
    PRECISE[(r['kingdom'], lv, 'Two-step vote')] = float(r['vote_wf1'])
    PRECISE[(r['kingdom'], lv, r['best_single'])] = float(r['single_wf1'])

fig, axes = plt.subplots(2, 3, figsize=(6.4, 3.7), sharex=True, sharey=True,
                         gridspec_kw={'hspace': 0.18, 'wspace': 0.06})
ys = {m: i for i, m in enumerate(METHODS)}
for i, (lv, lvname) in enumerate(LEVEL):
    for j, (k, kname) in enumerate(KING):
        ax = axes[i, j]
        ax.axhspan(4.5, 6.5, color=S.GRID, alpha=0.45, lw=0, zorder=0)
        singles = [by[(k, lv, m)] for m in METHODS[:5]]
        best_single = max(singles, key=lambda r: S.num(r['wf1']))['method']
        for m in METHODS:
            r = by[(k, lv, m)]
            col, mk, size = STYLE[m]
            w, lo, hi, mac = (S.num(r[c]) for c in ('wf1', 'wf1_lo', 'wf1_hi', 'macro_f1'))
            y = ys[m]
            ax.hlines(y, lo, hi, color=col, lw=1.4, zorder=3)
            ax.scatter(w, y, s=size, marker=mk, color=col, edgecolor='white', lw=0.5, zorder=4)
            ax.scatter(mac, y, s=12, marker='D', facecolor='white', edgecolor=col, lw=0.7, zorder=3)
            if m in ('Two-step vote', best_single) or (m == 'Weighted vote' and lv == 'sf'):
                ax.text(1.015, y, f'{PRECISE.get((k, lv, m), w):.3f}', va='center', ha='left',
                        fontsize=5.6, color=S.INK if m == 'Two-step vote' else S.INK2, clip_on=False)
        ax.set_xlim(0.0, 1.16)
        ax.set_ylim(len(METHODS) - 0.4, -0.6)
        ax.grid(axis='y', visible=False)
        ax.tick_params(axis='y', length=0)
        if i == 0:
            ax.set_title(kname, fontsize=7.5, fontweight='bold', pad=4)
        if j == 0:
            ax.set_ylabel(lvname, fontsize=7.5, fontweight='bold', labelpad=6)
axes[0, 0].set_yticks(range(len(METHODS)))
axes[0, 0].set_yticklabels(METHODS)
for ax in axes[1]:
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
fig.supxlabel('F1 on the held-out test half', fontsize=7, y=0.005)

handles = [Line2D([], [], color=S.INK2, marker='o', ms=4, lw=1.4, markeredgecolor='white',
                  label='Support-weighted F1, 95% CI'),
           Line2D([], [], ls='', marker='D', ms=3.4, markerfacecolor='white', markeredgecolor=S.INK2,
                  label='Macro-F1'),
           Line2D([], [], ls='', marker='s', ms=4.2, color=S.GLACIER, markeredgecolor='white', label='Weighted vote'),
           Line2D([], [], ls='', marker='o', ms=5, color=S.MIDNIGHT, markeredgecolor='white',
                  label='Two-step vote (pipeline default)')]
fig.legend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.04), handletextpad=0.3,
           columnspacing=1.2, fontsize=6.2)
S.save(fig, 'fig3_classifiers')
