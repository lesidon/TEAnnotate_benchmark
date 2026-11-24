#!/usr/bin/env python3
"""Nature-styled figures for the whole-pipeline TE benchmark.

Reads results/benchmark_metrics.tsv (from aggregate.py) and results/concordance_all.tsv,
and renders — GOLD and PROXY strictly on separate axes (their metric types are not
comparable, per the soft-vs-hard scoring rule):

  gold_overall.<fmt>            per gold genome: F1 / sensitivity / precision by method
  gold_F1_by_category.<fmt>     per gold genome: F1 heatmap (TE category x method)
  proxy_coverage.<fmt>          per kingdom: genome coverage % by method (grouped bars)
  proxy_concordance.<fmt>       per proxy genome: pairwise Jaccard heatmap (method x method)
  proxy_lai.<fmt>               whole-genome LAI by genome (context)

Headless (Agg). Missing inputs are skipped with a warning.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 400,
    "font.size": 7, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 8, "axes.titleweight": "regular", "axes.labelsize": 7.5,
    "legend.fontsize": 6, "legend.frameon": False,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

# Okabe-Ito colourblind-safe palette; pipeline is emphasised in black.
_OKABE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
GOLD_CATEGORIES = ["ltr", "line", "sine", "nonltr", "tir", "mite", "helitron", "others", "any"]


def method_colors(methods):
    colors, k = {}, 0
    for m in methods:
        if m == "pipeline":
            colors[m] = "#000000"
        else:
            colors[m] = _OKABE[k % len(_OKABE)]; k += 1
    return colors


def _save(fig, outdir, stem, formats):
    for fmt in formats:
        fig.savefig(os.path.join(outdir, f"{stem}.{fmt}"), bbox_inches="tight")
    plt.close(fig)
    sys.stderr.write(f"[plot] wrote {stem}.{{{','.join(formats)}}}\n")


def _to_num(s):
    return pd.to_numeric(s, errors="coerce")


# ── GOLD ────────────────────────────────────────────────────────────────────
def gold_overall(master, outdir, formats):
    g = master[(master.scoring_type == "gold") & (master.category == "any")]
    if g.empty:
        return
    species = sorted(g.species.unique())
    shown = ["sensitivity", "precision", "F1"]
    fig, axes = plt.subplots(1, len(species), figsize=(2.6 * len(species) + 0.4, 2.6),
                             squeeze=False, constrained_layout=True)
    for ax, sp in zip(axes[0], species):
        sub = g[g.species == sp]
        methods = sorted(sub.method.unique())
        cmap = method_colors(methods)
        x = np.arange(len(shown)); w = 0.8 / max(len(methods), 1)
        for j, m in enumerate(methods):
            vals = [_to_num(sub[(sub.method == m) & (sub.metric_name == mt)].value).mean()
                    for mt in shown]
            bars = ax.bar(x + j * w, vals, w, label=m, color=cmap[m],
                          edgecolor="white", linewidth=0.3)
            ax.bar_label(bars, fmt="%.2f", padding=1, fontsize=4.5, rotation=90)
        ax.set_xticks(x + 0.4 - w / 2); ax.set_xticklabels(shown)
        ax.set_ylim(0, 1.05); ax.set_title(sp.replace("_", " "))
        ax.set_ylabel("score")
    axes[0][-1].legend(title="method", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle("Gold-standard scoring (EDTA lib-test.py)", fontsize=8)
    _save(fig, outdir, "gold_overall", formats)


def gold_F1_by_category(master, outdir, formats):
    g = master[(master.scoring_type == "gold") & (master.metric_name == "F1")]
    if g.empty:
        return
    species = sorted(g.species.unique())
    fig, axes = plt.subplots(1, len(species), figsize=(2.4 * len(species) + 0.6, 2.8),
                             squeeze=False, constrained_layout=True)
    for ax, sp in zip(axes[0], species):
        sub = g[g.species == sp]
        methods = sorted(sub.method.unique())
        cats = [c for c in GOLD_CATEGORIES if c in set(sub.category)]
        mat = np.full((len(cats), len(methods)), np.nan)
        for r, c in enumerate(cats):
            for cc, m in enumerate(methods):
                v = _to_num(sub[(sub.category == c) & (sub.method == m)].value)
                if len(v):
                    mat[r, cc] = v.mean()
        im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats)
        ax.set_title(sp.replace("_", " "))
        for r in range(len(cats)):
            for cc in range(len(methods)):
                if not np.isnan(mat[r, cc]):
                    ax.text(cc, r, f"{mat[r, cc]:.2f}", ha="center", va="center",
                            fontsize=4.5, color="white" if mat[r, cc] < 0.6 else "black")
    fig.colorbar(im, ax=axes[0][-1], fraction=0.046, pad=0.04, label="F1")
    fig.suptitle("Gold-standard F1 by TE category", fontsize=8)
    _save(fig, outdir, "gold_F1_by_category", formats)


# ── PROXY ─────────────────────────────────────────────────────────────────────
def proxy_coverage(master, outdir, formats):
    p = master[(master.scoring_type == "proxy") & (master.metric_name == "coverage_pct")]
    if p.empty:
        return
    p = p.assign(value=_to_num(p.value))
    kingdoms = sorted(p.kingdom.unique())
    fig, axes = plt.subplots(1, len(kingdoms), figsize=(2.8 * len(kingdoms) + 0.4, 2.7),
                             squeeze=False, constrained_layout=True)
    all_methods = sorted(p.method.unique()); cmap = method_colors(all_methods)
    for ax, kg in zip(axes[0], kingdoms):
        sub = p[p.kingdom == kg]
        species = sorted(sub.species.unique())
        x = np.arange(len(species)); w = 0.8 / max(len(all_methods), 1)
        for j, m in enumerate(all_methods):
            vals = [sub[(sub.species == sp) & (sub.method == m)].value.mean()
                    if not sub[(sub.species == sp) & (sub.method == m)].empty else np.nan
                    for sp in species]
            ax.bar(x + j * w, vals, w, label=m, color=cmap[m], edgecolor="white", linewidth=0.3)
        ax.set_xticks(x + 0.4 - w / 2)
        ax.set_xticklabels([s.replace("_", " ") for s in species], rotation=30, ha="right")
        ax.set_title(kg); ax.set_ylabel("genome masked (%)")
    axes[0][-1].legend(title="method", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle("Proxy metric: genome coverage by method", fontsize=8)
    _save(fig, outdir, "proxy_coverage", formats)


def proxy_concordance(concord, outdir, formats):
    if concord is None or concord.empty:
        return
    species = sorted(concord.species.unique())
    n = len(species); cols = min(3, n); rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(2.5 * cols, 2.5 * rows),
                             squeeze=False, constrained_layout=True)
    im = None
    for idx, sp in enumerate(species):
        ax = axes[idx // cols][idx % cols]
        sub = concord[concord.species == sp]
        methods = sorted(set(sub.method_a) | set(sub.method_b))
        mat = np.full((len(methods), len(methods)), np.nan)
        pos = {m: i for i, m in enumerate(methods)}
        for _, r in sub.iterrows():
            i, j = pos[r.method_a], pos[r.method_b]
            mat[i, j] = mat[j, i] = r.jaccard
        im = ax.imshow(mat, cmap="magma", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.set_yticks(range(len(methods))); ax.set_yticklabels(methods)
        ax.set_title(sp.replace("_", " "))
        for i in range(len(methods)):
            for j in range(len(methods)):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                            fontsize=4.5, color="white" if mat[i, j] < 0.5 else "black")
    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].axis("off")
    if im is not None:
        fig.colorbar(im, ax=axes[-1][-1], fraction=0.046, pad=0.04, label="Jaccard")
    fig.suptitle("Proxy metric: pairwise masked-bp concordance", fontsize=8)
    _save(fig, outdir, "proxy_concordance", formats)


def proxy_lai(master, outdir, formats):
    p = master[(master.scoring_type == "proxy") & (master.metric_name == "LAI")].copy()
    p["value"] = _to_num(p.value)
    p = p.dropna(subset=["value"])
    if p.empty:
        return
    p = p.sort_values(["kingdom", "species"])
    fig, ax = plt.subplots(figsize=(0.5 * len(p) + 1.5, 2.4), constrained_layout=True)
    bars = ax.bar(range(len(p)), p.value, color="#0072B2", edgecolor="white", linewidth=0.3)
    ax.bar_label(bars, fmt="%.1f", padding=1, fontsize=5)
    ax.set_xticks(range(len(p)))
    ax.set_xticklabels([s.replace("_", " ") for s in p.species], rotation=30, ha="right")
    ax.set_ylabel("LAI"); ax.set_title("Whole-genome LTR Assembly Index (context)")
    _save(fig, outdir, "proxy_lai", formats)


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=os.path.join(here, "results"))
    ap.add_argument("--outdir", default=None, help="default: <results>/figures")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args(argv)

    master_path = os.path.join(args.results, "benchmark_metrics.tsv")
    if not os.path.isfile(master_path):
        raise SystemExit(f"[plot] {master_path} not found — run aggregate.py first")
    master = pd.read_csv(master_path, sep="\t", dtype=str)
    concord_path = os.path.join(args.results, "concordance_all.tsv")
    concord = pd.read_csv(concord_path, sep="\t") if os.path.isfile(concord_path) else None

    outdir = args.outdir or os.path.join(args.results, "figures")
    os.makedirs(outdir, exist_ok=True)

    gold_overall(master, outdir, args.formats)
    gold_F1_by_category(master, outdir, args.formats)
    proxy_coverage(master, outdir, args.formats)
    proxy_concordance(concord, outdir, args.formats)
    proxy_lai(master, outdir, args.formats)
    print(f"[plot] figures -> {outdir}")


if __name__ == "__main__":
    main()
