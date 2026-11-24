#!/usr/bin/env python3
"""
plot_benchmarks.py — publication-ready figures for the TE-classification
benchmark. Ingests the TSVs written by evaluate_classification.py and
calculate_curves.py and renders, per lineage tag:

  <tag>.roc_pr.<fmt>          ROC + PR overlays (all models, macro or micro avg)
  <tag>.confusion.<fmt>       row-normalised confusion-matrix heatmaps (grid)
  <tag>.metrics_bar.<fmt>     MCC / macro-F1 bar chart across models
  <tag>.benchmark_panel.<fmt> combined side-by-side panel (ROC | PR | confusions)

All confusion matrices for a tag share ONE common, canonically-ordered class
axis (union across models), so the grid is directly comparable. Rows with no
ground-truth support are shown grey (not annotated) instead of producing
division-by-zero artefacts.

Runs headless (Agg backend) for cluster nodes. Missing inputs are skipped with a
warning rather than aborting the whole figure set.

Example
-------
  python plot_benchmarks.py --indir results_eval results_curves \
      --tag dfam_fungi --curve-average macro --formats png pdf
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import te_label_utils as tu               # noqa: E402

LOG = tu.setup_logger("plot_benchmarks")

# ---- Nature-like style ---------------------------------------------------- #
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "font.size": 7,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 8,
    "axes.titleweight": "regular",
    "axes.labelsize": 7.5,
    "legend.fontsize": 6,
    "legend.frameon": False,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 1.1,
    # keep text editable in vector output (Nature submission requirement)
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})
_PALETTE = plt.get_cmap("tab10").colors
_CONF_CMAP = "Blues"        # sequential, zero -> white, prints cleanly


# --------------------------------------------------------------------------- #
# Discovery
# --------------------------------------------------------------------------- #
KINDS = {
    "roc_coords": "roc_coords.tsv",
    "pr_coords": "pr_coords.tsv",
    "curve_auc": "curve_auc.tsv",
    "confusion_matrix": "confusion_matrix.tsv",
    "summary": "all_models.summary_metrics.tsv",
}


def discover(indirs: List[str]) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Return {tag: {kind: {model: path}}}  (+ tag-level 'summary')."""
    found: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(
        lambda: defaultdict(dict))
    for d in indirs:
        for kind, suffix in KINDS.items():
            for path in glob.glob(os.path.join(d, f"*.{suffix}")):
                base = os.path.basename(path)[:-(len(suffix) + 1)]  # strip '.<suffix>'
                if kind == "summary":
                    tag = base[:-len(".all_models")] if base.endswith(
                        ".all_models") else base
                    found[tag]["summary"]["_"] = path
                    continue
                tag, _, model = base.rpartition(".")
                if not tag:
                    tag, model = model, "model"
                found[tag][kind][model] = path
    return found


def sort_models(models: List[str]) -> List[str]:
    """Ensemble last, otherwise alphabetical (stable legend/colour order)."""
    ens = [m for m in models if m.lower() in ("ensemble", "vote", "consensus")]
    rest = sorted(m for m in models if m not in ens)
    return rest + ens


def model_colors(models: List[str]) -> Dict[str, tuple]:
    return {m: _PALETTE[i % len(_PALETTE)] for i, m in enumerate(models)}


# --------------------------------------------------------------------------- #
# ROC / PR curves
# --------------------------------------------------------------------------- #
def _avg_rows(df: pd.DataFrame, average: str) -> pd.DataFrame:
    return df[df["class"] == average]


# How each curve was scored (from calculate_curves' 'score_type' column) -> the
# legend tag. A hard single-point curve must never sit unlabelled next to a soft
# vote-fraction curve (their AUC/AP are not comparable quantities).
_SCORE_LABEL = {"hard": "hard", "vote-fraction": "vote-frac", "proba": "prob"}


def _score_suffix(adf: pd.DataFrame) -> str:
    if "score_type" in adf.columns and len(adf):
        st = str(adf["score_type"].iloc[0])
        return " [%s]" % _SCORE_LABEL.get(st, st)
    return ""


def _warn_mixed_scores(tag: str, auc_data: Dict[str, str]) -> None:
    """Warn loudly if one figure would overlay curves scored different ways.

    A 'vote-fraction' curve integrates over thresholds; a 'hard' curve is a single
    operating point. Their AUC/AP are not comparable, so overlaying them (even with
    legend tags) is a red flag worth surfacing.
    """
    by_type: Dict[str, List[str]] = defaultdict(list)
    for model, path in auc_data.items():
        try:
            adf = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if "score_type" in adf.columns and len(adf):
            by_type[str(adf["score_type"].iloc[0])].append(model)
    if len(by_type) > 1:
        desc = "; ".join("%s=%s" % (k, ",".join(v)) for k, v in by_type.items())
        LOG.warning(
            "tag '%s' MIXES curve score types (%s). AUC/AP are NOT comparable "
            "across types — legend tags [hard]/[vote-frac]/[prob] mark this; do "
            "not read them as a head-to-head.", tag, desc)


def plot_roc(ax, tag_data: Dict[str, str], auc_data: Dict[str, str],
             average: str, colors: Dict[str, tuple]) -> None:
    ax.plot([0, 1], [0, 1], ls=":", c="0.6", lw=0.8, zorder=1)
    for model in sort_models(list(tag_data.keys())):
        df = pd.read_csv(tag_data[model], sep="\t")
        sub = _avg_rows(df, average).sort_values("fpr")
        if sub.empty:
            continue
        auc_val, suffix = np.nan, ""
        if model in auc_data:
            adf = pd.read_csv(auc_data[model], sep="\t")
            row = adf[adf["class"] == average]
            if not row.empty:
                auc_val = float(row["roc_auc"].iloc[0])
            suffix = _score_suffix(adf)
        ax.plot(sub["fpr"], sub["tpr"], color=colors[model],
                label=f"{model}{suffix} ({auc_val:.3f})")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(f"ROC ({average}-average)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", title="AUC", title_fontsize=6)


def plot_pr(ax, tag_data: Dict[str, str], auc_data: Dict[str, str],
            average: str, colors: Dict[str, tuple],
            baseline: Optional[float]) -> None:
    if baseline is not None:
        ax.axhline(baseline, ls=":", c="0.6", lw=0.8)
    for model in sort_models(list(tag_data.keys())):
        df = pd.read_csv(tag_data[model], sep="\t")
        sub = _avg_rows(df, average).sort_values("recall")
        if sub.empty:
            continue
        ap_val, suffix = np.nan, ""
        if model in auc_data:
            adf = pd.read_csv(auc_data[model], sep="\t")
            row = adf[adf["class"] == average]
            if not row.empty:
                ap_val = float(row["pr_auc"].iloc[0])
            suffix = _score_suffix(adf)
        ax.plot(sub["recall"], sub["precision"], color=colors[model],
                label=f"{model}{suffix} ({ap_val:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"Precision-Recall ({average}-average)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.legend(loc="lower left", title="AP", title_fontsize=6)


# --------------------------------------------------------------------------- #
# Confusion matrices (shared class axis)
# --------------------------------------------------------------------------- #
def _ordered_classes(all_labels) -> List[str]:
    """Canonical order (class/order first, Unknown/Other last), then extras."""
    preferred = list(tu.CANONICAL_CLASSES) + list(tu.CANONICAL_ORDERS)
    out: List[str] = []
    seen = set()
    for c in preferred:
        if c in all_labels and c not in seen:
            out.append(c)
            seen.add(c)
    for c in sorted(all_labels):
        if c not in seen:
            out.append(c)
            seen.add(c)
    return out


def _load_confusion(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = [str(x) for x in df.index]
    df.columns = [str(x) for x in df.columns]
    return df


def _row_normalise(counts: np.ndarray) -> np.ndarray:
    """Row-fraction matrix; zero-support rows -> NaN (rendered as 'no data')."""
    rs = counts.sum(axis=1, keepdims=True)
    frac = np.full(counts.shape, np.nan, dtype=float)
    np.divide(counts, rs, out=frac, where=rs > 0)   # out= avoids uninit garbage
    return frac


def _draw_confusion(ax, disp: np.ndarray, classes: List[str], title: str,
                    show_x: bool, show_y: bool, vmax: float, is_fraction: bool,
                    equal_aspect: bool = True):
    cmap = plt.get_cmap(_CONF_CMAP).copy()
    cmap.set_bad("0.9")                              # grey for no-support rows
    data = np.ma.masked_invalid(disp)
    im = ax.imshow(data, cmap=cmap, vmin=0, vmax=vmax,
                   aspect="equal" if equal_aspect else "auto")
    n = len(classes)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    if show_x:
        ax.set_xticklabels(classes, rotation=45, ha="right",
                           rotation_mode="anchor")
        ax.set_xlabel("Predicted")
    else:
        ax.set_xticklabels([])
    if show_y:
        ax.set_yticklabels(classes)
        ax.set_ylabel("True")
    else:
        ax.set_yticklabels([])
    ax.set_title(title, pad=3)

    # Crisp white cell borders, no axis spines/ticks (Nature heatmap look).
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    if n <= 16:
        thr = 0.55 * vmax
        for i in range(n):
            for j in range(n):
                v = disp[i, j]
                if np.isnan(v) or v <= (0.004 if is_fraction else 0.5):
                    continue
                txt = f"{v:.2f}" if is_fraction else f"{int(round(v))}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=5.5,
                        color="white" if v >= thr else "0.15")
    return im


def _confusion_grid(conf_paths: Dict[str, str], normalize: bool
                    ) -> Tuple[List[str], List[str], Dict[str, np.ndarray], bool, float]:
    """Load + align every model's matrix onto one common class axis."""
    models = sort_models(list(conf_paths))
    dfs = {m: _load_confusion(conf_paths[m]) for m in models}
    all_labels = set()
    for df in dfs.values():
        all_labels.update(df.index)
        all_labels.update(df.columns)
    classes = _ordered_classes(all_labels)

    disp: Dict[str, np.ndarray] = {}
    vmax = 1.0 if normalize else 0.0
    for m in models:
        counts = (dfs[m].reindex(index=classes, columns=classes)
                  .fillna(0).to_numpy(dtype=float))
        if normalize:
            disp[m] = _row_normalise(counts)
        else:
            disp[m] = counts
            vmax = max(vmax, float(np.nanmax(counts)) if counts.size else 0.0)
    return models, classes, disp, normalize, (vmax if vmax > 0 else 1.0)


def plot_confusion_grid(conf_paths: Dict[str, str], tag: str, normalize: bool):
    """Standalone confusion grid figure (own colorbar)."""
    models, classes, disp, is_frac, vmax = _confusion_grid(conf_paths, normalize)
    n = len(models)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    cell = max(1.7, 0.28 * len(classes) + 0.9)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(cell * ncol + 1.1, cell * nrow + 0.4),
                             squeeze=False, constrained_layout=True)
    im = None
    for k, m in enumerate(models):
        r, c = divmod(k, ncol)
        show_x = ((r + 1) * ncol + c) >= n         # bottom-most used axis in col
        show_y = (c == 0)
        im = _draw_confusion(axes[r][c], disp[m], classes, m, show_x, show_y,
                             vmax, is_frac)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    if im is not None:
        cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.65,
                          aspect=30, pad=0.015)
        cb.set_label("row fraction" if is_frac else "count")
        cb.outline.set_visible(False)
    # No suptitle: with the tight grid it overlaps the top-row subplot titles.
    # Identify the figure from its filename / the paper caption instead.
    return fig


# --------------------------------------------------------------------------- #
# Aggregate metric bars
# --------------------------------------------------------------------------- #
def plot_metrics_bar(ax, summary_path: str) -> None:
    df = pd.read_csv(summary_path, sep="\t").sort_values("mcc")
    x = np.arange(len(df))
    w = 0.4
    b1 = ax.barh(x - w / 2, df["mcc"], height=w, label="MCC", color=_PALETTE[0])
    b2 = ax.barh(x + w / 2, df["macro_f1"], height=w, label="Macro-F1",
                 color=_PALETTE[1])
    # Exact value at the end of every bar (matplotlib >= 3.4).
    for bars in (b1, b2):
        ax.bar_label(bars, fmt="%.3f", padding=2, fontsize=5.5, color="0.2")
    ax.set_yticks(x)
    ax.set_yticklabels(df["model"])
    ax.set_xlabel("Score")
    # Leave head-room on the right so value labels never clip the frame; extend
    # left only if any metric is negative (MCC can be).
    lo = min(0.0, float(df[["mcc", "macro_f1"]].to_numpy().min()))
    ax.set_xlim(lo - (0.08 if lo < 0 else 0.0), 1.18)
    ax.margins(y=0.02)
    ax.set_title("Aggregate performance")
    ax.legend(loc="lower right")


def _save(fig, outdir: str, stem: str, formats: List[str]) -> None:
    for fmt in formats:
        p = os.path.join(outdir, f"{stem}.{fmt}")
        fig.savefig(p, bbox_inches="tight")
        LOG.info("wrote %s", p)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot ROC/PR curves and confusion matrices for the TE "
                    "classification benchmark.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--indir", nargs="+", required=True,
                   help="Directory(ies) with evaluate_/calculate_ TSV outputs.")
    p.add_argument("--tag", default=None,
                   help="Only plot this lineage tag (default: all discovered).")
    p.add_argument("--curve-average", default="macro", choices=["macro", "micro"])
    p.add_argument("--formats", nargs="+", default=["png", "pdf"],
                   help="Output image formats.")
    p.add_argument("--outdir", default="benchmark_figures")
    p.add_argument("--no-normalize-confusion", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.ensure_dir(args.outdir)
    catalog = discover(args.indir)
    if not catalog:
        raise SystemExit(f"ERROR: no benchmark TSVs found under {args.indir}")

    normalize = not args.no_normalize_confusion
    tags = [args.tag] if args.tag else sorted(catalog.keys())
    for tag in tags:
        data = catalog.get(tag)
        if not data:
            LOG.warning("no data for tag '%s'", tag)
            continue
        roc = data.get("roc_coords", {})
        pr = data.get("pr_coords", {})
        auc = data.get("curve_auc", {})
        conf = data.get("confusion_matrix", {})
        summary = data.get("summary", {}).get("_")
        models = sort_models(list(set(roc) | set(pr) | set(conf)))
        colors = model_colors(models)
        avg = args.curve_average
        LOG.info("tag '%s': models=%s", tag, ", ".join(models))
        _warn_mixed_scores(tag, auc)

        # 1. ROC + PR side by side.
        if roc or pr:
            fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5),
                                     constrained_layout=True)
            if roc:
                plot_roc(axes[0], roc, auc, avg, colors)
            if pr:
                plot_pr(axes[1], pr, auc, avg, colors, baseline=None)
            fig.suptitle(f"{tag} — classification curves")
            _save(fig, args.outdir, f"{tag}.roc_pr", args.formats)

        # 2. Confusion heatmaps grid (shared class axis).
        if conf:
            fig = plot_confusion_grid(conf, tag, normalize)
            _save(fig, args.outdir, f"{tag}.confusion", args.formats)

        # 3. Aggregate metric bars.
        if summary and os.path.exists(summary):
            fig, ax = plt.subplots(figsize=(4.5, max(2.4, 0.5 * len(models))),
                                   constrained_layout=True)
            plot_metrics_bar(ax, summary)
            _save(fig, args.outdir, f"{tag}.metrics_bar", args.formats)

        # 4. Combined panel (ROC | PR on top; confusions below, shared axis).
        if (roc or pr) and conf:
            models_c, classes, disp, is_frac, vmax = _confusion_grid(
                conf, normalize)
            ncol = max(2, min(4, len(models_c)))
            fig = plt.figure(figsize=(2.9 * ncol, 6.6), constrained_layout=True)
            gs = fig.add_gridspec(2, ncol, height_ratios=[1.15, 1.0])
            half = max(1, ncol // 2)
            ax_roc = fig.add_subplot(gs[0, :half])
            ax_pr = fig.add_subplot(gs[0, half:])
            if roc:
                plot_roc(ax_roc, roc, auc, avg, colors)
            if pr:
                plot_pr(ax_pr, pr, auc, avg, colors, baseline=None)
            im = None
            for k, m in enumerate(models_c[:ncol]):
                ax = fig.add_subplot(gs[1, k])
                im = _draw_confusion(ax, disp[m], classes, m,
                                     show_x=True, show_y=(k == 0),
                                     vmax=vmax, is_fraction=is_frac,
                                     equal_aspect=False)
            if im is not None:
                cb = fig.colorbar(im, ax=fig.axes[-len(models_c[:ncol]):],
                                  shrink=0.6, aspect=25, pad=0.01)
                cb.set_label("row fraction" if is_frac else "count")
                cb.outline.set_visible(False)
            fig.suptitle(f"{tag} — TE classification benchmark")
            _save(fig, args.outdir, f"{tag}.benchmark_panel", args.formats)

    LOG.info("All figures written to %s/", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
