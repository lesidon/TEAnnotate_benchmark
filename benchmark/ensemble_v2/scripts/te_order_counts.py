#!/usr/bin/env python3
"""
te_order_counts.py — how many TE families by ORDER, per kingdom.

Reads the per-kingdom ground-truth tables written by panteon_labels.py
(panteon_<fungi|plant|animal>.heldout.truth.tsv), counts families by TE order
(LTR / LINE / SINE / DNA / RC / DIRS / PLE / Maverick …), and renders a
Nature-styled grouped horizontal bar chart plus a tidy counts table.

Form  : magnitude comparison over a categorical axis (orders) grouped by kingdom
        -> grouped horizontal bars.
Colour: identity = kingdom, assigned in fixed order from the colourblind-safe
        Okabe-Ito palette (never cycled); text stays in ink tokens, not the bar hue.

Example
-------
  python te_order_counts.py --indir truth --outdir figs --formats png pdf
  # or explicit files:
  python te_order_counts.py \
      --truth animal=truth/panteon_animal.heldout.truth.tsv \
              plant=truth/panteon_plant.heldout.truth.tsv \
              fungi=truth/panteon_fungi.heldout.truth.tsv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

import te_label_utils as tu               # noqa: E402

LOG = tu.setup_logger("te_order_counts")

# Okabe-Ito, colourblind-safe. Fixed hue order; colour follows the KINGDOM entity.
KINGDOM_COLORS = OrderedDict([
    ("animal", "#0072B2"),   # blue
    ("plant",  "#009E73"),   # green
    ("fungi",  "#E69F00"),   # orange
])
_EXTRA_HUES = ["#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999"]
_KINGDOM_ALIAS = {"fungus": "fungi", "fungal": "fungi", "animalia": "animal",
                  "metazoa": "animal", "plantae": "plant", "viridiplantae": "plant"}

# ---- Nature-like style (matches plot_benchmarks.py) ----------------------- #
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 400,
    "font.size": 7, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 8, "axes.labelsize": 8,
    "legend.fontsize": 7, "legend.frameon": False,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def canon_kingdom(name: str) -> str:
    k = name.strip().lower()
    return _KINGDOM_ALIAS.get(k, k)


def discover(indir: str) -> "OrderedDict[str, str]":
    """Find panteon_<kingdom>.heldout.truth.tsv (or *.truth.tsv) under indir."""
    out: "OrderedDict[str, str]" = OrderedDict()
    pats = ["panteon_*.heldout.truth.tsv", "panteon_*.truth.tsv", "*.truth.tsv"]
    seen = set()
    for pat in pats:
        for path in sorted(glob.glob(os.path.join(indir, pat))):
            if path in seen:
                continue
            seen.add(path)
            m = re.search(r"panteon_([A-Za-z]+)", os.path.basename(path)) or \
                re.match(r"([A-Za-z]+)", os.path.basename(path))
            king = canon_kingdom(m.group(1)) if m else os.path.basename(path)
            out.setdefault(king, path)
    return out


def parse_sources(items: List[str]) -> "OrderedDict[str, str]":
    out: "OrderedDict[str, str]" = OrderedDict()
    for it in items:
        if "=" not in it:
            raise SystemExit(f"ERROR: --truth entry '{it}' is not KINGDOM=FILE")
        name, path = it.split("=", 1)
        out[canon_kingdom(name)] = path
    return out


def count_orders(path: str, order_col: str, level: str, keep_nonte: bool
                 ) -> "Counter[str]":
    """Count families by normalised order in one truth TSV."""
    counts: "Counter[str]" = Counter()
    with tu.open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        # Prefer the RAW label and re-normalise live, so te_label_utils fixes
        # (e.g. Helitron -> RC) apply even to TSVs whose precomputed label_order
        # column is stale. Only fall back to label_order if no raw label exists.
        col = tu.detect_column(cols, [order_col, "label", "classification",
                                      "rm_type", "label_order", "class"])
        if col is None:
            raise SystemExit(f"ERROR: no order/label column in {path}: {cols}")
        for row in reader:
            order = tu.normalize_label(row.get(col, ""), level)
            if not keep_nonte and order in ("Unknown", "Other"):
                continue
            counts[order] += 1
    return counts


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #
def order_colors(kingdoms: List[str]) -> "OrderedDict[str, str]":
    colors: "OrderedDict[str, str]" = OrderedDict()
    extra = iter(_EXTRA_HUES)
    for k in kingdoms:
        colors[k] = KINGDOM_COLORS.get(k, next(extra, "#999999"))
    return colors


def plot(counts: "OrderedDict[str, Counter]", outdir: str, prefix: str,
         formats: List[str], log: bool, title: Optional[str]) -> None:
    kingdoms = list(counts.keys())
    colors = order_colors(kingdoms)

    # Orders present anywhere, sorted by total count (magnitude form).
    totals: "Counter[str]" = Counter()
    for c in counts.values():
        totals.update(c)
    orders = [o for o, _ in totals.most_common()]
    if not orders:
        raise SystemExit("ERROR: no TE orders counted.")

    x = np.arange(len(orders))
    nk = len(kingdoms)
    band = 0.82
    w = band / nk
    grand_max = max(totals.values())

    fig_w = max(4.2, 0.55 * len(orders) + 1.2)
    fig, ax = plt.subplots(figsize=(fig_w, 4.0), constrained_layout=True)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="0.88", linewidth=0.6)
    ax.xaxis.grid(False)

    for i, k in enumerate(kingdoms):
        vals = np.array([counts[k].get(o, 0) for o in orders], dtype=float)
        offset = (i - (nk - 1) / 2.0) * w
        ax.bar(x + offset, vals, width=w * 0.9, color=colors[k],
               label="%s (n=%d)" % (k.capitalize(), int(vals.sum())), zorder=3)
        # Direct count labels in muted ink (not the bar hue), rotated to fit.
        for xx, v in zip(x + offset, vals):
            if v > 0:
                ax.text(xx, v * (1.08 if log else 1.0) + (0 if log else grand_max * 0.012),
                        f"{int(v):,}", ha="center", va="bottom",
                        rotation=90, fontsize=5.2, color="0.25")

    ax.set_xticks(x)
    ax.set_xticklabels(orders, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Number of TE families")
    ax.tick_params(axis="x", length=0)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.9, top=grand_max * 3.2)
    else:
        ax.set_ylim(0, grand_max * 1.20)
    ax.legend(loc="upper right", title="Kingdom", title_fontsize=7)
    if title:
        ax.set_title(title)

    tu.ensure_dir(outdir)
    for fmt in formats:
        p = os.path.join(outdir, f"{prefix}.{fmt}")
        fig.savefig(p, bbox_inches="tight")
        LOG.info("wrote %s", p)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Count TE families by order per kingdom and plot (Nature style).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = p.add_mutually_exclusive_group()
    src.add_argument("--truth", nargs="+", metavar="KINGDOM=FILE",
                     help="Explicit truth TSVs, e.g. animal=... plant=... fungi=...")
    src.add_argument("--indir", default=".",
                     help="Directory to auto-discover panteon_*.heldout.truth.tsv.")
    p.add_argument("--order-col", default="label",
                   help="Column to read from, RE-NORMALISED live (default 'label', "
                        "the raw label — so label fixes apply). Use 'label_order' "
                        "to trust the precomputed column instead.")
    p.add_argument("--level", default="order",
                   choices=["class", "order", "superfamily"],
                   help="Granularity to count at (default order).")
    p.add_argument("--keep-nonte", action="store_true",
                   help="Also count Unknown/Other (dropped by default).")
    p.add_argument("--log", action="store_true", help="Log-scale the count axis.")
    p.add_argument("--title", default=None, help="Optional figure title.")
    p.add_argument("--formats", nargs="+", default=["png", "pdf"])
    p.add_argument("--prefix", default="te_order_counts")
    p.add_argument("--outdir", default="order_count_figures")
    p.add_argument("--label-map", help="JSON of extra label synonyms.")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)

    sources = parse_sources(args.truth) if args.truth else discover(args.indir)
    if not sources:
        raise SystemExit("ERROR: no truth TSVs found (use --truth or --indir).")

    # Kingdoms in the canonical animal/plant/fungi order, then any extras.
    ordered_kingdoms = [k for k in KINGDOM_COLORS if k in sources] + \
                       [k for k in sources if k not in KINGDOM_COLORS]

    counts: "OrderedDict[str, Counter]" = OrderedDict()
    for k in ordered_kingdoms:
        counts[k] = count_orders(sources[k], args.order_col, args.level,
                                 args.keep_nonte)
        LOG.info("%-8s %s -> %d families across %d %ss",
                 k, sources[k], sum(counts[k].values()), len(counts[k]), args.level)

    # Tidy counts table (order x kingdom + totals), rows sorted by total.
    totals: "Counter[str]" = Counter()
    for c in counts.values():
        totals.update(c)
    all_orders = [o for o, _ in totals.most_common()]
    tu.ensure_dir(args.outdir)
    tsv = os.path.join(args.outdir, f"{args.prefix}.tsv")
    with open(tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([args.level] + [k for k in counts] + ["total"])
        for o in all_orders:
            rowvals = [counts[k].get(o, 0) for k in counts]
            w.writerow([o] + rowvals + [sum(rowvals)])
        w.writerow(["total"] + [sum(counts[k].values()) for k in counts]
                   + [sum(sum(c.values()) for c in counts.values())])
    LOG.info("wrote %s", tsv)

    plot(counts, args.outdir, args.prefix, args.formats, args.log, args.title)
    return 0


if __name__ == "__main__":
    sys.exit(main())
