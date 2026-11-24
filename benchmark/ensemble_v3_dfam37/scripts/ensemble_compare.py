#!/usr/bin/env python3
"""
ensemble_compare.py — compare ensemble strategies for TE classification and rank
all methods with a Friedman + Nemenyi test (à la PanTEon, Orozco-Arias et al. 2026).

Your production pipeline uses hard majority voting. PanTEon reports that a
supervised STACKING meta-classifier beats simple/weighted voting by up to ~13% F1
and the best single model by ~7-8%. This script evaluates, on the *same* per-model
predictions you already produced (no re-running the classifiers):

  * each individual model
  * simple majority vote          (optionally requiring --min-votes agreement)
  * weighted majority vote        (PanTEon-style weights, override with --weights)
  * XGBoost stacking meta-model   (leakage-safe: stratified out-of-fold predictions;
                                   falls back to sklearn HistGradientBoosting if
                                   xgboost is not installed)

and ranks every method with the Friedman test + Nemenyi critical-difference over
per-superfamily F1 — the standard multi-classifier comparison.

Inputs mirror evaluate_classification.py (truth TSV + --voting-table and/or
--predictions), so the same --keep-list / --exclude-list leakage control applies.

Outputs (to --outdir):
  <tag>.ensemble_comparison.tsv    method -> accuracy, macro_f1, mcc, coverage
  <tag>.stacking_oof_predictions.tsv
  <tag>.method_ranks.tsv           mean Friedman rank per method
  <tag>.friedman_nemenyi.tsv       chi2, p-value, critical difference, N blocks

Example
-------
  python ensemble_compare.py --truth dfam_fungi.tsv --voting-table voting_table.tsv \
      --level superfamily --tag dfam_fungi --keep-list heldout.ids --outdir cmp_out
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np

import te_label_utils as tu
# Reuse the evaluator's loaders + confusion-matrix metrics (same directory).
from evaluate_classification import (load_truth, load_id_list, load_voting_table,
                                     load_prediction_source, metrics_from_confusion,
                                     order_classes)

LOG = tu.setup_logger("ensemble_compare")

_ABSTAIN = "ABSTAIN"

# PanTEon (2026) tool-specific weights for weighted voting; others default to 1.0.
DEFAULT_WEIGHTS = {
    "neuralte": 1.5, "terrier": 1.5, "create": 1.3, "teclass2": 1.2,
    "deepte": 1.1, "classifyte": 1.0, "terl": 1.0,
}

# Nemenyi critical values q_alpha (alpha=0.05), indexed by number of methods k
# (Demsar 2006, Table 5 — studentized range / sqrt(2)).
_NEMENYI_Q05 = {
    2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728, 6: 2.850, 7: 2.949, 8: 3.031,
    9: 3.102, 10: 3.164, 11: 3.219, 12: 3.268, 13: 3.313, 14: 3.354,
    15: 3.391, 16: 3.426, 17: 3.458, 18: 3.489, 19: 3.517, 20: 3.544,
}


# --------------------------------------------------------------------------- #
# Prediction assembly
# --------------------------------------------------------------------------- #
def norm_or_abstain(raw: object, level: str) -> str:
    lab = tu.normalize_label(raw, level)
    return _ABSTAIN if lab in ("Unknown", "Other") else lab


def assemble(truth: Dict[str, str], model_preds: Dict[str, Dict[str, str]],
             level: str) -> Tuple[List[str], List[str], List[str], np.ndarray]:
    """Return (ids, y_true, models, pred_matrix[N,M] of normalised labels)."""
    ids = sorted(truth)
    models = list(model_preds)
    y_true = [tu.normalize_label(truth[i], level) for i in ids]
    mat = np.empty((len(ids), len(models)), dtype=object)
    for j, m in enumerate(models):
        pm = model_preds[m]
        for i, sid in enumerate(ids):
            mat[i, j] = norm_or_abstain(pm.get(sid), level)
    return ids, y_true, models, mat


# --------------------------------------------------------------------------- #
# Ensemble strategies -> per-sample predicted label
# --------------------------------------------------------------------------- #
def vote_majority(row: np.ndarray, min_votes: int) -> str:
    votes = [v for v in row if v != _ABSTAIN]
    if not votes:
        return "Unknown"
    counts = Counter(votes)
    top, n = counts.most_common(1)[0]
    if n < max(min_votes, 1):
        return "Unknown"
    # tie at the top -> abstain
    if sum(1 for c in counts.values() if c == n) > 1:
        return "Unknown"
    return top


def vote_weighted(row: np.ndarray, weights: np.ndarray) -> str:
    acc: Dict[str, float] = {}
    for v, w in zip(row, weights):
        if v != _ABSTAIN:
            acc[v] = acc.get(v, 0.0) + w
    if not acc:
        return "Unknown"
    best = max(acc.values())
    winners = [k for k, s in acc.items() if abs(s - best) < 1e-9]
    return winners[0] if len(winners) == 1 else "Unknown"


def stacking_oof(mat: np.ndarray, y_true: List[str], models: List[str],
                 folds: int, seed: int) -> Optional[List[str]]:
    """Leakage-safe stacking: out-of-fold predictions from a gradient-boosted
    meta-classifier over one-hot-encoded per-model predictions."""
    from sklearn.preprocessing import OneHotEncoder, LabelEncoder
    from sklearn.model_selection import StratifiedKFold
    try:
        from xgboost import XGBClassifier
        def make_clf():
            # Multiclass objective is inferred from y by the sklearn wrapper.
            return XGBClassifier(
                n_estimators=300, max_depth=6, learning_rate=0.1,
                subsample=0.9, colsample_bytree=0.9, tree_method="hist",
                random_state=seed, verbosity=0)
        backend = "xgboost"
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        def make_clf():
            return HistGradientBoostingClassifier(random_state=seed)
        backend = "sklearn-HGB"
    LOG.info("Stacking meta-classifier backend: %s", backend)

    y = np.array(y_true)
    # StratifiedKFold needs >= folds samples per class; drop rarer classes.
    counts = Counter(y_true)
    eff_folds = max(2, min(folds, min(counts.values())))
    keep = np.array([counts[c] >= eff_folds for c in y])
    if keep.sum() < eff_folds or len(set(y[keep])) < 2:
        LOG.warning("Not enough per-class support for stacking (need >=%d/class); "
                    "skipping.", eff_folds)
        return None
    if eff_folds != folds:
        LOG.warning("Reduced stacking folds %d -> %d (rare-class support).",
                    folds, eff_folds)

    Xraw = mat[keep]
    yk = y[keep]
    # sparse_output (sklearn >=1.2) vs sparse (older) — support both.
    try:
        enc = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        enc = OneHotEncoder(handle_unknown="ignore", sparse=False)
    X = enc.fit_transform(Xraw)
    le = LabelEncoder()
    yi = le.fit_transform(yk)

    oof = np.empty(len(yi), dtype=int)
    skf = StratifiedKFold(n_splits=eff_folds, shuffle=True, random_state=seed)
    for tr, te in skf.split(X, yi):
        clf = make_clf()
        clf.fit(X[tr], yi[tr])
        oof[te] = clf.predict(X[te]).astype(int)
    oof_labels = le.inverse_transform(oof)

    # Re-expand to full length; dropped rare-class samples abstain.
    out = ["Unknown"] * len(y_true)
    ptr = 0
    for i in range(len(y_true)):
        if keep[i]:
            out[i] = oof_labels[ptr]
            ptr += 1
    return out


# --------------------------------------------------------------------------- #
# Scoring + statistics
# --------------------------------------------------------------------------- #
def score(y_true: List[str], y_pred: List[str]) -> Tuple[dict, Dict[str, float]]:
    """Return (summary metrics, per-class F1 dict) for one method."""
    conf: Counter = Counter()
    cs = set()
    attempted = 0
    for t, p in zip(y_true, y_pred):
        pp = p if p not in (_ABSTAIN,) else "Unknown"
        if pp != "Unknown":
            attempted += 1
        conf[(t, pp)] += 1
        cs.add(t); cs.add(pp)
    classes = order_classes(cs)
    summary, per_class = metrics_from_confusion(conf, classes)
    summary["coverage"] = attempted / len(y_true) if y_true else 0.0
    f1 = {r["class"]: r["f1"] for r in per_class if r["support"] > 0}
    return summary, f1


def friedman_nemenyi(method_f1: Dict[str, Dict[str, float]]
                     ) -> Tuple[dict, Dict[str, float]]:
    """Friedman test + Nemenyi CD over per-class F1 (blocks = classes)."""
    from scipy.stats import friedmanchisquare
    methods = list(method_f1)
    classes = sorted(set().union(*[set(d) for d in method_f1.values()]))
    if len(methods) < 3 or len(classes) < 3:
        return ({"note": "need >=3 methods and >=3 classes for Friedman"}, {})
    # matrix methods x classes
    M = np.array([[method_f1[m].get(c, 0.0) for c in classes] for m in methods])
    chi2, p = friedmanchisquare(*[M[i, :] for i in range(len(methods))])
    # per-class ranks (higher F1 = rank 1)
    ranks = np.zeros_like(M)
    for j in range(M.shape[1]):
        order = (-M[:, j]).argsort()
        r = np.empty(len(methods)); r[order] = np.arange(1, len(methods) + 1)
        # average ties
        col = M[:, j]
        for val in set(col):
            idx = np.where(col == val)[0]
            if len(idx) > 1:
                r[idx] = r[idx].mean()
        ranks[:, j] = r
    mean_ranks = {m: float(ranks[i, :].mean()) for i, m in enumerate(methods)}
    k, N = len(methods), len(classes)
    q = _NEMENYI_Q05.get(k, _NEMENYI_Q05[max(_NEMENYI_Q05)])
    cd = q * math.sqrt(k * (k + 1) / (6.0 * N))
    stats = {"friedman_chi2": float(chi2), "p_value": float(p),
             "n_methods": k, "n_blocks_classes": N, "critical_difference": cd}
    return stats, mean_ranks


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compare majority / weighted / stacking ensembles and rank all "
                    "methods (Friedman + Nemenyi).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--truth", required=True)
    p.add_argument("--truth-id-col", default="family_id")
    p.add_argument("--truth-label-col", default="label")
    p.add_argument("--voting-table")
    p.add_argument("--predictions", nargs="+", default=[], metavar="NAME=FILE")
    p.add_argument("--pred-format", default="auto",
                   choices=["auto", "tsv", "fasta", "gff3"])
    p.add_argument("--pred-id-col", default=None)
    p.add_argument("--pred-label-col", default=None)
    p.add_argument("--gff-id-attr", default="ID")
    p.add_argument("--gff-label-attr", default="classification")
    p.add_argument("--level", default="superfamily",
                   choices=["class", "order", "superfamily"])
    p.add_argument("--min-votes", type=int, default=1,
                   help="Min agreeing votes for the simple-majority winner.")
    p.add_argument("--weights", nargs="*", default=[], metavar="NAME=W",
                   help="Override weighted-vote weights (default: PanTEon values).")
    p.add_argument("--stacking-folds", type=int, default=5,
                   help="Stratified folds for leakage-safe out-of-fold stacking.")
    p.add_argument("--keep-list", default=None)
    p.add_argument("--exclude-list", default=None)
    p.add_argument("--keep-nonte-truth", action="store_true")
    p.add_argument("--label-map", default=None)
    p.add_argument("--tag", default=None)
    p.add_argument("--outdir", default="ensemble_cmp")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--no-stacking", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)
    tu.ensure_dir(args.outdir)
    tag = args.tag or os.path.splitext(os.path.basename(args.truth))[0]

    keep_ids = load_id_list(args.keep_list)
    exclude_ids = load_id_list(args.exclude_list)
    truth = load_truth(args.truth, args.truth_id_col, args.truth_label_col,
                       args.level, args.keep_nonte_truth, keep_ids, exclude_ids)
    if not truth:
        raise SystemExit("ERROR: no usable ground-truth rows.")
    LOG.info("Loaded %d ground-truth families (level=%s)", len(truth), args.level)

    model_preds: Dict[str, Dict[str, str]] = {}
    if args.voting_table:
        models, per_model = load_voting_table(args.voting_table)
        for m in models:
            if m.lower() != "ensemble":       # drop any precomputed ensemble column
                model_preds[m] = per_model[m]
    for item in args.predictions:
        if "=" not in item:
            raise SystemExit(f"ERROR: --predictions entry '{item}' is not NAME=FILE")
        name, path = item.split("=", 1)
        model_preds[name] = load_prediction_source(
            path, args.pred_format, args.gff_id_attr, args.gff_label_attr,
            args.pred_id_col, args.pred_label_col)
    if len(model_preds) < 2:
        raise SystemExit("ERROR: need >=2 model prediction sources for an ensemble.")
    LOG.info("Models: %s", ", ".join(model_preds))

    ids, y_true, models, mat = assemble(truth, model_preds, args.level)

    # Weights vector (normalised model name -> weight).
    wmap = dict(DEFAULT_WEIGHTS)
    for kv in args.weights:
        if "=" in kv:
            k, v = kv.split("=", 1)
            wmap[k.lower()] = float(v)
    weights = np.array([wmap.get(m.lower(), 1.0) for m in models])
    LOG.info("Weighted-vote weights: %s",
             ", ".join(f"{m}={w:g}" for m, w in zip(models, weights)))

    # ── Build each method's predictions ────────────────────────────────────
    methods: Dict[str, List[str]] = {}
    for j, m in enumerate(models):
        methods[m] = [x if x != _ABSTAIN else "Unknown" for x in mat[:, j]]
    methods["majority_vote"] = [vote_majority(mat[i, :], args.min_votes)
                                for i in range(len(ids))]
    methods["weighted_vote"] = [vote_weighted(mat[i, :], weights)
                                for i in range(len(ids))]
    if not args.no_stacking:
        stk = stacking_oof(mat, y_true, models, args.stacking_folds, args.seed)
        if stk is not None:
            methods["stacking_xgb"] = stk
            with open(os.path.join(args.outdir, f"{tag}.stacking_oof_predictions.tsv"),
                      "w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["family_id", "true_label", "stacking_pred"])
                for sid, t, p in zip(ids, y_true, stk):
                    w.writerow([sid, t, p])

    # ── Score every method ─────────────────────────────────────────────────
    rows = []
    method_f1: Dict[str, Dict[str, float]] = {}
    for name, preds in methods.items():
        summ, f1 = score(y_true, preds)
        method_f1[name] = f1
        rows.append((name, summ))
        LOG.info("%-16s acc=%.3f macroF1=%.3f MCC=%.3f cov=%.3f",
                 name, summ["accuracy"], summ["macro_f1"], summ["mcc"],
                 summ["coverage"])

    is_ens = {"majority_vote", "weighted_vote", "stacking_xgb"}
    with open(os.path.join(args.outdir, f"{tag}.ensemble_comparison.tsv"),
              "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["method", "kind", "n", "accuracy", "balanced_accuracy",
                    "macro_f1", "weighted_f1", "mcc", "cohen_kappa", "coverage"])
        for name, s in sorted(rows, key=lambda r: r[1]["macro_f1"], reverse=True):
            w.writerow([name, "ensemble" if name in is_ens else "model",
                        s["n"], f"{s['accuracy']:.6f}",
                        f"{s['balanced_accuracy']:.6f}", f"{s['macro_f1']:.6f}",
                        f"{s['weighted_f1']:.6f}", f"{s['mcc']:.6f}",
                        f"{s['cohen_kappa']:.6f}", f"{s['coverage']:.6f}"])

    # ── Friedman + Nemenyi ranking ─────────────────────────────────────────
    stats, mean_ranks = friedman_nemenyi(method_f1)
    with open(os.path.join(args.outdir, f"{tag}.friedman_nemenyi.tsv"),
              "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in stats.items():
            w.writerow([k, f"{v:.6f}" if isinstance(v, float) else v])
    if mean_ranks:
        with open(os.path.join(args.outdir, f"{tag}.method_ranks.tsv"),
                  "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["method", "mean_rank", "kind"])
            for m, r in sorted(mean_ranks.items(), key=lambda x: x[1]):
                w.writerow([m, f"{r:.4f}", "ensemble" if m in is_ens else "model"])
        LOG.info("Friedman chi2=%.3f p=%.3g | CD(0.05)=%.3f | best rank: %s",
                 stats.get("friedman_chi2", float("nan")),
                 stats.get("p_value", float("nan")),
                 stats.get("critical_difference", float("nan")),
                 min(mean_ranks, key=mean_ranks.get))

    LOG.info("Wrote comparison + ranks to %s/", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
