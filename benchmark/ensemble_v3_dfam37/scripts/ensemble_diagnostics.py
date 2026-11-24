#!/usr/bin/env python3
"""
ensemble_diagnostics.py — decide whether to downweight, prune, or keep each model
in the TE-classification ensemble, from evidence rather than intuition.

A voting/weighted ensemble beats its best member ONLY when the weaker members are
individually decent AND make independent, complementary errors (right where the
strong models are wrong). If they are weak *and* their errors coincide with the
majority, they just drag the consensus off the strong models' correct calls. This
script measures exactly that, per kingdom/tag:

  1. Oracle ceiling      the accuracy an ideal per-sequence selector could reach
                         (fraction where ANY model is correct) — the hard upper
                         bound for any voting or stacking combiner.
  2. Per-model contribution
                         accuracy, coverage, macro-F1, and crucially
                         'unique_correct' (# sequences a model gets right that NO
                         other model does) and 'oracle_loss_if_removed'
                         (= unique_correct / N). A model with unique_correct ~ 0
                         adds no reachable signal -> a removal candidate.
  3. Error correlation   pairwise Pearson correlation of the models' error vectors.
                         High positive = redundant, coincident failures.
  4. Best-subset search  the model subset that maximises the (majority- or
                         weighted-) voted macro-F1. Exhaustive for a handful of
                         models, greedy forward otherwise. Tells you whether ANY
                         subset beats the best single model — and which models to
                         drop if so.

Inputs mirror evaluate_classification.py (so predictions/truth are read and
normalised identically). Ensemble columns in the voting table (Ensemble*) are
excluded automatically — only base models are analysed.

Example
-------
  python ensemble_diagnostics.py \
      --truth truth/panteon_fungi.heldout.truth.tsv \
      --voting-table results/panteon_fungi.voting_table.tsv \
      --level order --tag panteon_fungi --outdir diag_out
"""
from __future__ import annotations

import argparse
import csv
import itertools
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

import te_label_utils as tu
import evaluate_classification as ec

LOG = tu.setup_logger("ensemble_diagnostics")
ABSTAIN = ec._ABSTAIN


# --------------------------------------------------------------------------- #
# Voting + scoring for a model subset
# --------------------------------------------------------------------------- #
def ensemble_vote(subset: List[str], pred_norm: Dict[str, Dict[str, str]],
                  sid: str, weights: Optional[Dict[str, float]]) -> str:
    """Majority (or weighted) vote of a subset on one sequence; tie -> Unknown."""
    tally: Dict[str, float] = defaultdict(float)
    for m in subset:
        lab = pred_norm[m].get(sid)
        if not lab or lab.lower() in ABSTAIN:
            continue
        tally[lab] += weights.get(m.lower(), 1.0) if weights else 1.0
    if not tally:
        return "Unknown"
    top = max(tally.values())
    winners = [c for c, v in tally.items() if abs(v - top) < 1e-9]
    return winners[0] if len(winners) == 1 else "Unknown"


def subset_metrics(subset: List[str], truth_norm: Dict[str, str],
                   pred_norm: Dict[str, Dict[str, str]], ids: List[str],
                   weights: Optional[Dict[str, float]]) -> dict:
    """Confusion-matrix metrics for a subset's voted predictions."""
    conf: Dict[Tuple[str, str], int] = Counter()
    classes = set()
    for sid in ids:
        t = truth_norm[sid]
        p = ensemble_vote(subset, pred_norm, sid, weights)
        conf[(t, p)] += 1
        classes.add(t)
        classes.add(p)
    summary, _ = ec.metrics_from_confusion(conf, ec.order_classes(classes))
    return summary


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def per_model_contribution(models: List[str], truth_norm: Dict[str, str],
                           pred_norm: Dict[str, Dict[str, str]], ids: List[str]
                           ) -> Tuple[dict, float]:
    """Per-model correctness stats + the oracle accuracy."""
    correct_total: "Counter[str]" = Counter()
    unique: "Counter[str]" = Counter()
    redundant: "Counter[str]" = Counter()
    attempted: "Counter[str]" = Counter()
    oracle_correct = 0
    n = len(ids)

    for sid in ids:
        t = truth_norm[sid]
        hits = []
        for m in models:
            p = pred_norm[m].get(sid)
            if p and p.lower() not in ABSTAIN:
                attempted[m] += 1
            if p == t:
                hits.append(m)
        if hits:
            oracle_correct += 1
        if len(hits) == 1:
            unique[hits[0]] += 1
        elif len(hits) > 1:
            for m in hits:
                redundant[m] += 1
        for m in hits:
            correct_total[m] += 1

    stats = {}
    for m in models:
        stats[m] = {
            "accuracy": correct_total[m] / n if n else 0.0,
            "coverage": attempted[m] / n if n else 0.0,
            "correct_total": correct_total[m],
            "unique_correct": unique[m],
            "redundant_correct": redundant[m],
            "oracle_loss_if_removed": unique[m] / n if n else 0.0,
        }
    return stats, (oracle_correct / n if n else 0.0)


def error_correlation(models: List[str], truth_norm: Dict[str, str],
                      pred_norm: Dict[str, Dict[str, str]], ids: List[str]
                      ) -> np.ndarray:
    """Pearson correlation of the models' binary error (1=wrong/abstain) vectors."""
    err = np.zeros((len(models), len(ids)), dtype=float)
    for j, sid in enumerate(ids):
        t = truth_norm[sid]
        for i, m in enumerate(models):
            p = pred_norm[m].get(sid)
            err[i, j] = 0.0 if p == t else 1.0
    with np.errstate(invalid="ignore"):
        corr = np.corrcoef(err)          # constant rows -> nan (handled on write)
    return np.atleast_2d(corr)


def search_subsets(models: List[str], truth_norm, pred_norm, ids,
                   weights, max_exhaustive: int
                   ) -> Tuple[List[Tuple[Tuple[str, ...], dict]], str]:
    """Score model subsets; exhaustive if small, else greedy forward."""
    scored: List[Tuple[Tuple[str, ...], dict]] = []
    if len(models) <= max_exhaustive:
        mode = "exhaustive"
        for r in range(1, len(models) + 1):
            for combo in itertools.combinations(models, r):
                scored.append((combo, subset_metrics(list(combo), truth_norm,
                                                      pred_norm, ids, weights)))
    else:
        mode = "greedy-forward"
        chosen: List[str] = []
        remaining = list(models)
        best_f1 = -1.0
        while remaining:
            trials = []
            for m in remaining:
                s = subset_metrics(chosen + [m], truth_norm, pred_norm, ids, weights)
                trials.append((m, s))
            m_best, s_best = max(trials, key=lambda x: x[1]["macro_f1"])
            chosen.append(m_best)
            remaining.remove(m_best)
            scored.append((tuple(chosen), s_best))
            if s_best["macro_f1"] <= best_f1 + 1e-6:
                # keep growing for the record, but note the plateau
                pass
            best_f1 = max(best_f1, s_best["macro_f1"])
    scored.sort(key=lambda x: x[1]["macro_f1"], reverse=True)
    return scored, mode


# --------------------------------------------------------------------------- #
# Loading (mirrors evaluate_classification)
# --------------------------------------------------------------------------- #
def load_models(args) -> Dict[str, Dict[str, str]]:
    model_preds: Dict[str, Dict[str, str]] = {}
    if args.voting_table:
        models, per_model = ec.load_voting_table(args.voting_table)
        for m in models:
            model_preds[m] = per_model[m]
    for item in args.predictions:
        if "=" not in item:
            raise SystemExit(f"ERROR: --predictions entry '{item}' is not NAME=FILE")
        name, path = item.split("=", 1)
        model_preds[name] = ec.load_prediction_source(
            path, args.pred_format, args.gff_id_attr, args.gff_label_attr,
            args.pred_id_col, args.pred_label_col)
    # Drop derived ensemble columns and any user-excluded models.
    drop = {x.lower() for x in args.exclude_models}
    for m in list(model_preds):
        if m.lower().startswith("ensemble") or m.lower() in drop:
            model_preds.pop(m)
            LOG.info("excluding non-base column '%s'", m)
    return model_preds


def parse_weights(spec: List[str]) -> Optional[Dict[str, float]]:
    """Parse weights given space- OR comma-separated (keys lower-cased)."""
    if not spec:
        return None
    out: Dict[str, float] = {}
    for item in spec:                        # each item may itself be 'a=1,b=2'
        for tok in item.replace(";", ",").split(","):
            tok = tok.strip()
            if not tok:
                continue
            if "=" not in tok:
                LOG.warning("ignoring bad weight '%s' (expected NAME=W)", tok)
                continue
            k, v = tok.split("=", 1)
            try:
                out[k.strip().lower()] = float(v)
            except ValueError:
                LOG.warning("ignoring bad weight '%s'", tok)
    return out or None


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
def _w(path: str, header: List[str], rows: List[list]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Oracle / complementarity / best-subset diagnostics for the "
                    "TE-classification ensemble.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--truth", required=True, help="Ground-truth TSV.")
    p.add_argument("--truth-id-col", default="family_id")
    p.add_argument("--truth-label-col", default="label")
    p.add_argument("--voting-table", help="Wide voting table (id + per-model cols).")
    p.add_argument("--predictions", nargs="+", default=[], metavar="NAME=FILE")
    p.add_argument("--pred-format", default="auto",
                   choices=["auto", "tsv", "fasta", "gff3"])
    p.add_argument("--pred-id-col", default=None)
    p.add_argument("--pred-label-col", default=None)
    p.add_argument("--gff-id-attr", default="ID")
    p.add_argument("--gff-label-attr", default="classification")
    p.add_argument("--level", default="order", choices=["class", "order", "superfamily"])
    p.add_argument("--exclude-models", nargs="*", default=[],
                   help="Base model columns to ignore.")
    p.add_argument("--weights", nargs="*", default=[], metavar="NAME=W",
                   help="If set, subset search uses WEIGHTED voting with these "
                        "weights; otherwise plain majority (isolates subset choice "
                        "from weight tuning).")
    p.add_argument("--keep-list", default=None,
                   help="Restrict to these ids (leakage-free set).")
    p.add_argument("--exclude-list", default=None)
    p.add_argument("--keep-nonte-truth", action="store_true")
    p.add_argument("--max-exhaustive", type=int, default=14,
                   help="Exhaustive subset search up to this many models; greedy above.")
    p.add_argument("--top", type=int, default=50,
                   help="How many top subsets to write.")
    p.add_argument("--tag", default=None, help="Output filename tag.")
    p.add_argument("--label-map", help="JSON of extra label synonyms.")
    p.add_argument("--outdir", default="diag_out")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)
    tu.ensure_dir(args.outdir)
    tag = args.tag or os.path.splitext(os.path.basename(args.truth))[0]

    keep_ids = ec.load_id_list(args.keep_list)
    exclude_ids = ec.load_id_list(args.exclude_list)
    truth_raw = ec.load_truth(args.truth, args.truth_id_col, args.truth_label_col,
                              args.level, args.keep_nonte_truth, keep_ids, exclude_ids)
    if not truth_raw:
        raise SystemExit("ERROR: no usable ground-truth rows.")

    model_preds = load_models(args)
    if not model_preds:
        raise SystemExit("ERROR: no base-model predictions (use --voting-table / --predictions).")
    models = ec.sort_models(list(model_preds)) if hasattr(ec, "sort_models") \
        else sorted(model_preds)
    weights = parse_weights(args.weights)

    # Normalise once.
    truth_norm = {sid: tu.normalize_label(lab, args.level) for sid, lab in truth_raw.items()}
    pred_norm = {m: {sid: tu.normalize_label(lab, args.level)
                     for sid, lab in preds.items()} for m, preds in model_preds.items()}
    ids = sorted(truth_norm)
    n = len(ids)
    LOG.info("tag '%s': %d sequences, %d base models: %s",
             tag, n, len(models), ", ".join(models))

    # 1-2. Oracle + per-model contribution.
    stats, oracle = per_model_contribution(models, truth_norm, pred_norm, ids)
    single = {m: subset_metrics([m], truth_norm, pred_norm, ids, None) for m in models}
    best_single = max(models, key=lambda m: single[m]["macro_f1"])

    contrib_rows = []
    for m in sorted(models, key=lambda x: single[x]["macro_f1"], reverse=True):
        s = stats[m]
        contrib_rows.append([
            m, f"{single[m]['macro_f1']:.4f}", f"{single[m]['mcc']:.4f}",
            f"{s['accuracy']:.4f}", f"{s['coverage']:.4f}", s["correct_total"],
            s["unique_correct"], s["redundant_correct"],
            f"{s['oracle_loss_if_removed']:.4f}"])
    _w(os.path.join(args.outdir, f"{tag}.per_model_contribution.tsv"),
       ["model", "macro_f1", "mcc", "accuracy", "coverage", "correct_total",
        "unique_correct", "redundant_correct", "oracle_loss_if_removed"], contrib_rows)

    # 3. Error correlation.
    corr = error_correlation(models, truth_norm, pred_norm, ids)
    corr_rows = []
    for i, m in enumerate(models):
        row = [m]
        for j in range(len(models)):
            v = corr[i, j]
            row.append("NA" if not np.isfinite(v) else f"{v:.3f}")
        corr_rows.append(row)
    _w(os.path.join(args.outdir, f"{tag}.error_correlation.tsv"),
       ["model"] + models, corr_rows)

    # 4. Best-subset search.
    scored, mode = search_subsets(models, truth_norm, pred_norm, ids, weights,
                                  args.max_exhaustive)
    full = subset_metrics(models, truth_norm, pred_norm, ids, weights)
    best_combo, best_s = scored[0]
    sub_rows = []
    for combo, s in scored[:args.top]:
        sub_rows.append(["+".join(combo), len(combo), f"{s['macro_f1']:.4f}",
                         f"{s['mcc']:.4f}", f"{s['accuracy']:.4f}",
                         f"{s['balanced_accuracy']:.4f}"])
    _w(os.path.join(args.outdir, f"{tag}.subset_ranking.tsv"),
       ["subset", "n_models", "macro_f1", "mcc", "accuracy", "balanced_accuracy"],
       sub_rows)

    # Headline summary.
    vote_mode = "weighted" if weights else "majority"
    beats = best_s["macro_f1"] > single[best_single]["macro_f1"] + 1e-6
    summ = [
        ["n_sequences", n],
        ["level", args.level],
        ["vote_mode_for_search", vote_mode],
        ["subset_search", mode],
        ["oracle_accuracy", f"{oracle:.4f}"],
        ["best_single_model", best_single],
        ["best_single_macro_f1", f"{single[best_single]['macro_f1']:.4f}"],
        ["best_single_mcc", f"{single[best_single]['mcc']:.4f}"],
        ["full_ensemble_macro_f1", f"{full['macro_f1']:.4f}"],
        ["full_ensemble_mcc", f"{full['mcc']:.4f}"],
        ["best_subset", "+".join(best_combo)],
        ["best_subset_macro_f1", f"{best_s['macro_f1']:.4f}"],
        ["best_subset_mcc", f"{best_s['mcc']:.4f}"],
        ["best_subset_beats_best_single", str(beats)],
        ["oracle_headroom_over_best_single",
         f"{oracle - single[best_single]['accuracy']:.4f}"],
    ]
    _w(os.path.join(args.outdir, f"{tag}.diagnostics_summary.tsv"),
       ["metric", "value"], summ)

    # Console verdict.
    LOG.info("oracle=%.3f | best single=%s (F1 %.3f) | full ensemble F1 %.3f | "
             "best subset=%s (F1 %.3f) | beats best single: %s",
             oracle, best_single, single[best_single]["macro_f1"],
             full["macro_f1"], "+".join(best_combo), best_s["macro_f1"], beats)
    dead = [m for m in models if stats[m]["unique_correct"] == 0]
    if dead:
        LOG.info("removal candidates (unique_correct==0, no oracle loss): %s",
                 ", ".join(dead))
    LOG.info("Wrote diagnostics to %s/ (%s.*.tsv)", args.outdir, tag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
