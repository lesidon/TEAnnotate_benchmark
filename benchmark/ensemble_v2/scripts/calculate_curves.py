#!/usr/bin/env python3
"""
calculate_curves.py — ROC and Precision-Recall curves (per-class one-vs-rest,
plus micro- and macro-averages) for a TE classifier or the ensemble vote.

Because TE superfamily distributions are highly imbalanced, PR curves / Average
Precision (AP) are the PRIMARY summary; ROC/AUC are reported alongside for
completeness (they are optimistic under imbalance).

Three input modes:

  --proba FILE   Wide TSV with a soft-score column per class. Columns:
                   <id>, <true_label>, and one probability column per class,
                   prefixed (default "prob:"), e.g. "prob:LTR/Gypsy".
                 Use for any model that emits class probabilities.

  --preds FILE   Per-model HARD-prediction TSV (id + predicted label), plus
                 --truth for the labels. Most classifiers emit only a class, no
                 probabilities: each becomes a single operating point in ROC/PR
                 space (predicted class scored 1.0, others 0; abstentions all-0).
                 The per-class "AUC" of that degenerate curve equals the model's
                 balanced accuracy for that class — a fair point-estimate to plot
                 several label-only models on the same axes.

  --votes FILE   Wide hard-label table (id + one column per constituent model),
                 e.g. the pipeline's voting_table.tsv, plus --truth for labels.
                 --model selects what to score:
                   --model Ensemble  -> vote-fraction curve (uses ALL columns):
                       s_c = (# models predicting c) / (denominator),
                       the quantity the 3-of-5 rule thresholds.
                   --model <column>  -> that single model's hard-prediction
                       curve (one operating point), read straight from its
                       column — no separate per-model file needed.

Outputs (per --model / --tag) to --outdir:
  <tag>.<model>.roc_coords.tsv   class, fpr, tpr, threshold
  <tag>.<model>.pr_coords.tsv    class, recall, precision, threshold   (PRIMARY)
  <tag>.<model>.curve_auc.tsv    class, n_pos, roc_auc, pr_auc(AP)

Labels are normalised (te_label_utils) at the requested --level so probability
columns / votes for synonymous classes are merged consistently.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

import te_label_utils as tu

LOG = tu.setup_logger("calculate_curves")

_ABSTAIN = {"unknown", "other", "na", "n/a", "none", "", "."}


# --------------------------------------------------------------------------- #
# Input assembly  ->  (classes, y_true[str], scores[N,K])
# --------------------------------------------------------------------------- #
def load_proba(path: str, level: str, id_col: Optional[str], true_col: Optional[str],
               prob_prefix: str) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """Assemble a normalised score matrix from a wide probability TSV (chunked)."""
    ids: List[str] = []
    y_true: List[str] = []
    class_cols: Dict[str, List[str]] = {}   # normalised class -> raw prob columns
    norm_classes: List[str] = []
    score_rows: List[np.ndarray] = []
    cid = ctrue = None
    prob_cols: List[str] = []

    for chunk in tu.iter_tsv_chunks(path, dtype=None):
        if cid is None:
            cols = list(chunk.columns)
            cid = id_col or tu.detect_column(cols, ["sample_id", "family_id", "id"])
            ctrue = true_col or tu.detect_column(
                cols, ["true_label", "truth", "label", "class", "y_true"])
            prob_cols = [c for c in cols if str(c).startswith(prob_prefix)]
            if cid is None or ctrue is None or not prob_cols:
                raise SystemExit(
                    f"ERROR: need id, true-label and '{prob_prefix}*' columns. "
                    f"Found: {cols}")
            for c in prob_cols:
                nc = tu.normalize_label(c[len(prob_prefix):], level)
                class_cols.setdefault(nc, []).append(c)
            norm_classes = list(class_cols.keys())
            LOG.info("Probability columns collapsed to %d classes: %s",
                     len(norm_classes), ", ".join(norm_classes))
        # Sum synonymous probability columns into each normalised class.
        mat = np.zeros((len(chunk), len(norm_classes)), dtype=np.float64)
        for j, nc in enumerate(norm_classes):
            for rawc in class_cols[nc]:
                mat[:, j] += chunk[rawc].astype(float).to_numpy()
        score_rows.append(mat)
        ids.extend(chunk[cid].astype(str).tolist())
        y_true.extend(tu.normalize_label(v, level) for v in chunk[ctrue].astype(str))

    scores = np.vstack(score_rows) if score_rows else np.zeros((0, len(norm_classes)))
    # Row-normalise so each row sums to 1 where possible (defensive).
    row_sums = scores.sum(axis=1, keepdims=True)
    np.divide(scores, row_sums, out=scores, where=row_sums > 0)
    return norm_classes, np.asarray(y_true), scores


def load_votes(path: str, truth_path: str, level: str,
               id_col: Optional[str], truth_id_col: str, truth_label_col: str,
               exclude_cols: List[str], denominator: str
               ) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """Build ensemble vote-fraction scores from a wide hard-label voting table."""
    # 1. Truth map.
    truth = _load_truth_map(truth_path, level, truth_id_col, truth_label_col)

    # 2. Per-id vote counts.
    exclude = {c.lower() for c in exclude_cols}
    votes: Dict[str, Dict[str, int]] = {}
    cast: Dict[str, int] = {}
    model_cols: List[str] = []
    vid = None
    for chunk in tu.iter_tsv_chunks(path):
        if vid is None:
            cols = list(chunk.columns)
            vid = id_col or tu.detect_column(cols, ["repeat_id", "sample_id", "id"])
            if vid is None:
                vid = cols[0]
            model_cols = [c for c in cols
                          if c != vid and c.lower() not in exclude
                          and not c.lower().startswith("ensemble")]
            LOG.info("Vote model columns: %s", ", ".join(model_cols))
        for _, row in chunk.iterrows():
            sid = str(row[vid]).split("#", 1)[0]
            if sid not in truth:
                continue
            d = votes.setdefault(sid, {})
            c = cast.get(sid, 0)
            for m in model_cols:
                lab = tu.normalize_label(str(row[m]), level)
                if lab.lower() in _ABSTAIN:
                    continue
                d[lab] = d.get(lab, 0) + 1
                c += 1
            cast[sid] = c

    ids = [s for s in votes if s in truth]
    class_set = set()
    for d in votes.values():
        class_set.update(d.keys())
    class_set.update(truth[s] for s in ids)
    class_set = {c for c in class_set if c.lower() not in _ABSTAIN}
    preferred = tu.CANONICAL_CLASSES + tu.CANONICAL_ORDERS
    classes = [c for c in preferred if c in class_set] + \
              sorted(c for c in class_set if c not in preferred)

    n_models = len(model_cols)
    y_true = np.array([truth[s] for s in ids])
    scores = np.zeros((len(ids), len(classes)), dtype=np.float64)
    cidx = {c: j for j, c in enumerate(classes)}
    for i, s in enumerate(ids):
        denom = n_models if denominator == "fixed" else max(cast.get(s, 0), 1)
        for c, n in votes[s].items():
            if c in cidx:
                scores[i, cidx[c]] = n / denom
    return classes, y_true, scores


def _load_truth_map(truth_path: str, level: str,
                    truth_id_col: str, truth_label_col: str) -> Dict[str, str]:
    """Read {id: normalised_label} from a ground-truth TSV (streaming)."""
    truth: Dict[str, str] = {}
    for chunk in tu.iter_tsv_chunks(truth_path):
        cols = list(chunk.columns)
        cid = tu.detect_column(cols, [truth_id_col, "family_id", "id"])
        clab = tu.detect_column(cols, [truth_label_col, "label", "label_class",
                                        "label_order", "rm_type", "class"])
        if cid is None or clab is None:
            raise SystemExit(f"ERROR: truth id/label columns not in {cols}")
        for sid, lab in zip(chunk[cid].astype(str), chunk[clab].astype(str)):
            truth[sid.split("#", 1)[0]] = tu.normalize_label(lab, level)
    return truth


def load_preds(path: str, truth_path: str, level: str,
               id_col: Optional[str], label_col: Optional[str],
               truth_id_col: str, truth_label_col: str
               ) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """One-hot scores from a per-model HARD-prediction TSV (id + label).

    A label-only model is a single point in ROC/PR space: its predicted class is
    scored 1.0 and every other class 0. Abstentions / missing predictions leave
    the row all-zero (a "no-call" that counts against recall for every class).
    """
    truth = _load_truth_map(truth_path, level, truth_id_col, truth_label_col)

    preds: Dict[str, str] = {}
    pid = plab = None
    for chunk in tu.iter_tsv_chunks(path):
        if pid is None:
            cols = list(chunk.columns)
            pid = id_col or tu.detect_column(
                cols, ["sample_id", "family_id", "id", "seq_id", "name", "query"])
            plab = label_col or tu.detect_column(
                cols, ["prediction", "pred", "label", "class", "classification",
                       "superfamily", "order", "final", "consensus"])
            if pid is None or plab is None:
                raise SystemExit(
                    f"ERROR: cannot locate id/label columns in {path}. "
                    f"Columns: {cols}")
        for sid, lab in zip(chunk[pid].astype(str), chunk[plab].astype(str)):
            preds[sid.split("#", 1)[0]] = tu.normalize_label(lab, level)

    return curves_from_preds(preds, truth)


def curves_from_preds(preds: Dict[str, str], truth: Dict[str, str]
                      ) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """One-hot score matrix from {id: label} hard predictions + truth.

    Shared by --preds (a 2-column file) and --votes --model <col> (one column of
    the wide voting table). Predicted class -> 1.0, others 0; abstain/out-of-set
    -> all-zero row.
    """
    ids = [s for s in preds if s in truth]
    if not ids:
        return [], np.asarray([]), np.zeros((0, 0))

    class_set = set()
    for s in ids:
        if preds[s].lower() not in _ABSTAIN:
            class_set.add(preds[s])
        class_set.add(truth[s])
    class_set = {c for c in class_set if c.lower() not in _ABSTAIN}
    preferred = tu.CANONICAL_CLASSES + tu.CANONICAL_ORDERS
    classes = [c for c in preferred if c in class_set] + \
              sorted(c for c in class_set if c not in preferred)

    y_true = np.array([truth[s] for s in ids])
    scores = np.zeros((len(ids), len(classes)), dtype=np.float64)
    cidx = {c: j for j, c in enumerate(classes)}
    n_abstain = sum(1 for s in ids if cidx.get(preds[s]) is None)
    for i, s in enumerate(ids):
        j = cidx.get(preds[s])
        if j is not None:                 # abstain / out-of-set -> all-zero row
            scores[i, j] = 1.0
    LOG.info("Hard predictions: %d ids matched truth, %d classes, %d abstentions",
             len(ids), len(classes), n_abstain)
    return classes, y_true, scores


def _read_header(path: str) -> List[str]:
    """Return the column names of a (optionally gzipped) TSV, cheaply."""
    with tu.open_text(path) as fh:
        first = fh.readline().rstrip("\n").rstrip("\r")
    return first.split("\t")


def read_voting_column(path: str, id_col: str, model_col: str, level: str
                       ) -> Dict[str, str]:
    """Extract one model column of a wide voting table as {id: norm_label}."""
    preds: Dict[str, str] = {}
    for chunk in tu.iter_tsv_chunks(path):
        for sid, lab in zip(chunk[id_col].astype(str), chunk[model_col].astype(str)):
            preds[sid.split("#", 1)[0]] = tu.normalize_label(lab, level)
    return preds


# --------------------------------------------------------------------------- #
# Curve maths
# --------------------------------------------------------------------------- #
def _thin(x: np.ndarray, max_points: int) -> np.ndarray:
    if max_points <= 0 or len(x) <= max_points:
        return np.arange(len(x))
    return np.unique(np.linspace(0, len(x) - 1, max_points).astype(int))


def compute_curves(classes: List[str], y_true: np.ndarray, scores: np.ndarray,
                   max_points: int) -> Tuple[List[dict], List[dict], List[dict]]:
    from sklearn.metrics import (roc_curve, auc, precision_recall_curve,
                                  average_precision_score)

    K = len(classes)
    # Build the one-vs-rest indicator directly (column j <-> classes[j]). We do
    # NOT use sklearn.label_binarize: for the binary case (K==2) it returns a
    # single column, which silently misaligns with the K-column `scores` matrix
    # and raises IndexError at --level class.
    cls_index = {c: j for j, c in enumerate(classes)}
    Y = np.zeros((len(y_true), K), dtype=int)
    for i, lab in enumerate(y_true):
        j = cls_index.get(lab)
        if j is not None:            # truth labels outside `classes` stay all-0
            Y[i, j] = 1
    roc_rows: List[dict] = []
    pr_rows: List[dict] = []
    auc_rows: List[dict] = []

    fpr_grid = np.linspace(0, 1, 200)
    rec_grid = np.linspace(0, 1, 200)
    tpr_interp: List[np.ndarray] = []
    prec_interp: List[np.ndarray] = []
    per_class_roc_auc: List[float] = []
    per_class_ap: List[float] = []

    for j, c in enumerate(classes):
        yj = Y[:, j]
        sj = scores[:, j]
        n_pos = int(yj.sum())
        if n_pos == 0 or n_pos == len(yj):
            auc_rows.append({"class": c, "n_pos": n_pos,
                             "roc_auc": float("nan"), "pr_auc": float("nan")})
            continue
        fpr, tpr, thr = roc_curve(yj, sj)
        roc_auc = auc(fpr, tpr)
        prec, rec, pthr = precision_recall_curve(yj, sj)
        ap = average_precision_score(yj, sj)
        per_class_roc_auc.append(roc_auc)
        per_class_ap.append(ap)

        # Interpolated for macro averaging.
        tpr_interp.append(np.interp(fpr_grid, fpr, tpr))
        # precision_recall_curve returns recall descending; reverse for interp.
        prec_interp.append(np.interp(rec_grid, rec[::-1], prec[::-1]))

        ridx = _thin(fpr, max_points)
        for i in ridx:
            roc_rows.append({"class": c, "fpr": fpr[i], "tpr": tpr[i],
                             "threshold": thr[i] if i < len(thr) else float("nan")})
        pidx = _thin(rec, max_points)
        pthr_full = np.append(pthr, float("nan"))
        for i in pidx:
            pr_rows.append({"class": c, "recall": rec[i], "precision": prec[i],
                            "threshold": pthr_full[i]})
        auc_rows.append({"class": c, "n_pos": n_pos,
                         "roc_auc": roc_auc, "pr_auc": ap})

    # ---- micro-average (pool all one-vs-rest decisions) ------------------- #
    y_flat = Y.ravel()
    s_flat = scores.ravel()
    if y_flat.sum() > 0:
        fpr, tpr, thr = roc_curve(y_flat, s_flat)
        micro_roc_auc = auc(fpr, tpr)
        prec, rec, pthr = precision_recall_curve(y_flat, s_flat)
        micro_ap = average_precision_score(y_flat, s_flat)
        for i in _thin(fpr, max_points):
            roc_rows.append({"class": "micro", "fpr": fpr[i], "tpr": tpr[i],
                             "threshold": thr[i] if i < len(thr) else float("nan")})
        pthr_full = np.append(pthr, float("nan"))
        for i in _thin(rec, max_points):
            pr_rows.append({"class": "micro", "recall": rec[i],
                            "precision": prec[i], "threshold": pthr_full[i]})
        auc_rows.append({"class": "micro", "n_pos": int(y_flat.sum()),
                         "roc_auc": micro_roc_auc, "pr_auc": micro_ap})

    # ---- macro-average (mean of interpolated per-class curves) ------------ #
    if tpr_interp:
        mean_tpr = np.mean(tpr_interp, axis=0)
        macro_roc_auc = auc(fpr_grid, mean_tpr)
        mean_prec = np.mean(prec_interp, axis=0)
        macro_ap = float(np.mean(per_class_ap))
        for i in _thin(fpr_grid, max_points):
            roc_rows.append({"class": "macro", "fpr": fpr_grid[i],
                             "tpr": mean_tpr[i], "threshold": float("nan")})
        for i in _thin(rec_grid, max_points):
            pr_rows.append({"class": "macro", "recall": rec_grid[i],
                            "precision": mean_prec[i], "threshold": float("nan")})
        auc_rows.append({"class": "macro", "n_pos": len(per_class_ap),
                         "roc_auc": macro_roc_auc, "pr_auc": macro_ap})
        auc_rows.append({"class": "macro_auc_mean", "n_pos": len(per_class_roc_auc),
                         "roc_auc": float(np.mean(per_class_roc_auc)),
                         "pr_auc": float(np.mean(per_class_ap))})

    return roc_rows, pr_rows, auc_rows


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
def _write(path: str, header: List[str], rows: List[dict]) -> None:
    import csv
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            out = []
            for h in header:
                v = r.get(h, "")
                out.append(f"{v:.6f}" if isinstance(v, float) else v)
            w.writerow(out)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute micro/macro ROC and PR (AP) curves + coordinates for "
                    "a TE classifier or the ensemble vote.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--proba", help="Wide probability TSV (prob:<class> columns).")
    src.add_argument("--preds", help="Per-model hard-prediction TSV, id + "
                                     "predicted label (needs --truth).")
    src.add_argument("--votes", help="Wide hard-label voting table (needs --truth).")
    p.add_argument("--truth", help="Ground-truth TSV (required for --preds/--votes).")
    p.add_argument("--model", default="model", help="Model name for output files.")
    p.add_argument("--tag", default="benchmark", help="Lineage/run tag for filenames.")
    p.add_argument("--level", default="order",
                   choices=["class", "order", "superfamily"])
    p.add_argument("--prob-prefix", default="prob:",
                   help="Column prefix identifying probability columns (proba mode).")
    p.add_argument("--id-col", default=None)
    p.add_argument("--true-col", default=None)
    p.add_argument("--pred-label-col", default=None,
                   help="Predicted-label column name (--preds mode; auto-detected).")
    p.add_argument("--truth-id-col", default="family_id")
    p.add_argument("--truth-label-col", default="label")
    p.add_argument("--vote-exclude-cols", nargs="*", default=["Ensemble"],
                   help="Voting-table columns to exclude from vote counting.")
    p.add_argument("--vote-denominator", default="fixed",
                   choices=["fixed", "cast"],
                   help="'fixed' divides by #model columns (abstentions lower the "
                        "score); 'cast' divides by votes actually cast.")
    p.add_argument("--max-points", type=int, default=500,
                   help="Max coordinates per curve (thinned for plotting; 0=all).")
    p.add_argument("--label-map", help="JSON of extra label synonyms.")
    p.add_argument("--outdir", default="curves_out")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)
    tu.ensure_dir(args.outdir)

    if args.votes:
        if not args.truth:
            raise SystemExit("ERROR: --votes mode requires --truth.")
        header = _read_header(args.votes)
        vid = args.id_col or tu.detect_column(
            header, ["repeat_id", "sample_id", "id", "family_id"]) or header[0]
        model_cols = [c for c in header if c != vid]
        # Strip stray whitespace/commas (e.g. a shell 'for x in a, b' loop).
        req = args.model.strip().strip(",").strip()
        match = next((c for c in model_cols if c.lower() == req.lower()), None)

        if req.lower() == "ensemble":
            # Continuous ensemble vote-fraction curve (uses ALL model columns).
            classes, y_true, scores = load_votes(
                args.votes, args.truth, args.level, args.id_col,
                args.truth_id_col, args.truth_label_col,
                args.vote_exclude_cols, args.vote_denominator)
            score_type = "vote-fraction"
        elif match is not None:
            # One specific model's hard-prediction curve (single operating point).
            args.model = match          # clean name for output files
            LOG.info("Per-model curve from voting-table column '%s'", match)
            truth = _load_truth_map(args.truth, args.level,
                                    args.truth_id_col, args.truth_label_col)
            preds = read_voting_column(args.votes, vid, match, args.level)
            classes, y_true, scores = curves_from_preds(preds, truth)
            score_type = "hard"
        else:
            raise SystemExit(
                f"ERROR: --model '{args.model}' is not a column in {args.votes}. "
                f"Use --model Ensemble for the vote-fraction curve, or one of: "
                f"{', '.join(model_cols)}.")
    elif args.preds:
        if not args.truth:
            raise SystemExit("ERROR: --preds mode requires --truth.")
        classes, y_true, scores = load_preds(
            args.preds, args.truth, args.level, args.id_col,
            args.pred_label_col, args.truth_id_col, args.truth_label_col)
        score_type = "hard"
    else:
        classes, y_true, scores = load_proba(
            args.proba, args.level, args.id_col, args.true_col, args.prob_prefix)
        score_type = "proba"

    if scores.shape[0] == 0:
        raise SystemExit("ERROR: no samples matched between scores and truth.")
    LOG.info("Curves for model '%s': %d samples x %d classes",
             args.model, scores.shape[0], len(classes))

    roc_rows, pr_rows, auc_rows = compute_curves(
        classes, y_true, scores, args.max_points)

    safe = args.model.replace("/", "_").replace(" ", "_")
    _write(os.path.join(args.outdir, f"{args.tag}.{safe}.roc_coords.tsv"),
           ["class", "fpr", "tpr", "threshold"], roc_rows)
    _write(os.path.join(args.outdir, f"{args.tag}.{safe}.pr_coords.tsv"),
           ["class", "recall", "precision", "threshold"], pr_rows)
    # Stamp the model name + how it was scored into the AUC summary. 'score_type'
    # (hard | vote-fraction | proba) lets plot_benchmarks tag the legend so a soft
    # (vote-fraction) curve is never compared to hard single-point curves unlabelled.
    for r in auc_rows:
        r["model"] = args.model
        r["score_type"] = score_type
    _write(os.path.join(args.outdir, f"{args.tag}.{safe}.curve_auc.tsv"),
           ["model", "score_type", "class", "n_pos", "roc_auc", "pr_auc"], auc_rows)
    LOG.info("Scored '%s' as %s", args.model, score_type)

    for r in auc_rows:
        if r["class"] in ("micro", "macro"):
            LOG.info("  %-6s ROC-AUC=%.4f  PR-AUC(AP)=%.4f",
                     r["class"], r["roc_auc"], r["pr_auc"])
    LOG.info("Coordinates + AUC summary written to %s/", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
