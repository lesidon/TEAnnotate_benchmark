#!/usr/bin/env python3
"""
evaluate_classification.py — multi-class evaluation of TE classifiers and the
ensemble vote against a Dfam ground-truth partition.

Given a ground-truth TSV (from dfam_partitioner.py) and one or more prediction
sources, it computes — per model AND for the ensemble — a confusion matrix and
the full metric panel:

  per-class : support, TP, FP, FN, precision, recall, F1
  aggregate : accuracy, balanced accuracy, macro/micro/weighted precision-recall-F1,
              Matthews correlation coefficient (MCC), Cohen's kappa, coverage

All aggregate metrics are derived directly from the confusion matrix, which is
accumulated incrementally, so memory is O(#classes^2) regardless of dataset size.

Prediction sources (any combination):
  --voting-table FILE            wide table: id + one column per model (+ Ensemble)
  --predictions NAME=FILE ...     per-model file (.fasta/.tsv/.gff3, auto-detected)

Outputs (per lineage tag) to --outdir:
  <tag>.<model>.confusion_matrix.tsv
  <tag>.<model>.per_class_metrics.tsv
  <tag>.<model>.summary_metrics.tsv
  <tag>.all_models.summary_metrics.tsv     (tidy, one row per model — for plotting)

Example
-------
  python evaluate_classification.py \
      --truth dfam_partitions/dfam_fungi.tsv \
      --voting-table results/sampleX.voting_table.tsv \
      --predictions teclass2=preds/tc2.tsv deepte=preds/deepte.fasta \
      --level superfamily --lineage fungi --outdir eval_out
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import te_label_utils as tu

LOG = tu.setup_logger("evaluate_classification")

_ABSTAIN = {"unknown", "other", "na", "n/a", "none", "", "."}


# --------------------------------------------------------------------------- #
# Prediction loaders  ->  {sample_id: raw_label}
# --------------------------------------------------------------------------- #
def load_predictions_tsv(path: str, id_col: Optional[str],
                         label_col: Optional[str]) -> Dict[str, str]:
    preds: Dict[str, str] = {}
    first = True
    cid = idx_label = None
    for chunk in tu.iter_tsv_chunks(path):
        if first:
            cols = list(chunk.columns)
            cid = id_col or tu.detect_column(
                cols, ["sample_id", "family_id", "id", "seq_id", "name", "query"])
            idx_label = label_col or tu.detect_column(
                cols, ["prediction", "pred", "label", "class", "classification",
                       "superfamily", "final", "consensus"])
            if cid is None or idx_label is None:
                raise SystemExit(
                    f"ERROR: cannot locate id/label columns in {path}. "
                    f"Columns: {cols}")
            first = False
        for sid, lab in zip(chunk[cid].astype(str), chunk[idx_label].astype(str)):
            preds[sid.split("#", 1)[0]] = lab
    return preds


def load_predictions_fasta(path: str) -> Dict[str, str]:
    from Bio import SeqIO
    preds: Dict[str, str] = {}
    with tu.open_text(path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            header = rec.description if rec.description else rec.id
            sid = rec.id.split("#", 1)[0]
            # Label = the classified header (parse_raw_label handles '#Class/Sub').
            preds[sid] = header
    return preds


def load_predictions_gff3(path: str, id_attr: str, label_attr: str) -> Dict[str, str]:
    preds: Dict[str, str] = {}
    with tu.open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = {}
            for kv in f[8].split(";"):
                kv = kv.strip()
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip().lower()] = v.strip()
            sid = attrs.get(id_attr.lower()) or attrs.get("id") or attrs.get("name")
            lab = (attrs.get(label_attr.lower()) or attrs.get("classification")
                   or attrs.get("class") or attrs.get("target") or attrs.get("name"))
            if sid and lab:
                preds[sid.split("#", 1)[0]] = lab
    return preds


def load_prediction_source(path: str, fmt: str,
                           gff_id_attr: str, gff_label_attr: str,
                           id_col: Optional[str],
                           label_col: Optional[str]) -> Dict[str, str]:
    if fmt == "auto":
        low = path.lower().rstrip(".gz")
        if low.endswith((".fa", ".fasta", ".fna")):
            fmt = "fasta"
        elif low.endswith((".gff", ".gff3")):
            fmt = "gff3"
        else:
            fmt = "tsv"
    if fmt == "fasta":
        return load_predictions_fasta(path)
    if fmt == "gff3":
        return load_predictions_gff3(path, gff_id_attr, gff_label_attr)
    return load_predictions_tsv(path, id_col, label_col)


def load_voting_table(path: str) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """Return (model_names, {model: {id: raw_label}}) from a wide voting table."""
    models: List[str] = []
    per_model: Dict[str, Dict[str, str]] = {}
    id_col = None
    for chunk in tu.iter_tsv_chunks(path):
        if id_col is None:
            cols = list(chunk.columns)
            id_col = tu.detect_column(cols, ["repeat_id", "sample_id", "id", "family_id"])
            if id_col is None:
                id_col = cols[0]
            models = [c for c in cols if c != id_col]
            per_model = {m: {} for m in models}
        for _, row in chunk.iterrows():
            sid = str(row[id_col]).split("#", 1)[0]
            for m in models:
                per_model[m][sid] = str(row[m])
    return models, per_model


# --------------------------------------------------------------------------- #
# Metrics from a confusion matrix
# --------------------------------------------------------------------------- #
def _safe_div(a: float, b: float) -> float:
    return a / b if b else 0.0


def metrics_from_confusion(conf: Dict[Tuple[str, str], int],
                           classes: List[str]) -> Tuple[dict, List[dict]]:
    idx = {c: i for i, c in enumerate(classes)}
    k = len(classes)
    C = [[0] * k for _ in range(k)]
    for (t, p), n in conf.items():
        if t in idx and p in idx:
            C[idx[t]][idx[p]] += n

    row_t = [sum(C[i]) for i in range(k)]          # true totals (support)
    col_p = [sum(C[i][j] for i in range(k)) for j in range(k)]  # predicted totals
    diag = [C[i][i] for i in range(k)]
    total = sum(row_t)

    per_class: List[dict] = []
    macro_p = macro_r = macro_f = 0.0
    weighted_f = 0.0
    bal_acc_num = 0
    supported = 0
    for i, c in enumerate(classes):
        tp = diag[i]
        fp = col_p[i] - tp
        fn = row_t[i] - tp
        prec = _safe_div(tp, tp + fp)
        rec = _safe_div(tp, tp + fn)
        f1 = _safe_div(2 * prec * rec, prec + rec)
        per_class.append({
            "class": c, "support": row_t[i], "tp": tp, "fp": fp, "fn": fn,
            "precision": prec, "recall": rec, "f1": f1,
        })
        macro_p += prec
        macro_r += rec
        macro_f += f1
        weighted_f += _safe_div(row_t[i], total) * f1
        if row_t[i] > 0:
            bal_acc_num += rec
            supported += 1

    macro_p /= k or 1
    macro_r /= k or 1
    macro_f /= k or 1
    accuracy = _safe_div(sum(diag), total)          # == micro-F1
    balanced_accuracy = _safe_div(bal_acc_num, supported)

    # Multi-class MCC (Gorodkin) from the confusion matrix.
    s = total
    c_correct = sum(diag)
    sum_pt = sum(col_p[i] * row_t[i] for i in range(k))
    sum_p2 = sum(v * v for v in col_p)
    sum_t2 = sum(v * v for v in row_t)
    mcc_den = math.sqrt((s * s - sum_p2) * (s * s - sum_t2))
    mcc = _safe_div(c_correct * s - sum_pt, mcc_den) if mcc_den else 0.0

    # Cohen's kappa vs. the (implicit) ground truth marginals.
    p_o = accuracy
    p_e = _safe_div(sum(col_p[i] * row_t[i] for i in range(k)), s * s)
    kappa = _safe_div(p_o - p_e, 1 - p_e)

    summary = {
        "n": total,
        "accuracy": accuracy,
        "balanced_accuracy": balanced_accuracy,
        "macro_precision": macro_p,
        "macro_recall": macro_r,
        "macro_f1": macro_f,
        "micro_f1": accuracy,
        "weighted_f1": weighted_f,
        "mcc": mcc,
        "cohen_kappa": kappa,
    }
    return summary, per_class


def order_classes(classes) -> List[str]:
    preferred = tu.CANONICAL_CLASSES + tu.CANONICAL_ORDERS
    canon = [c for c in preferred if c in classes]
    rest = sorted(c for c in classes if c not in preferred)
    return canon + rest


# --------------------------------------------------------------------------- #
# Evaluate one model
# --------------------------------------------------------------------------- #
def evaluate_model(truth: Dict[str, str], preds_raw: Dict[str, str],
                   level: str, abstain_policy: str) -> Tuple[dict, List[dict],
                                                             Dict[Tuple[str, str], int],
                                                             List[str], dict]:
    conf: Dict[Tuple[str, str], int] = Counter()
    class_set = set()
    n_total = n_attempted = n_missing = 0
    for sid, t_raw in truth.items():
        t = tu.normalize_label(t_raw, level)
        n_total += 1
        raw_pred = preds_raw.get(sid)
        if raw_pred is None:
            n_missing += 1
            p = "Unknown"
            is_abstain = True
        else:
            p = tu.normalize_label(raw_pred, level)
            is_abstain = p.lower() in _ABSTAIN
        if not is_abstain:
            n_attempted += 1
        if is_abstain and abstain_policy == "drop":
            continue
        if is_abstain:
            p = "Unknown"
        conf[(t, p)] += 1
        class_set.add(t)
        class_set.add(p)

    classes = order_classes(class_set)
    summary, per_class = metrics_from_confusion(conf, classes)
    cov = _safe_div(n_attempted, n_total)
    summary["coverage"] = cov
    summary["n_total_truth"] = n_total
    summary["n_attempted"] = n_attempted
    summary["n_missing"] = n_missing
    counts = {"n_total": n_total, "n_attempted": n_attempted, "n_missing": n_missing}
    return summary, per_class, conf, classes, counts


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
def write_confusion(path: str, conf: Dict[Tuple[str, str], int],
                    classes: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["true\\pred"] + classes)
        for t in classes:
            w.writerow([t] + [conf.get((t, p), 0) for p in classes])


def write_per_class(path: str, per_class: List[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["class", "support", "tp", "fp", "fn",
                    "precision", "recall", "f1"])
        for r in per_class:
            w.writerow([r["class"], r["support"], r["tp"], r["fp"], r["fn"],
                        f"{r['precision']:.6f}", f"{r['recall']:.6f}",
                        f"{r['f1']:.6f}"])


SUMMARY_FIELDS = [
    "n", "n_total_truth", "n_attempted", "n_missing", "coverage",
    "accuracy", "balanced_accuracy", "macro_precision", "macro_recall",
    "macro_f1", "micro_f1", "weighted_f1", "mcc", "cohen_kappa",
]


def write_summary(path: str, summary: dict) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for key in SUMMARY_FIELDS:
            val = summary.get(key, "")
            w.writerow([key, f"{val:.6f}" if isinstance(val, float) else val])


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Evaluate individual TE classifiers and the ensemble vote "
                    "against a Dfam ground-truth partition.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--truth", required=True,
                   help="Ground-truth TSV from dfam_partitioner.py.")
    p.add_argument("--truth-id-col", default="family_id")
    p.add_argument("--truth-label-col", default="label",
                   help="Raw label column in the truth TSV (normalised at eval).")
    p.add_argument("--voting-table",
                   help="Wide ensemble voting table (id + per-model columns).")
    p.add_argument("--predictions", nargs="+", default=[],
                   metavar="NAME=FILE",
                   help="Per-model prediction files, e.g. teclass2=preds.tsv "
                        "deepte=preds.fasta.")
    p.add_argument("--pred-format", default="auto",
                   choices=["auto", "tsv", "fasta", "gff3"])
    p.add_argument("--pred-id-col", default=None)
    p.add_argument("--pred-label-col", default=None)
    p.add_argument("--gff-id-attr", default="ID")
    p.add_argument("--gff-label-attr", default="classification")
    p.add_argument("--level", default="order",
                   choices=["class", "order", "superfamily"],
                   help="TE hierarchy level to score: class (ClassI/ClassII), "
                        "order (LTR/LINE/DNA/...) or superfamily.")
    p.add_argument("--abstain-policy", default="count",
                   choices=["count", "drop"],
                   help="'count' penalises abstentions (pred=Unknown); 'drop' "
                        "scores only attempted calls (report coverage separately).")
    p.add_argument("--keep-nonte-truth", action="store_true",
                   help="Keep truth families that normalise to Unknown/Other "
                        "(dropped by default).")
    p.add_argument("--keep-list", default=None,
                   help="File of family ids (or a FASTA) to RESTRICT evaluation to "
                        "— e.g. leakage_filter.py's *.heldout.ids for a leakage-free "
                        "score without re-running the classifiers.")
    p.add_argument("--exclude-list", default=None,
                   help="File of family ids (or a FASTA) to DROP from evaluation "
                        "(e.g. leaked/train-overlapping families).")
    p.add_argument("--lineage", default=None,
                   help="Tag for output filenames (default: truth file stem).")
    p.add_argument("--label-map", help="JSON of extra label synonyms.")
    p.add_argument("--outdir", default="eval_out")
    p.add_argument("--verbose", action="store_true")
    return p


def load_id_list(path: Optional[str]) -> Optional[set]:
    """Read a set of family ids from a plain id list or a FASTA (>id ...)."""
    if not path:
        return None
    ids: set = set()
    with tu.open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            token = line[1:] if line.startswith(">") else line
            token = token.split()[0].split("#", 1)[0]
            if token:
                ids.add(token)
    return ids


def load_truth(path: str, id_col: str, label_col: str, level: str,
               keep_nonte: bool, keep_ids: Optional[set] = None,
               exclude_ids: Optional[set] = None) -> Dict[str, str]:
    truth: Dict[str, str] = {}
    first = True
    cid = clab = None
    for chunk in tu.iter_tsv_chunks(path):
        if first:
            cols = list(chunk.columns)
            cid = tu.detect_column(cols, [id_col, "family_id", "id", "sample_id"])
            clab = tu.detect_column(cols, [label_col, "label", "label_order",
                                           "rm_type", "classification", "class"])
            if cid is None or clab is None:
                raise SystemExit(f"ERROR: truth id/label columns not found in {cols}")
            first = False
        for sid, lab in zip(chunk[cid].astype(str), chunk[clab].astype(str)):
            base = sid.split("#", 1)[0]
            # Leakage / subset control: restrict to keep_ids, drop exclude_ids.
            if keep_ids is not None and base not in keep_ids:
                continue
            if exclude_ids is not None and base in exclude_ids:
                continue
            norm = tu.normalize_label(lab, level)
            if not keep_nonte and norm in ("Unknown", "Other"):
                continue
            truth[base] = lab
    return truth


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)
    tu.ensure_dir(args.outdir)

    tag = args.lineage or os.path.splitext(os.path.basename(args.truth))[0]

    keep_ids = load_id_list(args.keep_list)
    exclude_ids = load_id_list(args.exclude_list)
    if keep_ids is not None:
        LOG.info("Restricting evaluation to %d ids from %s",
                 len(keep_ids), args.keep_list)
    if exclude_ids is not None:
        LOG.info("Excluding %d ids from %s", len(exclude_ids), args.exclude_list)

    truth = load_truth(args.truth, args.truth_id_col, args.truth_label_col,
                       args.level, args.keep_nonte_truth, keep_ids, exclude_ids)
    if not truth:
        raise SystemExit("ERROR: no usable ground-truth rows loaded.")
    LOG.info("Loaded %d ground-truth families (level=%s, tag=%s)",
             len(truth), args.level, tag)

    # Assemble {model_name: {id: raw_pred}} from all sources.
    model_preds: Dict[str, Dict[str, str]] = {}
    if args.voting_table:
        models, per_model = load_voting_table(args.voting_table)
        for m in models:
            model_preds[m] = per_model[m]
        LOG.info("Voting table contributed models: %s", ", ".join(models))
    for item in args.predictions:
        if "=" not in item:
            raise SystemExit(f"ERROR: --predictions entry '{item}' is not NAME=FILE")
        name, path = item.split("=", 1)
        model_preds[name] = load_prediction_source(
            path, args.pred_format, args.gff_id_attr, args.gff_label_attr,
            args.pred_id_col, args.pred_label_col)
        LOG.info("Loaded %d predictions for model '%s'", len(model_preds[name]), name)

    if not model_preds:
        raise SystemExit("ERROR: no prediction sources given "
                         "(use --voting-table and/or --predictions).")

    combined_rows: List[dict] = []
    for model, preds in model_preds.items():
        summary, per_class, conf, classes, counts = evaluate_model(
            truth, preds, args.level, args.abstain_policy)
        safe = model.replace("/", "_").replace(" ", "_")
        write_confusion(os.path.join(args.outdir, f"{tag}.{safe}.confusion_matrix.tsv"),
                        conf, classes)
        write_per_class(os.path.join(args.outdir, f"{tag}.{safe}.per_class_metrics.tsv"),
                        per_class)
        write_summary(os.path.join(args.outdir, f"{tag}.{safe}.summary_metrics.tsv"),
                      summary)
        row = {"model": model, "lineage": tag, "level": args.level}
        row.update({k: summary.get(k, "") for k in SUMMARY_FIELDS})
        combined_rows.append(row)
        LOG.info("%-12s | n=%d acc=%.3f macroF1=%.3f MCC=%.3f cov=%.3f",
                 model, summary["n"], summary["accuracy"], summary["macro_f1"],
                 summary["mcc"], summary["coverage"])

    # Tidy combined table (one row per model) for downstream plotting/stats.
    combined = os.path.join(args.outdir, f"{tag}.all_models.summary_metrics.tsv")
    with open(combined, "w", encoding="utf-8", newline="") as fh:
        header = ["model", "lineage", "level"] + SUMMARY_FIELDS
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in combined_rows:
            w.writerow([r.get(h, "") if not isinstance(r.get(h), float)
                        else f"{r[h]:.6f}" for h in header])
    LOG.info("Wrote combined summary -> %s", combined)
    return 0


if __name__ == "__main__":
    sys.exit(main())
