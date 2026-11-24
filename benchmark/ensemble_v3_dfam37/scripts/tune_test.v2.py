#!/usr/bin/env python3
"""tune_test.py - per-kingdom ensemble tuning on TRAIN and evaluation on TEST.

Every candidate ensemble is decided by the pipeline's own voter (v9 bin/ensemble_vote.py:
call_with_agreement, with the abstention gate computed as resolve_gates does), so a
tuned setting transfers to the pipeline unchanged. Labels of truth and predictions are
resolved by ens_labels.order_sf (one convention for every source).

Search space (TRAIN only):
  weighted voting: every non-empty model subset x weights {1, 1.5, 2} per active model
                   x min_frac {0.3, 0.4, 0.5, 0.6} x min_agree {1, 2}
  hard voting:     every subset x min_frac x min_agree
Objective: order-level macro-F1 over truth orders with >= --min-support sequences;
an abstention (no call) counts as a miss. Ties: higher MCC, then fewer models.

TEST: single models, the current v9 configuration, all five models with v9 weights, the
tuned configuration, and a stacking reference (gradient boosting on one-hot model calls,
fit on TRAIN; it always answers). 95% CIs from 1000 cluster bootstraps; McNemar tests of
the tuned ensemble against the best single model and against v9.

usage: tune_test.py --kingdom fungi --workdir ensemble_v2 --voter pipeline_head_v9/bin
"""
import argparse
import csv
import itertools
import json
import math
import os
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ens_labels import load_label_overrides, order_sf  # noqa: E402

MODELS = ["teclass2", "deepte", "create", "neuralte", "terrier"]
DISPLAY = {"teclass2": "TEclass2", "deepte": "DeepTE", "create": "CREATE",
           "neuralte": "NeuralTE", "terrier": "Terrier"}
V9_WEIGHTS = {"neuralte": 1.5, "terrier": 1.5, "create": 1.3, "teclass2": 1.2, "deepte": 1.1}
V9_SUBSET = {"animal": ["deepte", "terrier"], "plant": ["create", "deepte", "teclass2", "terrier"],
             "fungi": ["terrier"]}
V9_MIN_FRAC, V9_MIN_AGREE = 0.5, 2
EV = None   # pipeline voter module, imported in main / workers


def load_voter(path):
    global EV
    sys.path.insert(0, path)
    import ensemble_vote as ev  # noqa: E402
    EV = ev


def decide(labels, cfg):
    """One ensemble call with the pipeline voter; labels = tuple aligned with MODELS."""
    include = set(cfg["subset"])
    ml = {m: l for m, l in zip(MODELS, labels) if l}
    weights = {m: cfg["weights"].get(m, 1.0) for m in MODELS}
    n = len(include)
    req_votes = max(1, int(math.ceil(cfg["min_frac"] * n)))
    req_w = cfg["min_frac"] * sum(weights[m] for m in include)
    return EV.call_with_agreement(cfg["method"], ml, req_votes, weights, req_w, include,
                                  cfg["min_agree"])


def to_order(label):
    return order_sf(label)[0] if label else "Unknown"


def metrics(truth, pred, classes):
    """macro-F1 over `classes`, multiclass MCC, accuracy, coverage. Abstain = 'Unknown'."""
    labs = sorted(set(truth) | set(pred))
    idx = {l: i for i, l in enumerate(labs)}
    cm = np.zeros((len(labs), len(labs)))
    for t, p in zip(truth, pred):
        cm[idx[t], idx[p]] += 1
    return metrics_from_cm(cm, labs, classes)


def metrics_from_cm(cm, labs, classes):
    tot = cm.sum()
    tp = np.diag(cm)
    f1s = []
    for c in classes:
        if c not in labs:
            continue
        i = labs.index(c)
        p = tp[i] / cm[:, i].sum() if cm[:, i].sum() else 0.0
        r = tp[i] / cm[i, :].sum() if cm[i, :].sum() else 0.0
        f1s.append(2 * p * r / (p + r) if p + r else 0.0)
    t, pk = cm.sum(1), cm.sum(0)
    c_ = tp.sum()
    num = c_ * tot - (t * pk).sum()
    den = math.sqrt(max(tot ** 2 - (pk ** 2).sum(), 0) * max(tot ** 2 - (t ** 2).sum(), 0))
    unk = labs.index("Unknown") if "Unknown" in labs else None
    called = tot - (cm[:, unk].sum() if unk is not None else 0)
    return {"macro_f1": float(np.mean(f1s)) if f1s else 0.0, "mcc": float(num / den) if den else 0.0,
            "accuracy": float(c_ / tot) if tot else 0.0, "coverage": float(called / tot) if tot else 0.0,
            "precision_called": float(c_ / called) if called else 0.0}


def eval_cfg(args):
    cfg, patterns, classes = args
    cm_rows = defaultdict(Counter)
    for labels, tcount in patterns:
        p = to_order(decide(labels, cfg))
        for t, n in tcount.items():
            cm_rows[t][p] += n
    labs = sorted(set(cm_rows) | {p for r in cm_rows.values() for p in r})
    cm = np.array([[cm_rows[t][p] for p in labs] for t in labs], dtype=float)
    return cfg, metrics_from_cm(cm, labs, classes)


def search_space():
    cfgs = []
    for k in range(1, len(MODELS) + 1):
        for sub in itertools.combinations(MODELS, k):
            for mf in (0.3, 0.4, 0.5, 0.6):
                for ma in (1, 2):
                    if ma > len(sub):
                        continue
                    cfgs.append({"method": "hard", "subset": list(sub), "weights": {},
                                 "min_frac": mf, "min_agree": ma})
                    grids = [(1.0,)] if len(sub) == 1 else [(1.0, 1.5, 2.0)] * len(sub)
                    for ws in itertools.product(*grids):
                        if len(sub) > 1 and min(ws) != 1.0:      # weights are scale-free: anchor min at 1
                            continue
                        cfgs.append({"method": "weighted", "subset": list(sub),
                                     "weights": dict(zip(sub, ws)), "min_frac": mf, "min_agree": ma})
    return cfgs


def cfg_name(cfg):
    w = ",".join("%s=%g" % (DISPLAY[m], cfg["weights"].get(m, 1.0)) for m in cfg["subset"])
    return "%s[%s] frac=%.1f agree=%d" % (cfg["method"], w if cfg["method"] == "weighted"
                                           else ",".join(DISPLAY[m] for m in cfg["subset"]),
                                           cfg["min_frac"], cfg["min_agree"])


def bootstrap_ci(truth, preds, clusters, classes, n_boot=1000, seed=0):
    """95% percentile CIs of macro-F1 / MCC per method, resampling clusters."""
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(clusters, return_inverse=True)
    C = len(uniq)
    labs = sorted(set(truth) | {p for v in preds.values() for p in v})
    li = {l: i for i, l in enumerate(labs)}
    t_idx = np.array([li[t] for t in truth])
    p_idx = {m: np.array([li[p] for p in v]) for m, v in preds.items()}
    out = {m: {"macro_f1": [], "mcc": []} for m in preds}
    L = len(labs)
    for _ in range(n_boot):
        # drawing clusters with replacement == weighting each sequence by its cluster's draw count
        w = np.bincount(rng.integers(0, C, C), minlength=C)[inv].astype(float)
        for m, pv in p_idx.items():
            cm = np.bincount(t_idx * L + pv, weights=w, minlength=L * L).reshape(L, L)
            r = metrics_from_cm(cm, labs, classes)
            out[m]["macro_f1"].append(r["macro_f1"]); out[m]["mcc"].append(r["mcc"])
    return {m: {k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) for k, v in d.items()}
            for m, d in out.items()}


def mcnemar(truth, a, b):
    ca = np.array([x == y for x, y in zip(truth, a)]); cb = np.array([x == y for x, y in zip(truth, b)])
    n01 = int((ca & ~cb).sum()); n10 = int((~ca & cb).sum())
    if n01 + n10 == 0:
        return n01, n10, 1.0
    from scipy.stats import chi2
    stat = (abs(n01 - n10) - 1) ** 2 / (n01 + n10)
    return n01, n10, float(chi2.sf(stat, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kingdom", required=True, choices=["plant", "animal", "fungi"])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--voter", required=True, help="pipeline bin/ dir holding ensemble_vote.py")
    ap.add_argument("--label-map", default=None)
    ap.add_argument("--min-support", type=int, default=20)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--n-boot", type=int, default=1000)
    a = ap.parse_args()
    load_label_overrides(a.label_map)
    load_voter(a.voter)
    k, W = a.kingdom, a.workdir

    truth = {r["family_id"]: r for r in csv.DictReader(open("%s/pools/%s.truth.tsv" % (W, k)), delimiter="\t")}
    split = {r["family_id"]: r for r in csv.DictReader(open("%s/split/%s.split.tsv" % (W, k)), delimiter="\t")}
    vt = list(csv.DictReader(open("%s/preds/%s/voting_table.tsv" % (W, k)), delimiter="\t"))
    rows = []
    for r in vt:
        fid = r["repeat_id"]
        if fid not in truth or fid not in split:
            continue
        labels = tuple("" if r.get(DISPLAY[m], "NA") in ("NA", "") else r[DISPLAY[m]] for m in MODELS)
        rows.append((fid, labels, truth[fid]["order"], truth[fid]["superfamily"],
                     split[fid]["split"], split[fid]["cluster"]))
    tr = [x for x in rows if x[4] == "train"]
    te = [x for x in rows if x[4] == "test"]
    print("%s: %d with predictions (train %d, test %d)" % (k, len(rows), len(tr), len(te)), flush=True)

    # ---- TRAIN: search -------------------------------------------------------
    sup = Counter(x[2] for x in tr)
    classes = sorted(c for c, n in sup.items() if n >= a.min_support)
    pat = defaultdict(Counter)
    for _, labels, o, _, _, _ in tr:
        pat[labels][o] += 1
    patterns = list(pat.items())
    cfgs = search_space()
    print("train: %d unique prediction patterns, %d candidate ensembles, classes %s"
          % (len(patterns), len(cfgs), classes), flush=True)
    with Pool(a.threads, initializer=load_voter, initargs=(a.voter,)) as pool:
        res = pool.map(eval_cfg, [(c, patterns, classes) for c in cfgs], chunksize=8)
    res.sort(key=lambda x: (-round(x[1]["macro_f1"], 4), -round(x[1]["mcc"], 4), len(x[0]["subset"])))
    os.makedirs("%s/tune" % W, exist_ok=True)
    with open("%s/tune/%s.top_configs.tsv" % (W, k), "w") as fh:
        fh.write("rank\tconfig\tn_models\tmacro_f1\tmcc\taccuracy\tcoverage\tprecision_called\n")
        for i, (c, m) in enumerate(res[:100], 1):
            fh.write("%d\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (
                i, cfg_name(c), len(c["subset"]), m["macro_f1"], m["mcc"], m["accuracy"],
                m["coverage"], m["precision_called"]))
    best = res[0][0]
    json.dump(best, open("%s/tune/%s.best.json" % (W, k), "w"), indent=1)
    print("best on train: %s  %s" % (cfg_name(best), res[0][1]), flush=True)

    # ---- TEST ------------------------------------------------------------------
    t_truth = [x[2] for x in te]
    t_sf = [x[3] for x in te]
    clusters = np.array([x[5] for x in te])
    sup_te = Counter(t_truth)
    classes_te = sorted(c for c, n in sup_te.items() if n >= a.min_support)
    v9 = {"method": "weighted", "subset": V9_SUBSET[k], "weights": V9_WEIGHTS,
          "min_frac": V9_MIN_FRAC, "min_agree": V9_MIN_AGREE}
    all5 = {"method": "weighted", "subset": MODELS, "weights": V9_WEIGHTS,
            "min_frac": V9_MIN_FRAC, "min_agree": V9_MIN_AGREE}
    raw = {}
    for m in MODELS:
        raw[DISPLAY[m]] = [x[1][MODELS.index(m)] or "" for x in te]
    for name, cfg in (("v9 (current)", v9), ("all five, v9 weights", all5), ("tuned on train", best)):
        raw[name] = [decide(x[1], cfg) or "" for x in te]
    # stacking reference (fit on TRAIN)
    try:
        from sklearn.ensemble import HistGradientBoostingClassifier
        vocab = sorted({to_order(l) for x in rows for l in x[1]})
        def feats(x):
            v = np.zeros(len(MODELS) * len(vocab))
            for j, l in enumerate(x[1]):
                v[j * len(vocab) + vocab.index(to_order(l))] = 1
            return v
        Xtr = np.array([feats(x) for x in tr]); ytr = [x[2] for x in tr]
        clf = HistGradientBoostingClassifier(max_iter=300, random_state=0).fit(Xtr, ytr)
        raw["stacking (reference)"] = list(clf.predict(np.array([feats(x) for x in te])))
    except Exception as e:  # pragma: no cover
        print("stacking skipped:", e)
    preds = {m: [to_order(l) for l in v] for m, v in raw.items()}
    sf_preds = {m: [order_sf(l)[1] if l else "Unknown" for l in v] for m, v in raw.items()}
    # superfamilies no classifier ever emits in this kingdom (e.g. LARD) are unreachable
    # for every method, so they are left out of the superfamily score and reported
    emitted = {order_sf(l)[1] for x in rows for l in x[1] if l}
    sf_cand = sorted(c for c, n in Counter(t_sf).items() if n >= a.min_support and c not in ("LTR", "DNA", "LINE", "SINE", "RC"))
    sf_classes = [c for c in sf_cand if c in emitted]
    sf_dropped = ["%s(n=%d)" % (c, Counter(t_sf)[c]) for c in sf_cand if c not in emitted]
    ci = bootstrap_ci(t_truth, preds, clusters, classes_te, n_boot=a.n_boot)
    singles = [DISPLAY[m] for m in MODELS]
    best_single = max(singles, key=lambda m: metrics(t_truth, preds[m], classes_te)["macro_f1"])
    os.makedirs("%s/test" % W, exist_ok=True)
    with open("%s/test/%s.results.tsv" % (W, k), "w") as fh:
        fh.write("method\tmacro_f1\tmacro_f1_ci\tmcc\tmcc_ci\taccuracy\tcoverage\tprecision_called\tsf_macro_f1\n")
        for m in preds:
            r = metrics(t_truth, preds[m], classes_te)
            # stacking predicts orders only, so it has no superfamily score
            s = (metrics(t_sf, sf_preds[m], sf_classes) if sf_classes and not m.startswith("stacking")
                 else {"macro_f1": float("nan")})
            fh.write("%s\t%.4f\t%.3f-%.3f\t%.4f\t%.3f-%.3f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (
                m, r["macro_f1"], *ci[m]["macro_f1"], r["mcc"], *ci[m]["mcc"], r["accuracy"],
                r["coverage"], r["precision_called"], s["macro_f1"]))
        for ref in (best_single, "v9 (current)"):
            n01, n10, p = mcnemar(t_truth, preds["tuned on train"], preds[ref])
            fh.write("# McNemar tuned vs %s: tuned-only correct %d, %s-only correct %d, p=%.2e\n" % (ref, n01, ref, n10, p))
        fh.write("# tuned config: %s\n# test classes (support >= %d): %s\n# superfamily classes: %s\n" % (
            cfg_name(best), a.min_support, ",".join(classes_te), ",".join(sf_classes)))
        fh.write("# superfamilies dropped (never emitted by any classifier): %s\n" % (",".join(sf_dropped) or "none"))
    print(open("%s/test/%s.results.tsv" % (W, k)).read(), flush=True)


if __name__ == "__main__":
    main()

