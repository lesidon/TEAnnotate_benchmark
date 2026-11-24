#!/usr/bin/env python3
"""two_step.py - two-step ensemble vote: Order first, then a Superfamily INSIDE that Order.

Step 1 (Order): the pipeline voter (call_with_agreement) on order-only votes. Every
classifier label is first resolved to (order, superfamily) with ens_labels.order_sf, so
DIRS/Ngaro/VIPER vote for the DIRS order and Penelope for PLE (the truth convention),
not for LTR / LINE by their RepeatMasker prefix. The order configuration is searched on
TRAIN with the same space and objective as tune_test.py.

Step 2 (Superfamily): only classifiers whose own call lies in the winning order vote,
each with its own superfamily weight (0 = no say). A superfamily is eligible only if it
belongs to that order: a label whose superfamily the truth taxonomy places in another
order (e.g. 'LINE/Gypsy') keeps its order vote but loses its superfamily vote. The top
superfamily is adopted if untied and carrying >= sf_frac of the superfamily weight of
the classifiers voting for that order (sf_frac 0 = plurality); otherwise the call stays
at order level. Searched on TRAIN for superfamily macro-F1 with step 1 fixed, so every
order-level call and metric is identical with or without step 2.

usage: two_step.py --kingdom plant --workdir ensemble_v2 --voter pipeline_head_v9/bin
"""
import argparse
import csv
import itertools
import json
import os
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tune_test as TT  # noqa: E402
from ens_labels import load_label_overrides, order_sf  # noqa: E402

MODELS, DISPLAY = TT.MODELS, TT.DISPLAY
ORDER_ONLY = {"LTR", "DNA", "LINE", "SINE", "RC"}      # truth 'superfamilies' that are just an order
NOT_SF = {"", "unknown", "unspecified", "unclassified", "na", "none", "other"}


def real_sf(sf):
    return sf if sf and sf.lower() not in NOT_SF and sf not in ORDER_ONLY else ""


def normalise(label, sf2order, dropped):
    """(order, superfamily) of one classifier call; superfamily '' if absent or outside its order."""
    if not label:
        return ("", "")
    o, sf = order_sf(label)
    if not o or o.lower() in NOT_SF:
        return ("", "")
    sf = real_sf(sf)
    if sf and sf2order.get(sf, o) != o:
        dropped[label] += 1
        sf = ""
    return (o, sf)


def decide_sf(order, norm, cfg):
    """Step 2: superfamily inside `order`, or None (call stays at order level)."""
    if not order:
        return None
    tot, w = 0.0, defaultdict(float)
    for m, (o, sf) in zip(MODELS, norm):
        wm = cfg["weights"].get(m, 0.0)
        if not wm or o != order:
            continue
        tot += wm
        if sf:
            w[sf] += wm
    if not w:
        return None
    ranked = sorted(w.values(), reverse=True)
    if len(ranked) > 1 and abs(ranked[0] - ranked[1]) < 1e-9:
        return None
    top = max(w.items(), key=lambda kv: kv[1])
    return top[0] if top[1] >= cfg["frac"] * tot - 1e-9 else None


def cm_metrics(pairs, classes):
    rows = defaultdict(Counter)
    for (t, p), n in pairs.items():
        rows[t][p] += n
    labs = sorted(set(rows) | {p for r in rows.values() for p in r})
    cm = np.array([[rows[t][p] for p in labs] for t in labs], dtype=float)
    return TT.metrics_from_cm(cm, labs, classes)


def eval_order(args):
    cfg, patterns, classes = args
    pairs = Counter()
    for otuple, tcount in patterns:
        p = TT.decide(otuple, cfg) or "Unknown"
        for t, n in tcount.items():
            pairs[(t, p)] += n
    return cfg, cm_metrics(pairs, classes)


def eval_sf(args):
    cfg, patterns, classes = args
    pairs = Counter()
    for (order, norm), tcount in patterns:
        p = decide_sf(order, norm, cfg) or "Unknown"
        for t, n in tcount.items():
            pairs[(t, p)] += n
    return cfg, cm_metrics(pairs, classes)


def sf_space():
    out = []
    for k in range(1, len(MODELS) + 1):
        for sub in itertools.combinations(MODELS, k):
            for ws in itertools.product((1.0, 2.0, 3.0), repeat=k):
                if min(ws) != 1.0:                       # scale-free: anchor the smallest weight at 1
                    continue
                for fr in (0.0, 0.5):
                    out.append({"weights": dict(zip(sub, ws)), "frac": fr})
    return out


def sf_name(cfg):
    w = ",".join("%s=%g" % (DISPLAY[m], v) for m, v in cfg["weights"].items())
    return "sf[%s] %s" % (w, "plurality" if cfg["frac"] == 0 else "frac=%.1f" % cfg["frac"])


def per_class_f1(truth, pred, classes):
    out = {}
    for c in classes:
        tp = sum(1 for t, p in zip(truth, pred) if t == c and p == c)
        fp = sum(1 for t, p in zip(truth, pred) if t != c and p == c)
        fn = sum(1 for t, p in zip(truth, pred) if t == c and p != c)
        out[c] = 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kingdom", required=True, choices=["plant", "animal", "fungi"])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--voter", required=True)
    ap.add_argument("--label-map", default=None)
    ap.add_argument("--min-support", type=int, default=20)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--n-boot", type=int, default=1000)
    a = ap.parse_args()
    load_label_overrides(a.label_map)
    TT.load_voter(a.voter)
    k, W = a.kingdom, a.workdir

    truth = {r["family_id"]: r for r in csv.DictReader(open("%s/pools/%s.truth.tsv" % (W, k)), delimiter="\t")}
    split = {r["family_id"]: r for r in csv.DictReader(open("%s/split/%s.split.tsv" % (W, k)), delimiter="\t")}
    vt = list(csv.DictReader(open("%s/preds/%s/voting_table.tsv" % (W, k)), delimiter="\t"))

    # canonical superfamily -> order from the curated truth of this kingdom
    so = defaultdict(Counter)
    for r in truth.values():
        sf = real_sf(r["superfamily"])
        if sf and r["order"]:
            so[sf][r["order"]] += 1
    sf2order = {sf: c.most_common(1)[0][0] for sf, c in so.items()}

    dropped = Counter()
    rows = []
    for r in vt:
        fid = r["repeat_id"]
        if fid not in truth or fid not in split:
            continue
        raw = tuple("" if r.get(DISPLAY[m], "NA") in ("NA", "") else r[DISPLAY[m]] for m in MODELS)
        norm = tuple(normalise(l, sf2order, dropped) for l in raw)
        rows.append({"raw": raw, "norm": norm, "order": truth[fid]["order"],
                     "sf": real_sf(truth[fid]["superfamily"]) or "Unknown",
                     "split": split[fid]["split"], "cluster": split[fid]["cluster"]})
    tr = [x for x in rows if x["split"] == "train"]
    te = [x for x in rows if x["split"] == "test"]
    n_calls = sum(1 for x in rows for l in x["raw"] if l)
    print("%s: %d rows (train %d, test %d); superfamily votes dropped as outside their order: %d of %d calls %s"
          % (k, len(rows), len(tr), len(te), sum(dropped.values()), n_calls, dict(dropped.most_common(8))), flush=True)

    emitted = {sf for x in rows for (_, sf) in x["norm"] if sf}
    def sf_classes_of(xs):
        c = Counter(x["sf"] for x in xs)
        return sorted(s for s, n in c.items() if n >= a.min_support and s != "Unknown" and s in emitted)
    def order_classes_of(xs):
        c = Counter(x["order"] for x in xs)
        return sorted(s for s, n in c.items() if n >= a.min_support)

    # ---- step 1: order search on TRAIN (normalised order votes) -------------------
    ocls = order_classes_of(tr)
    pat = defaultdict(Counter)
    for x in tr:
        pat[tuple(o for o, _ in x["norm"])][x["order"]] += 1
    cfgs = TT.search_space()
    with Pool(a.threads, initializer=TT.load_voter, initargs=(a.voter,)) as pool:
        res = pool.map(eval_order, [(c, list(pat.items()), ocls) for c in cfgs], chunksize=8)
    res.sort(key=lambda x: (-round(x[1]["macro_f1"], 4), -round(x[1]["mcc"], 4), len(x[0]["subset"])))
    ocfg = res[0][0]
    print("step 1 (order) best on train: %s  macro-F1 %.4f" % (TT.cfg_name(ocfg), res[0][1]["macro_f1"]), flush=True)

    def order_call(x):
        return TT.decide(tuple(o for o, _ in x["norm"]), ocfg) or ""

    # ---- step 2: superfamily search on TRAIN, step 1 fixed ------------------------
    scls = sf_classes_of(tr)
    spat = defaultdict(Counter)
    for x in tr:
        spat[(order_call(x), x["norm"])][x["sf"]] += 1
    scfgs = sf_space()
    with Pool(a.threads) as pool:
        sres = pool.map(eval_sf, [(c, list(spat.items()), scls) for c in scfgs], chunksize=16)
    sres.sort(key=lambda x: (-round(x[1]["macro_f1"], 4), -round(x[1]["accuracy"], 4), len(x[0]["weights"])))
    scfg = sres[0][0]
    print("step 2 (superfamily) best on train: %s  sf macro-F1 %.4f  (%d configs, classes %s)"
          % (sf_name(scfg), sres[0][1]["macro_f1"], len(scfgs), scls), flush=True)
    os.makedirs("%s/tune" % W, exist_ok=True)
    json.dump({"order": ocfg, "superfamily": scfg}, open("%s/tune/%s.two_step.best.json" % (W, k), "w"), indent=1)
    with open("%s/tune/%s.two_step.top_sf.tsv" % (W, k), "w") as fh:
        fh.write("rank\tconfig\tsf_macro_f1\tsf_accuracy\tsf_coverage\tsf_precision_called\n")
        for i, (c, m) in enumerate(sres[:50], 1):
            fh.write("%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n" % (i, sf_name(c), m["macro_f1"], m["accuracy"],
                                                         m["coverage"], m["precision_called"]))
    # best single classifier at superfamily level, chosen on TRAIN
    best_single = max(MODELS, key=lambda m: eval_sf(({"weights": {m: 1.0}, "frac": 0.0},
                                                      list(spat.items()), scls))[1]["macro_f1"])

    # ---- TEST ---------------------------------------------------------------------
    t_order = [x["order"] for x in te]
    t_sf = [x["sf"] for x in te]
    clusters = np.array([x["cluster"] for x in te])
    ocls_te, scls_te = order_classes_of(te), sf_classes_of(te)
    old_best = json.load(open("%s/tune/%s.best.json" % (W, k)))
    P_o, P_s = {}, {}
    for m in MODELS:
        i = MODELS.index(m)
        P_o[DISPLAY[m]] = [x["norm"][i][0] or "Unknown" for x in te]
        P_s[DISPLAY[m]] = [x["norm"][i][1] or "Unknown" for x in te]
    old = [TT.decide(x["raw"], old_best) or "" for x in te]
    P_o["tuned, current voter"] = [TT.to_order(l) for l in old]
    P_s["tuned, current voter"] = [real_sf(order_sf(l)[1]) or "Unknown" if l else "Unknown" for l in old]
    ocalls = [order_call(x) for x in te]
    variants = [("two-step, sf by %s alone" % DISPLAY[best_single], {"weights": {best_single: 1.0}, "frac": 0.0}),
                ("two-step, sf all five equal", {"weights": {m: 1.0 for m in MODELS}, "frac": 0.0}),
                ("two-step tuned", scfg)]
    for name, c in variants:
        P_o[name] = [o or "Unknown" for o in ocalls]
        P_s[name] = [decide_sf(o, x["norm"], c) or "Unknown" for o, x in zip(ocalls, te)]
    ci_o = TT.bootstrap_ci(t_order, P_o, clusters, ocls_te, n_boot=a.n_boot)
    ci_s = TT.bootstrap_ci(t_sf, P_s, clusters, scls_te, n_boot=a.n_boot)

    os.makedirs("%s/test" % W, exist_ok=True)
    out = "%s/test/%s.two_step.tsv" % (W, k)
    with open(out, "w") as fh:
        fh.write("method\torder_f1\torder_f1_ci\tsf_f1\tsf_f1_ci\tsf_accuracy\tsf_coverage\tsf_precision_called\n")
        for m in P_o:
            ro, rs = TT.metrics(t_order, P_o[m], ocls_te), TT.metrics(t_sf, P_s[m], scls_te)
            fh.write("%s\t%.4f\t%.3f-%.3f\t%.4f\t%.3f-%.3f\t%.4f\t%.4f\t%.4f\n" % (
                m, ro["macro_f1"], *ci_o[m]["macro_f1"], rs["macro_f1"], *ci_s[m]["macro_f1"],
                rs["accuracy"], rs["coverage"], rs["precision_called"]))
        for ref in (DISPLAY[best_single], "tuned, current voter"):
            n01, n10, p = TT.mcnemar(t_sf, P_s["two-step tuned"], P_s[ref])
            fh.write("# McNemar (superfamily) two-step tuned vs %s: two-step-only correct %d, %s-only correct %d, p=%.2e\n"
                     % (ref, n01, ref, n10, p))
        fh.write("# step 1 order config: %s\n# step 2 superfamily config: %s\n" % (TT.cfg_name(ocfg), sf_name(scfg)))
        fh.write("# tuned order config of tune_test.py (current voter): %s\n" % TT.cfg_name(old_best))
        fh.write("# order classes: %s\n# superfamily classes: %s\n" % (",".join(ocls_te), ",".join(scls_te)))
        fh.write("# per-superfamily F1 (test):\n")
        show = [DISPLAY[best_single], "tuned, current voter", "two-step tuned"]
        pcs = {m: per_class_f1(t_sf, P_s[m], scls_te) for m in show}
        fh.write("#   %-16s %s\n" % ("superfamily", "  ".join("%22s" % m for m in show)))
        n_te = Counter(t_sf)
        for c in scls_te:
            fh.write("#   %-16s %s\n" % ("%s(%d)" % (c, n_te[c]), "  ".join("%22.2f" % pcs[m][c] for m in show)))
    print(open(out).read(), flush=True)


if __name__ == "__main__":
    main()

