#!/usr/bin/env python3
"""Paired cluster-bootstrap CI of the support-weighted F1 difference between the tuned two-step
vote and the best single classifier on the test half (order and superfamily level).

Rebuilds the test predictions exactly as two_step.py does, from the saved tuned configurations
(tune/<k>.wf1.two_step.best.json), without re-tuning. Each bootstrap draw resamples whole
clusters once and scores both methods on that same draw, so the CI is of the paired difference.

usage: paired_boot.py <workdir> <voter bin dir> [n_boot] [out.tsv] [tune dir]
  workdir   holds pools/, split/, preds/ and scripts/ (ensemble_v2 layout)
  tune dir  where <k>.wf1.two_step.best.json is read (default <workdir>/tune); pass the tune dir of
            another run to score fixed (e.g. the deployed) settings on this workdir's consensi
"""
import csv
import json
import sys
from collections import Counter

import numpy as np

W, VOTER = sys.argv[1], sys.argv[2]
N_BOOT = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
TUNE = sys.argv[5] if len(sys.argv) > 5 else W + "/tune"
sys.path.insert(0, W + "/scripts")
import tune_test as TT  # noqa: E402
import two_step as TS  # noqa: E402
from ens_labels import load_label_overrides, order_sf  # noqa: E402

load_label_overrides(W + "/scripts/label_map.json")
TT.load_voter(VOTER)
MODELS, DISPLAY = TT.MODELS, TT.DISPLAY


def wf1_from_cm(cm, labs, classes):
    return TT.metrics_from_cm(cm, labs, classes)["weighted_f1"]


def paired(truth, pa, pb, clusters, classes, n_boot, seed=0):
    labs = sorted(set(truth) | set(pa) | set(pb))
    li = {l: i for i, l in enumerate(labs)}
    L = len(labs)
    t = np.array([li[x] for x in truth])
    a = np.array([li[x] for x in pa])
    b = np.array([li[x] for x in pb])
    uniq, inv = np.unique(clusters, return_inverse=True)
    C = len(uniq)
    rng = np.random.default_rng(seed)
    full = lambda p, w: wf1_from_cm(np.bincount(t * L + p, weights=w, minlength=L * L).reshape(L, L), labs, classes)
    one = np.ones(len(t))
    obs = full(a, one) - full(b, one)
    d = np.empty(n_boot)
    for i in range(n_boot):
        w = np.bincount(rng.integers(0, C, C), minlength=C)[inv].astype(float)
        d[i] = full(a, w) - full(b, w)
    return full(a, one), full(b, one), obs, np.percentile(d, 2.5), np.percentile(d, 97.5), float((d <= 0).mean()), C


out = []
for k in ("plant", "animal", "fungi"):
    truth = {r["family_id"]: r for r in csv.DictReader(open("%s/pools/%s.truth.tsv" % (W, k)), delimiter="\t")}
    split = {r["family_id"]: r for r in csv.DictReader(open("%s/split/%s.split.tsv" % (W, k)), delimiter="\t")}
    vt = list(csv.DictReader(open("%s/preds/%s/voting_table.tsv" % (W, k)), delimiter="\t"))
    best = json.load(open("%s/%s.wf1.two_step.best.json" % (TUNE, k)))
    ocfg, scfg = best["order"], best["superfamily"]
    for r in truth.values():
        r["sf"] = TS.real_sf(order_sf(r["label_raw"])[1]) if r["label_raw"] else ""
    so = {}
    cnt = {}
    for r in truth.values():
        if r["sf"] and r["order"]:
            cnt.setdefault(r["sf"], Counter())[r["order"]] += 1
    sf2order = {sf: c.most_common(1)[0][0] for sf, c in cnt.items()}
    dropped = Counter()
    te = []
    all_rows = []
    for r in vt:
        fid = r["repeat_id"]
        if fid not in truth or fid not in split:
            continue
        raw = tuple("" if r.get(DISPLAY[m], "NA") in ("NA", "") else r[DISPLAY[m]] for m in MODELS)
        norm = tuple(TS.normalise(l, sf2order, dropped) for l in raw)
        x = {"norm": norm, "order": truth[fid]["order"], "sf": truth[fid]["sf"],
             "split": split[fid]["split"], "cluster": split[fid]["cluster"]}
        all_rows.append(x)
        if x["split"] == "test":
            te.append(x)
    emitted = {sf for x in all_rows for (_, sf) in x["norm"] if sf}
    keep = [i for i, x in enumerate(te) if x["sf"]]
    t_order = [x["order"] for x in te]
    t_sf = [te[i]["sf"] for i in keep]
    oc = Counter(t_order)
    ocls = sorted(s for s, n in oc.items() if n >= 20)
    sc = Counter(t_sf)
    scls = sorted(s for s, n in sc.items() if n >= 20 and s and s in emitted)
    clusters = np.array([x["cluster"] for x in te])
    ocalls = [TT.decide(tuple(o for o, _ in x["norm"]), ocfg) or "" for x in te]
    vote_o = [o or "Unknown" for o in ocalls]
    vote_s = [TS.decide_sf(o, x["norm"], scfg) or "Unknown" for o, x in zip(ocalls, te)]
    vote_s = [vote_s[i] for i in keep]
    for level, truth_l, vote, cl, classes in (
            ("order", t_order, vote_o, clusters, ocls),
            ("superfamily", t_sf, vote_s, clusters[keep], scls)):
        singles = {}
        for m in MODELS:
            j = MODELS.index(m)
            if level == "order":
                singles[m] = [x["norm"][j][0] or "Unknown" for x in te]
            else:
                singles[m] = [te[i]["norm"][j][1] or "Unknown" for i in keep]
        sc_ = {m: TT.metrics(truth_l, p, classes)["weighted_f1"] for m, p in singles.items()}
        bm = max(sc_, key=sc_.get)
        fa, fb, obs, lo, hi, p, C = paired(truth_l, vote, singles[bm], cl, classes, N_BOOT)
        line = "%s\t%s\t%d\t%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4g" % (
            k, level, len(truth_l), C, DISPLAY[bm], fa, fb, obs, lo, hi, p)
        print(line, flush=True)
        out.append(line)
with open(sys.argv[4] if len(sys.argv) > 4 else "paired_boot.tsv", "w") as fh:
    fh.write("kingdom\tlevel\tn_test\tn_clusters\tbest_single\tvote_wf1\tsingle_wf1\tdelta\tci_lo\tci_hi\tfrac_boot_delta_le0\n")
    fh.write("\n".join(out) + "\n")
