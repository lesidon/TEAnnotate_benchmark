#!/usr/bin/env python3
"""verify_voter.py <pipeline_bin> - does the pipeline's two-step voter (bin/ensemble_vote.py)
reproduce the benchmark's two-step decisions (ensemble_v2 two_step.py, weighted-F1 tuning)?
Compares the calls on every TEST consensus and recomputes the headline metrics."""
import csv, json, sys
from collections import Counter

X = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v2"
sys.path.insert(0, X + "/scripts")
import tune_test as TT                     # noqa: E402
import two_step as TS                      # noqa: E402
from ens_labels import load_label_overrides, order_sf  # noqa: E402

sys.path.insert(0, sys.argv[1])
import ensemble_vote as PV                 # noqa: E402  (pipeline voter under test)

load_label_overrides(X + "/scripts/label_map.json")
TT.load_voter(X + "/../pipeline_head_v9/bin")   # benchmark step 1 uses the v9 voter's weighted rule
KMAP = {"plant": "plant", "animal": "animal", "fungi": "fungus"}

for k, pk in KMAP.items():
    truth = {r["family_id"]: r for r in csv.DictReader(open(f"{X}/pools/{k}.truth.tsv"), delimiter="\t")}
    split = {r["family_id"]: r for r in csv.DictReader(open(f"{X}/split/{k}.split.tsv"), delimiter="\t")}
    best = json.load(open(f"{X}/tune/{k}.wf1.two_step.best.json"))
    so = {}
    for r in truth.values():
        r["sf"] = TS.real_sf(order_sf(r["label_raw"])[1]) if r["label_raw"] else ""
    for r in truth.values():
        if r["sf"] and r["order"]:
            so.setdefault(r["sf"], Counter())[r["order"]] += 1
    sf2order = {s: c.most_common(1)[0][0] for s, c in so.items()}
    cfg = PV.load_two_step(pk)
    ran = list(TT.MODELS)
    te, diff, ex = [], Counter(), []
    for r in csv.DictReader(open(f"{X}/preds/{k}/voting_table.tsv"), delimiter="\t"):
        fid = r["repeat_id"]
        if fid not in truth or split.get(fid, {}).get("split") != "test":
            continue
        raw = tuple("" if r.get(TT.DISPLAY[m], "NA") in ("NA", "") else r[TT.DISPLAY[m]] for m in TT.MODELS)
        norm = tuple(TS.normalise(l, sf2order, Counter()) for l in raw)
        # benchmark decision
        o = TT.decide(tuple(x for x, _ in norm), best["order"]) or ""
        s = TS.decide_sf(o, norm, best["superfamily"]) or ""
        # pipeline decision
        call = PV.decide_two_step({m: l for m, l in zip(TT.MODELS, raw) if l}, cfg, ran) or ""
        po = TT.to_order(call) if call else ""
        ps = TS.real_sf(order_sf(call)[1]) if call else ""
        if (po, ps) != (o, s):
            diff[("order" if po != o else "superfamily")] += 1
            if len(ex) < 5:
                ex.append((fid, raw, (o, s), call))
        te.append((truth[fid]["order"], truth[fid]["sf"], po or "Unknown", ps or "Unknown"))
    t_o = [x[0] for x in te]; p_o = [x[2] for x in te]
    ocls = sorted(c for c, n in Counter(t_o).items() if n >= 20)
    keep = [x for x in te if x[1]]
    t_s = [x[1] for x in keep]; p_s = [x[3] for x in keep]
    ref = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{X}/test/{k}.wf1.two_step.tsv") if not l.startswith("#")}
    sf_classes = [c for c in open(f"{X}/test/{k}.wf1.two_step.tsv").read().split("# superfamily classes: ")[1].split("\n")[0].split(",")]
    mo, ms = TT.metrics(t_o, p_o, ocls), TT.metrics(t_s, p_s, sf_classes)
    r = ref["two-step tuned"]
    print(f"{k}: {len(te)} test consensi; calls differing from the benchmark: {dict(diff) or 'none'}")
    print(f"   order wF1 pipeline {mo['weighted_f1']:.4f} vs benchmark {float(r[1]):.4f};"
          f"  superfamily wF1 pipeline {ms['weighted_f1']:.4f} vs benchmark {float(r[6]):.4f}")
    for e in ex:
        print("   e.g.", e)

