"""Tiered evidence: E1 TE protein domain (TEsorter/REXdb); E2 no domain but another tool
found the same element (>=50% at >=80% id); E3 neither, but >=2 full-length genomic copies;
E4 none of these. Junk FP / real TP per tier, and estimated lib-test 'any' if E4 is dropped."""
import csv, sys, collections
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
best = {"arabidopsis_thaliana": "HiTE 0.81", "drosophila_melanogaster": "HiTE 0.88",
        "zymoseptoria_tritici": "EarlGrey 0.86", "septoria_linicola": "HiTE 0.91"}
NONTE = ("Satellite", "rRNA", "tRNA", "snRNA", "Simple", "Low_complexity", "Cent")
for sp in sys.argv[1:]:
    dom = set()
    for l in open(f"{E}/junk_analysis/tesorter/{sp}/{sp}.cls.tsv"):
        if not l.startswith("#"): dom.add(l.split("\t")[0].split("#")[0])
    rows = list(csv.DictReader(open(f"{E}/junk_analysis/{sp}.features.tsv"), delimiter="\t"))
    def tier(r):
        if r["label"].startswith(NONTE): return "E0 non-TE label"
        if r["name"] in dom: return "E1 domain"
        if int(r["n_support"]) >= 1: return "E2 cross-tool"
        if int(r["n_full"]) >= 2: return "E3 >=2 full copies"
        return "E4 no evidence"
    agg = collections.defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        tp, fp = int(r["tp"]), int(r["fp"])
        t = tier(r); a = agg[t]; a[0] += 1
        if fp / max(1, tp + fp) > 0.5: a[1] += fp
        else: a[2] += tp
        a[3] += tp
    print(f"\n{sp}")
    for t in sorted(agg):
        a = agg[t]
        print(f"  {t:20s} n={a[0]:4d}  junk FP {a[1]/1e6:5.2f} Mb   real TP {a[2]/1e6:5.2f} Mb")
    # lib-test 'any' totals for the w/o library
    nice = sp[0].upper() + sp[1:]
    lt = [x for x in csv.DictReader(open(f"{E}/results/{nice}_bench_edta_manual.csv"))
          if x["Method"] == "teannotate w/o tetrimmer" and x["TE type"] == "any"][0]
    TP, FP, FN = float(lt["TP"]), float(lt["FP"]), float(lt["FN"])
    stp = sum(int(r["tp"]) for r in rows); sfp = sum(int(r["fp"]) for r in rows)
    for drop in (("E0 non-TE label", "E4 no evidence"), ("E0 non-TE label", "E3 >=2 full copies", "E4 no evidence")):
        d = [r for r in rows if tier(r) in drop]
        dtp = sum(int(r["tp"]) for r in d) / stp; dfp = sum(int(r["fp"]) for r in d) / sfp
        tp, fp, fn = TP * (1 - dtp), FP * (1 - dfp), FN + TP * dtp
        s, p = tp / (tp + fn), tp / (tp + fp)
        print(f"  drop {'+'.join(x.split()[0] for x in drop):9s}: {len(d):4d} seqs -> sens {s:.2f} prec {p:.2f} F1 {2*s*p/(s+p):.2f}"
              f"  (now {2*TP/(2*TP+FP+FN):.2f}; best standalone {best[sp]})")
