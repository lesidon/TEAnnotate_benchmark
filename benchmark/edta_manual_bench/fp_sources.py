"""Where do a library's false-positive bp come from? Uses lib-test's per-sequence
sequence_assessment.report (tp_all_cats / fp_all_cats per consensus) + library headers."""
import re, sys, collections

E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"

def source(name):
    n = name.replace("_", "-")
    if re.match(r"TE-\d+", n): return "EDTA"
    if re.match(r"(LTR|TIR|Helitron|Homology-Non-LTR|Denovo-Non-LTR)-\d+", n): return "HiTE"
    if n.startswith("insertion-HLE"): return "HELIANO"
    if n.startswith("rnd-"): return "EarlGrey"
    return "other"

def lib_meta(path):
    meta = {}
    name = None
    for l in open(path):
        if l.startswith(">"):
            h = l[1:].split()[0]
            name, _, label = h.partition("#")
            meta[name] = [label or "NA", 0]
        elif name:
            meta[name][1] += len(l.strip())
    return meta

for sp in sys.argv[1:]:
    for m in ("teannotate_wo_ref4", "teannotate_ref4", "hite", "earlgrey"):
        D = f"{E}/{sp}/{m}"
        try:
            meta = lib_meta(f"{D}/library.fa")
            rep = open(f"{D}/{sp}.fna.out.sequence_assessment.report").read().splitlines()
        except FileNotFoundError:
            continue
        hdr = rep[0].split("\t")
        rows = []
        for l in rep[1:]:
            f = dict(zip(hdr, l.split("\t")))
            name = f["sequence_id"]
            lab, ln = meta.get(name, ["?", 0])
            rows.append((name, lab, ln, int(f["number_of_matches"]), int(f["tp_all_cats"]), int(f["fp_all_cats"])))
        TP = sum(r[4] for r in rows); FP = sum(r[5] for r in rows)
        print(f"\n######## {sp} / {m}: {len(meta)} lib seqs, {len(rows)} with hits; TP {TP/1e6:.2f} Mb, FP {FP/1e6:.2f} Mb")
        by = collections.defaultdict(lambda: [0, 0, 0])
        for r in rows:
            k = (source(r[0]), r[1].split("/")[0])
            by[k][0] += 1; by[k][1] += r[4]; by[k][2] += r[5]
        print("  source/class        n   TP Mb  FP Mb  FP/(TP+FP)")
        for k, v in sorted(by.items(), key=lambda x: -x[1][2])[:12]:
            print(f"  {k[0]:9s}{k[1]:10s} {v[0]:4d} {v[1]/1e6:6.2f} {v[2]/1e6:6.2f}   {v[2]/max(1,v[1]+v[2]):.2f}")
        for thr in (0.9, 0.75, 0.5):
            sel = [r for r in rows if r[5] / max(1, r[4] + r[5]) > thr]
            print(f"  seqs with FP share >{thr}: {len(sel):4d}  carry FP {sum(r[5] for r in sel)/1e6:.2f} Mb "
                  f"({sum(r[5] for r in sel)/max(1,FP):.0%} of FP), TP {sum(r[4] for r in sel)/1e6:.2f} Mb ({sum(r[4] for r in sel)/max(1,TP):.0%} of TP)")
        print("  top FP consensi: name | label | len | hits | TP | FP")
        for r in sorted(rows, key=lambda r: -r[5])[:10]:
            print(f"    {r[0][:34]:34s} {r[1][:18]:18s} {r[2]:6d} {r[3]:5d} {r[4]:8d} {r[5]:8d}")
