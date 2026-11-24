"""TEsorter domain evidence vs junk, and domain class vs our label vs the reference class.
reference class of a consensus = the lib-test category holding most of its TP bp
(tp_ltr / tp_nonltr / tp_tir / tp_helitron from sequence_assessment.report)."""
import collections, sys, re
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
M = "teannotate_wo_ref4"

def cat(label):
    l = label.lower()
    if re.search(r"helitron|^rc", l): return "helitron"
    if re.search(r"ltr|gypsy|copia|bel|pao|erv|retrovirus|dirs|ngaro", l) and "non" not in l: return "ltr"
    if re.search(r"line|sine|penelope|nonltr|^l1|jockey|cr1|i-jockey|rte|tad1", l): return "nonltr"
    if re.search(r"tir|dna|mite|hat|mule|mutator|cacta|enspm|pif|harbinger|tc1|mariner|piggybac|transib|^p$|merlin", l): return "tir"
    return "other"

for sp in sys.argv[1:]:
    T = f"{E}/junk_analysis/tesorter/{sp}/{sp}.cls.tsv"
    dom = {}
    for l in open(T):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        name = f[0].split("#")[0]
        order, sf = f[1], f[2]
        dom[name] = ("helitron" if order == "Helitron" else "ltr" if order == "LTR" else
                     "nonltr" if order in ("LINE", "SINE", "pararetrovirus") else
                     "tir" if order in ("TIR", "Maverick", "Crypton") else "other", order + "/" + sf)
    rep = open(f"{E}/{sp}/{M}/{sp}.fna.out.sequence_assessment.report").read().splitlines()
    hdr = rep[0].split("\t")
    feats = {}
    for l in open(f"{E}/junk_analysis/{sp}.features.tsv"):
        f = l.rstrip("\n").split("\t")
        if f[0] == "species": fh = f; continue
        d = dict(zip(fh, f)); feats[d["name"]] = d
    junk_fp = {True: 0, False: 0}; real_tp = {True: 0, False: 0}
    lab_ok = collections.Counter(); confusion = collections.Counter()
    for l in rep[1:]:
        d = dict(zip(hdr, l.split("\t")))
        n = d["sequence_id"]
        if n not in feats: continue
        tp, fp = int(d["tp_all_cats"]), int(d["fp_all_cats"])
        has = n in dom
        if fp / max(1, tp + fp) > 0.5: junk_fp[has] += fp
        else: real_tp[has] += tp
        # reference class by TP bp
        byc = {c: int(d[f"tp_{c}"]) for c in ("ltr", "nonltr", "tir", "helitron")}
        refc, refbp = max(byc.items(), key=lambda x: x[1])
        if refbp < 5000 or not has: continue
        ours, theirs = cat(feats[n]["label"]), dom[n][0]
        if ours == theirs: lab_ok["agree"] += refbp
        else:
            confusion[(feats[n]["label"].split("/")[0], dom[n][1], refc)] += refbp
            lab_ok["ours right" if ours == refc else "domain right" if theirs == refc else "both wrong"] += refbp
    J = junk_fp[True] + junk_fp[False]; G = real_tp[True] + real_tp[False]
    print(f"\n{sp}: {len(dom)} consensi with TE domains")
    print(f"  junk FP bp in consensi WITH domains: {junk_fp[True]/1e6:.2f} Mb ({junk_fp[True]/max(1,J):.0%}); WITHOUT: {junk_fp[False]/1e6:.2f} Mb")
    print(f"  real TP bp in consensi WITH domains: {real_tp[True]/1e6:.2f} Mb ({real_tp[True]/max(1,G):.0%}); WITHOUT: {real_tp[False]/1e6:.2f} Mb")
    tot = sum(lab_ok.values())
    print("  label check on consensi with domains (bp of TP, by reference class):",
          ", ".join(f"{k} {v/1e6:.2f} Mb ({v/max(1,tot):.0%})" for k, v in lab_ok.most_common()))
    print("  disagreements (our label | TEsorter | reference class) by TP bp:")
    for k, v in confusion.most_common(8):
        print(f"    {k[0]:10s} | {k[1]:28s} | {k[2]:8s} {v/1e6:.2f} Mb")
