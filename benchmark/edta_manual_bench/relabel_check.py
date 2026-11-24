"""Did our pipeline relabel EarlGrey consensi, and are the relabeled ones junk?
Compares EarlGrey's own label with the label the same consensus carries in teannotate."""
import collections, re
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
for sp in ("drosophila_melanogaster", "zymoseptoria_tritici", "septoria_linicola"):
    eg = {}
    for l in open(f"{E}/{sp}/earlgrey/library.fa"):
        if l.startswith(">"):
            n, _, lab = l[1:].split()[0].partition("#")
            eg[n] = lab
    feats = [l.rstrip("\n").split("\t") for l in open(f"{E}/junk_analysis/{sp}.features.tsv")]
    hdr = feats[0]; rows = [dict(zip(hdr, f)) for f in feats[1:]]
    agg = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        if r["source"] != "EarlGrey": continue
        base = re.sub(r"_part\d+$", "", r["name"])
        if base not in eg: continue
        a, b = eg[base].split("/")[0], r["label"].split("/")[0]
        key = "kept " + a if a == b else f"{a} -> {b}"
        agg[key][0] += 1; agg[key][1] += int(r["tp"]); agg[key][2] += int(r["fp"])
    print(f"\n{sp}: EarlGrey consensi in the teannotate library, by label change")
    for k, v in sorted(agg.items(), key=lambda x: -x[1][2]):
        print(f"  {k:28s} n={v[0]:3d}  TP {v[1]/1e6:5.2f} Mb  FP {v[2]/1e6:5.2f} Mb  FP share {v[2]/max(1,v[1]+v[2]):.2f}")
