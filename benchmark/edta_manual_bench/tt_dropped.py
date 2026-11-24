"""What does TEtrimmer drop (skipped and not written back), and how much of it is real?
Joins TEtrimmer's summary.txt with the w/o-library features (same input names) and
TEsorter domains; TP/FP come from the w/o library's lib-test (approximate)."""
import csv, sys, glob, collections
B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E = f"{B}/edta_manual_bench"
NONTE = ("Satellite", "rRNA", "tRNA", "snRNA", "Simple", "Low_complexity", "Cent")
for sp in sys.argv[1:]:
    summ = None
    for cmd in glob.glob(f"{B}/work/refined_dedup_v4/{sp}/work/*/*/.command.sh"):
        t = open(cmd).read()
        if "TETRIMMER_CMD" in t and "--input_file" in t:
            cand = glob.glob(cmd.replace(".command.sh", "*_tetrimmer/summary.txt"))
            if cand: summ = cand[0]; break
    if not summ:
        print(sp, "no TEtrimmer summary found"); continue
    out_dir = summ.rsplit("/", 1)[0]
    written = set(l[1:].split("#")[0].strip() for l in open(f"{out_dir}/TEtrimmer_consensus.fasta") if l.startswith(">"))
    status = {r["input_name"]: r["status"] for r in csv.DictReader(open(summ))}
    dom = set(l.split("\t")[0].split("#")[0] for l in open(f"{E}/junk_analysis/tesorter/{sp}/{sp}.cls.tsv") if not l.startswith("#"))
    feats = {r["name"]: r for r in csv.DictReader(open(f"{E}/junk_analysis/{sp}.features.tsv"), delimiter="\t")}
    def tier(n, f):
        if f["label"].startswith(NONTE): return "E0 non-TE"
        if n in dom: return "E1 domain"
        if int(f["n_support"]) >= 1: return "E2 cross-tool"
        if int(f["n_full"]) >= 2: return "E3 full copies"
        return "E4 none"
    agg = collections.defaultdict(lambda: [0, 0, 0])
    for n, st in status.items():
        if st != "skipped" or n in written or n not in feats: continue
        f = feats[n]; t = tier(n, f)
        agg[t][0] += 1; agg[t][1] += int(f["tp"]); agg[t][2] += int(f["fp"])
    tot_tp = sum(int(f["tp"]) for f in feats.values())
    print(f"\n{sp}: TEtrimmer dropped {sum(v[0] for v in agg.values())} inputs (skipped, not written back)")
    for t in sorted(agg):
        v = agg[t]
        print(f"  {t:15s} n={v[0]:4d}  TP {v[1]/1e6:5.2f} Mb ({v[1]/tot_tp:4.0%} of library TP)  FP {v[2]/1e6:5.2f} Mb")
