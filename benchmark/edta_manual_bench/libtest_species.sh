#!/bin/bash
# libtest_species.sh <species>... - lib-test.py for every method vs the reference, then the CSV.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
for SP in "$@"; do
    KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
    G=$B/genomes/$KING/$SP/$SP.fna
    REF=$E/$SP/reference/$SP.fna.out
    for M in teannotate_ref9 teannotate_wo_ref9 teannotate_ref8 teannotate_wo_ref8 teannotate_ref7 teannotate_wo_ref7 teannotate_ref6 teannotate_wo_ref6 teannotate_ref5 teannotate_wo_ref5 teannotate teannotate_wo teannotate_ref teannotate_wo_ref teannotate_ref2 teannotate_wo_ref2 teannotate_ref3 teannotate_wo_ref3 teannotate_ref4 teannotate_wo_ref4 hite edta edta2 heliano earlgrey; do
        T=$E/$SP/$M/$SP.fna.out
        [[ -s $T ]] || continue
        echo "[$(date)] lib-test $SP/$M"
        ( cd "$E/$SP/$M" && $PY $B/tools/EDTA/lib-test.py --genome "$G" --reference "$REF" --test "$T" --extended_report ) \
            > "$E/$SP/$M/libtest.log" 2>&1 || echo "  lib-test FAILED for $SP/$M (see libtest.log)"
    done
done
$PY - "$@" <<'PY'
import os, re, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
label = {"teannotate_ref9": "teannotate with tetrimmer (v9)", "teannotate_wo_ref9": "teannotate w/o tetrimmer (v9)",
         "teannotate_ref8": "teannotate with tetrimmer (v8)", "teannotate_wo_ref8": "teannotate w/o tetrimmer (v8)",
         "teannotate_ref7": "teannotate with tetrimmer (v7)", "teannotate_wo_ref7": "teannotate w/o tetrimmer (v7)",
         "teannotate_ref6": "teannotate with tetrimmer (v6)", "teannotate_wo_ref6": "teannotate w/o tetrimmer (v6)",
         "teannotate_ref5": "teannotate with tetrimmer", "teannotate_wo_ref5": "teannotate w/o tetrimmer",
         "teannotate_ref4": "teannotate with tetrimmer (v4)", "teannotate_wo_ref4": "teannotate w/o tetrimmer (v4)",
         "teannotate_ref3": "teannotate with tetrimmer (refine v3)", "teannotate_wo_ref3": "teannotate w/o tetrimmer (refine v3)",
         "teannotate_ref2": "teannotate with tetrimmer (refine v2)", "teannotate_wo_ref2": "teannotate w/o tetrimmer (refine v2)",
         "edta2": "edta", "edta": "edta (rice fallback)",
         "teannotate_ref": "teannotate with tetrimmer (refine v1)", "teannotate_wo_ref": "teannotate w/o tetrimmer (refine v1)",
         "teannotate": "teannotate with tetrimmer (old dedup)", "teannotate_wo": "teannotate w/o tetrimmer (old dedup)"}
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    rows = []
    for m in ["teannotate_ref9", "teannotate_wo_ref9", "teannotate_ref8", "teannotate_wo_ref8", "teannotate_ref7", "teannotate_wo_ref7", "teannotate_ref6", "teannotate_wo_ref6", "teannotate_ref5", "teannotate_wo_ref5", "teannotate_ref4", "teannotate_wo_ref4", "teannotate_ref3", "teannotate_wo_ref3", "teannotate_ref2", "teannotate_wo_ref2", "teannotate_ref", "teannotate_wo_ref", "teannotate", "teannotate_wo", "hite", "edta2", "edta", "heliano", "earlgrey"]:
        rep = f"{E}/{sp}/{m}/{sp}.fna.out.repeatmasker.lib.report"
        if not os.path.exists(rep):
            continue
        for line in open(rep):
            f = line.rstrip("\n").split("\t")
            mm = re.search(r"\.(\w+)\.lib\.report$", f[0])
            if mm and len(f) >= 11:
                rows.append([label.get(m, m), mm.group(1)] + f[1:11])
    out = f"{E}/results/{nice}_bench_edta_manual.csv"
    with open(out, "w") as fh:
        fh.write("Method,TE type,sens,spec,accu,prec,FDR,F1,TP,TN,FP,FN\n")
        for r in rows:
            fh.write(",".join(r) + "\n")
    print(out, len(rows), "rows")
PY
