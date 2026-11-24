#!/bin/bash
# lib-test.py (EDTA) for every test vs its reference, then your CSV format per species.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
for SP in arabidopsis_thaliana drosophila_melanogaster; do
    KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
    G=$B/genomes/$KING/$SP/$SP.fna
    REF=$E/$SP/reference/$SP.fna.out
    for M in teannotate hite edta heliano earlgrey; do
        T=$E/$SP/$M/$SP.fna.out
        [[ -s $T ]] || continue
        echo "[$(date)] lib-test $SP/$M"
        ( cd "$E/$SP/$M" && $PY $B/tools/EDTA/lib-test.py --genome "$G" --reference "$REF" --test "$T" --extended_report ) \
            > "$E/$SP/$M/libtest.log" 2>&1 || echo "  lib-test FAILED for $SP/$M (see libtest.log)"
    done
done
# CSV in the user's format: Method,TE type,sens,spec,accu,prec,FDR,F1,TP,TN,FP,FN
$PY - <<'PY'
import glob, os, re
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
names = {"arabidopsis_thaliana": "Arabidopsis_thaliana", "drosophila_melanogaster": "Drosophila_melanogaster"}
order = ["teannotate", "hite", "edta", "heliano", "earlgrey"]
for sp, nice in names.items():
    rows = []
    for m in order:
        rep = f"{E}/{sp}/{m}/{sp}.fna.out.repeatmasker.lib.report"
        if not os.path.exists(rep):
            continue
        for line in open(rep):
            f = line.rstrip("\n").split("\t")
            mm = re.search(r"\.(\w+)\.lib\.report$", f[0])
            if mm and len(f) >= 11:
                rows.append([m, mm.group(1)] + f[1:11])
    out = f"{E}/results/{nice}_bench_edta_manual.csv"
    with open(out, "w") as fh:
        fh.write("Method,TE type,sens,spec,accu,prec,FDR,F1,TP,TN,FP,FN\n")
        for r in rows:
            fh.write(",".join(r) + "\n")
    print(out, len(rows), "rows")
PY
