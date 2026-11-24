#!/bin/bash
# thr_bench.sh - copy-threshold re-check: build the gate variants of the v10 w/o-TEtrimmer
# libraries (thr_variants.py) and score them like the dev bench: specs row -> rm_one.sh ->
# harmonised lib-test (+ any_dontcare) -> family level. Results go to
# results/<Sp>_bench_harmonised_thr.csv (a copy of the scorer), next to the default
# (teannotate_wo_ref10 = F2/H5) and the standalone tools; shared result files are untouched.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; T=$B/thresholds
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
SPECIES="zymoseptoria_tritici septoria_linicola drosophila_melanogaster arabidopsis_thaliana"
COMPARE="teannotate_wo_ref10 teannotate_ref10 earlgrey hite edta2"
$PY $T/thr_variants.py $WORK_DIR/refined_dedup_v10 $T/benchlibs $SPECIES
sed 's/_bench_harmonised.csv/_bench_harmonised_thr.csv/' $E/libtest_harmonised.py > $T/libtest_harmonised_thr.py
cd $E
ROWS=()
for SP in $SPECIES; do
    for L in $T/benchlibs/$SP/thr_*.fa; do
        [[ -s $L ]] || continue
        M=$(basename $L .fa)
        mkdir -p $SP/$M
        if cmp -s $L $SP/$M/library.fa && [[ -s $SP/$M/$SP.fna.out ]]; then echo "$SP/$M: already masked"; continue; fi
        cp $L $SP/$M/library.fa
        n=$(grep -nP "^$SP\t$M\t" specs.tsv | cut -d: -f1)
        if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$SP" "$M" "$L" >> specs.tsv; n=$(grep -c . specs.tsv); fi
        ROWS+=($(( n - 1 )))
    done
done
echo "[$(date)] RepeatMasker rows (${#ROWS[@]}): ${ROWS[*]}"
# four at a time (each RepeatMasker uses -pa 4)
for r in "${ROWS[@]}"; do
    SLURM_ARRAY_TASK_ID=$r bash $E/rm_one.sh > $E/logs/rm_${r}_thr.out 2>&1 &
    while (( $(jobs -rp | wc -l) >= 4 )); do sleep 20; done
done
wait
echo "[$(date)] harmonised lib-test + family level"
for SP in $SPECIES; do
    VAR=$(ls $T/benchlibs/$SP/thr_*.fa 2>/dev/null | xargs -n1 basename 2>/dev/null | sed 's/\.fa$//' | tr '\n' ' ')
    $PY $T/libtest_harmonised_thr.py $SP $COMPARE $VAR > $T/harm_$SP.log 2>&1 &
done
wait
SLURM_CPUS_PER_TASK=16 bash $E/family_level.sh $SPECIES > $T/famlevel.log 2>&1
echo "[$(date)] summary"
$PY - $SPECIES <<'PY'
import csv, glob, os, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    p = f"{E}/results/{nice}_bench_harmonised_thr.csv"
    if not os.path.exists(p):
        continue
    harm = list(csv.DictReader(open(p)))
    fam = list(csv.DictReader(open(f"{E}/results/{nice}_family_level.csv")))
    methods = sorted({x["Method"] for x in harm}, key=lambda m: (not m.startswith("teannotate"), m))
    print(f"\n{nice}")
    for m in methods:
        h = {x["TE type"]: x for x in harm if x["Method"] == m}
        g = lambda c: f"{float(h[c]['F1']):.3f} (s{float(h[c]['sens']):.2f}/p{float(h[c]['prec']):.2f})" if c in h else "-"
        fa = [x for x in fam if x["Method"] == m and x["class"] == "all"]
        fs = f"fam {fa[0]['perfect']}/{fa[0]['good']}/{fa[0]['present_only']} ({fa[0]['n_library']})" if fa else ""
        print(f"  {m:22s} any {g('any')}  any_dc {g('any_dontcare')}  tir {g('tir')} | {fs}")
PY
echo "[$(date)] all done"

