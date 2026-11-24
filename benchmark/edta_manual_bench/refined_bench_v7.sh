#!/bin/bash
# refined_bench_v5.sh - v5 teannotate arms (pipeline 859baad: evidence gate, TEtrimmer
# re-injection, 80-80-80 CD-HIT, ensemble agreement) into the benchmark, all in one job:
# specs rows -> RepeatMasker -> lib-test (default + harmonised) -> family level -> summary.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
SPECIES="zymoseptoria_tritici septoria_linicola drosophila_melanogaster arabidopsis_thaliana"
cd $E
ROWS=(); DONE=()
for SP in $SPECIES; do
    ok=0
    for arm in teannotate_ref7:pipeline_refined_v7 teannotate_wo_ref7:pipeline_wo_tetrimmer_refined_v7; do
        M=${arm%%:*}; L=$B/libraries/$SP/${arm#*:}/pipeline_library.fa
        [[ -s $L ]] || { echo "SKIP $SP/$M: no library"; continue; }
        ok=1
        mkdir -p $SP/$M; cp "$L" $SP/$M/library.fa
        n=$(grep -nP "^$SP\t$M\t" specs.tsv | cut -d: -f1)
        if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$SP" "$M" "$L" >> specs.tsv; n=$(grep -c . specs.tsv); fi
        ROWS+=($(( n - 1 )))
        echo "$SP/$M: $(grep -c '>' $SP/$M/library.fa) seqs (specs row $(( n - 1 )))"
    done
    (( ok )) && DONE+=($SP)
done
(( ${#DONE[@]} )) || { echo "nothing to benchmark"; exit 1; }
echo "[$(date)] RepeatMasker rows: ${ROWS[*]}"
for r in "${ROWS[@]}"; do SLURM_ARRAY_TASK_ID=$r bash $E/rm_one.sh > $E/logs/rm_${r}_v7.out 2>&1 & done
wait
for r in "${ROWS[@]}"; do tail -n 1 $E/logs/rm_${r}_v7.out; done
echo "[$(date)] lib-test (default)"
bash $E/libtest_species.sh "${DONE[@]}"
echo "[$(date)] lib-test (harmonised)"
for SP in "${DONE[@]}"; do
    $PY $E/libtest_harmonised.py $SP teannotate_ref7 teannotate_wo_ref7 teannotate_ref5 teannotate_wo_ref5 teannotate_ref4 teannotate_wo_ref4 hite earlgrey edta2 heliano &
done
wait
echo "[$(date)] family level"
SLURM_CPUS_PER_TASK=16 bash $E/family_level.sh "${DONE[@]}"
echo "[$(date)] summary"
$PY - "${DONE[@]}" <<'PY'
import csv, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    print(f"\n{nice}")
    for kind in ("bench_edta_manual", "bench_harmonised"):
        rows = list(csv.DictReader(open(f"{E}/results/{nice}_{kind}.csv")))
        for m in ["teannotate with tetrimmer (v7)", "teannotate w/o tetrimmer (v7)", "teannotate with tetrimmer", "teannotate w/o tetrimmer", "teannotate with tetrimmer (v4)",
                  "teannotate w/o tetrimmer (v4)", "hite", "earlgrey", "edta"]:
            r = {x["TE type"]: x for x in rows if x["Method"] == m}
            if not r: continue
            cols = ["any", "any_dontcare", "ltr", "nonltr", "tir", "helitron"] if kind == "bench_harmonised" else ["any", "ltr", "nonltr", "tir", "helitron"]
            print(f"  {kind[6:]:15s} {m:34s} " + " ".join(f"{c}={float(r[c]['F1']):.3f}" for c in cols if c in r))
    fam = list(csv.DictReader(open(f"{E}/results/{nice}_family_level.csv")))
    print("  family (perfect/good/present; library):",
          "; ".join(f"{x['Method']} {x['perfect']}/{x['good']}/{x['present_only']} ({x['n_library']})"
                    for x in fam if x["class"] == "all" and x["Method"] in
                    ("teannotate_ref7", "teannotate_wo_ref7", "teannotate_ref5", "teannotate_wo_ref5", "teannotate_ref4", "teannotate_wo_ref4", "hite", "earlgrey", "edta2")))
PY
echo "[$(date)] all done"
