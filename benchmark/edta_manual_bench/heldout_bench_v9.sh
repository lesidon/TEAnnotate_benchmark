#!/bin/bash
# heldout_bench_v9.sh <species>... - score the held-out genomes exactly like the
# development genomes: specs rows -> RepeatMasker (reference = curated library, plus every
# method) -> lib-test (default + harmonised) -> family level (Flynn et al. 2020) -> summary.
# Methods: frozen v9 pipeline with / without TEtrimmer, HiTE, EarlGrey, EDTA (v9 module
# script, race-free), HELIANO (plant route only). Missing libraries are skipped.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
L=$B/libraries
cd $E
ROWS=(); DONE=()
row() {  # row <sp> <method> <library>
    local sp=$1 m=$2 lib=$3
    [[ -s $lib ]] || { echo "SKIP $sp/$m: no library ($lib)"; return 1; }
    mkdir -p $sp/$m; cp "$lib" $sp/$m/library.fa
    local n; n=$(grep -nP "^$sp\t$m\t" specs.tsv | cut -d: -f1)
    if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$sp" "$m" "$lib" >> specs.tsv; n=$(grep -c . specs.tsv); fi
    ROWS+=($(( n - 1 )))
    echo "$sp/$m: $(grep -c '>' $sp/$m/library.fa) seqs (specs row $(( n - 1 )))"
}
for SP in "$@"; do
    row $SP reference          $B/curated/$SP/$SP.curated.lib.fa || continue
    ok=0
    row $SP teannotate_ref9    $L/$SP/pipeline_v9/pipeline_library.fa && ok=1
    row $SP teannotate_wo_ref9 $L/$SP/pipeline_wo_tetrimmer_v9/pipeline_library.fa && ok=1
    row $SP hite               $L/$SP/hite/hite_library.fa
    row $SP earlgrey           $L/$SP/earlgrey/earlgrey_library.fa
    if [[ -s $L/$SP/edta2/edta2_library.fa ]]; then row $SP edta2 $L/$SP/edta2/edta2_library.fa
    else row $SP edta2 $L/$SP/edta/edta_library.fa; fi      # plant route: the pipeline's own (v9) EDTA
    row $SP heliano            $L/$SP/heliano/heliano_library.fa
    (( ok )) && DONE+=($SP)
done
(( ${#DONE[@]} )) || { echo "nothing to benchmark"; exit 1; }
echo "[$(date)] RepeatMasker rows: ${ROWS[*]}"
printf '%s\n' "${ROWS[@]}" | xargs -P ${RM_PARALLEL:-8} -I{} bash -c "SLURM_ARRAY_TASK_ID={} bash $E/rm_one.sh > $E/logs/rm_{}_heldout.out 2>&1"
for r in "${ROWS[@]}"; do tail -n 1 $E/logs/rm_${r}_heldout.out; done
echo "[$(date)] lib-test (default)"
bash $E/libtest_species.sh "${DONE[@]}"
echo "[$(date)] lib-test (harmonised)"
for SP in "${DONE[@]}"; do
    $PY $E/libtest_harmonised.py $SP teannotate_ref9 teannotate_wo_ref9 hite earlgrey edta2 heliano &
done
wait
echo "[$(date)] family level"
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16} bash $E/family_level.sh "${DONE[@]}"
echo "[$(date)] summary"
$PY - "${DONE[@]}" <<'PY'
import csv, os, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
M = ["teannotate with tetrimmer (v9)", "teannotate w/o tetrimmer (v9)", "hite", "earlgrey", "edta", "heliano"]
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    print(f"\n{nice}")
    for kind in ("bench_edta_manual", "bench_harmonised"):
        f = f"{E}/results/{nice}_{kind}.csv"
        if not os.path.exists(f): print("  missing", f); continue
        rows = list(csv.DictReader(open(f)))
        for m in M:
            r = {x["TE type"]: x for x in rows if x["Method"] == m}
            if not r: continue
            cols = ["any", "any_dontcare", "ltr", "nonltr", "tir", "helitron"] if kind == "bench_harmonised" else ["any", "ltr", "nonltr", "tir", "helitron"]
            print(f"  {kind[6:]:15s} {m:32s} " + " ".join(f"{c}={float(r[c]['F1']):.3f}" for c in cols if c in r))
    f = f"{E}/results/{nice}_family_level.csv"
    if os.path.exists(f):
        fam = list(csv.DictReader(open(f)))
        print("  family (perfect/good/present; library):",
              "; ".join(f"{x['Method']} {x['perfect']}/{x['good']}/{x['present_only']} ({x['n_library']})"
                        for x in fam if x["class"] == "all"))
PY
echo "[$(date)] all done"

