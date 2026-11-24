#!/bin/bash
# heldout_bench_v10.sh [--register] <species>... - score the held-out genomes with the final
# (v10 library stages) pipeline libraries, exactly like the development genomes:
# RepeatMasker (reference = curated library, plus every method) -> harmonised lib-test
# (+ any_dontcare) -> family level (Flynn et al. 2020) -> summary.
# Methods: teannotate_ref10 / teannotate_wo_ref10 (libraries/<sp>/pipeline{,_wo_tetrimmer}_v10),
# HiTE, Earl Grey, EDTA (standalone edta2; plant route: the pipeline's own EDTA), HELIANO.
# --register only adds the specs.tsv rows and copies the libraries (run once, serially, so
# parallel per-species jobs never append to specs.tsv at the same time). A row whose library
# is unchanged and already masked is not masked again.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
L=$B/libraries
REGISTER=0; [[ $1 == --register ]] && { REGISTER=1; shift; }
cd $E
ROWS=(); DONE=()
row() {  # row <sp> <method> <library>
    local sp=$1 m=$2 lib=$3 n
    [[ -s $lib ]] || { echo "SKIP $sp/$m: no library ($lib)"; return 1; }
    n=$(grep -nP "^$sp\t$m\t" specs.tsv | head -n 1 | cut -d: -f1)
    if (( REGISTER )); then
        mkdir -p $sp/$m
        cmp -s "$lib" $sp/$m/library.fa || { cp "$lib" $sp/$m/library.fa; rm -f $sp/$m/$sp.fna.out; }
        if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$sp" "$m" "$lib" >> specs.tsv; n=$(grep -c . specs.tsv); fi
        echo "$sp/$m: $(grep -c '>' $sp/$m/library.fa) seqs (specs row $(( n - 1 )))"
        return 0
    fi
    [[ -n $n ]] || { echo "ERROR $sp/$m not registered (run --register first)"; return 1; }
    [[ -s $sp/$m/$sp.fna.out ]] || ROWS+=($(( n - 1 )))
}
for SP in "$@"; do
    row $SP reference           $B/curated/$SP/$SP.curated.lib.fa || continue
    ok=0
    row $SP teannotate_ref10    $L/$SP/pipeline_v10/pipeline_library.fa && ok=1
    row $SP teannotate_wo_ref10 $L/$SP/pipeline_wo_tetrimmer_v10/pipeline_library.fa && ok=1
    row $SP hite                $L/$SP/hite/hite_library.fa
    row $SP earlgrey            $L/$SP/earlgrey/earlgrey_library.fa
    if [[ -s $L/$SP/edta2/edta2_library.fa ]]; then row $SP edta2 $L/$SP/edta2/edta2_library.fa
    else row $SP edta2 $L/$SP/edta/edta_library.fa; fi      # plant route: the pipeline's own EDTA
    row $SP heliano             $L/$SP/heliano/heliano_library.fa
    (( ok )) && DONE+=($SP)
done
(( REGISTER )) && exit 0
(( ${#DONE[@]} )) || { echo "nothing to benchmark"; exit 1; }
echo "[$(date)] RepeatMasker rows to mask: ${ROWS[*]:-none}"
(( ${#ROWS[@]} )) && printf '%s\n' "${ROWS[@]}" | xargs -P ${RM_PARALLEL:-6} -I{} bash -c "SLURM_ARRAY_TASK_ID={} bash $E/rm_one.sh > $E/logs/rm_{}_heldout10.out 2>&1"
for r in "${ROWS[@]}"; do tail -n 1 $E/logs/rm_${r}_heldout10.out; done
echo "[$(date)] lib-test (harmonised)"
for SP in "${DONE[@]}"; do
    $PY $E/libtest_harmonised.py $SP teannotate_ref10 teannotate_wo_ref10 hite earlgrey edta2 heliano > $E/logs/harm10_$SP.log 2>&1 &
done
wait
echo "[$(date)] family level"
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16} bash $E/family_level.sh "${DONE[@]}"
echo "[$(date)] summary"
$PY - "${DONE[@]}" <<'PY'
import csv, os, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
M = ["teannotate_ref10", "teannotate_wo_ref10", "hite", "earlgrey", "edta", "edta2", "heliano"]
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    print(f"\n{nice}")
    f = f"{E}/results/{nice}_bench_harmonised.csv"
    rows = list(csv.DictReader(open(f))) if os.path.exists(f) else []
    fam = list(csv.DictReader(open(f"{E}/results/{nice}_family_level.csv"))) if os.path.exists(f"{E}/results/{nice}_family_level.csv") else []
    for m in M:
        r = {x["TE type"]: x for x in rows if x["Method"] == m}
        fa = [x for x in fam if x["Method"] == m and x["class"] == "all"]
        if not r and not fa:
            continue
        g = lambda t: f"{float(r[t]['F1']):.3f}" if t in r else "  -  "
        print(f"  {m:20s} any_dc {g('any_dontcare')} any {g('any')} ltr {g('ltr')} nonltr {g('nonltr')} tir {g('tir')} hel {g('helitron')}"
              + (f" | perfect {fa[0]['perfect']} good {fa[0]['good']} lib {fa[0]['n_library']}" if fa else ""))
PY
echo "[$(date)] all done"
