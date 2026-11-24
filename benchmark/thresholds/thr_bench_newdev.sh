#!/bin/bash
#SBATCH --job-name=thr_newdev
#SBATCH --partition=mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --output=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/thresholds/thr_newdev_%j.out
# thr_bench_newdev.sh - the copy-threshold grid of thr_bench.sh (Table S3) for the two later
# development genomes, C. elegans and rice. The gate variants are rebuilt by the unchanged
# thr_variants.py from the v10 w/o-TEtrimmer runs of the route combination in the ablation
# (ablation/reclass/<sp>/<route>/out_wo; their final libraries are the scored abl_wo_<route>),
# masked with rm_one.sh, scored with the unchanged libtest_harmonised_thr.py against the
# default (abl_wo_<route> = F2/H5), and evaluated at family level. family_level.py rewrites
# results/<Sp>_family_level.csv, so it is backed up first and every pre-existing row is checked.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; T=$B/thresholds
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
SPECIES="caenorhabditis_elegans oryza_sativa"
declare -A ROUTE=([caenorhabditis_elegans]=earlgrey+hite [oryza_sativa]=edta+heliano+hite)
nice() { echo "$(tr '[:lower:]' '[:upper:]' <<< ${1:0:1})${1:1}"; }
echo "[$(date)] start on $(hostname)"

# run outputs in the layout thr_variants.py expects: <work>/<sp>/out_wo/<sp>
S=$T/newdev_runs; mkdir -p $S
for SP in $SPECIES; do ln -sfn $B/ablation/reclass/$SP/${ROUTE[$SP]} $S/$SP; done
$PY $T/thr_variants.py $S $T/benchlibs $SPECIES || { echo "thr_variants.py failed"; exit 1; }

mkdir -p $T/backup_famlevel
for SP in $SPECIES; do cp -p $E/results/$(nice $SP)_family_level.csv $T/backup_famlevel/; done

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
# eight at a time (each RepeatMasker uses -pa 4)
for r in "${ROWS[@]}"; do
    SLURM_ARRAY_TASK_ID=$r bash $E/rm_one.sh > $E/logs/rm_${r}_thr.out 2>&1 &
    while (( $(jobs -rp | wc -l) >= 8 )); do sleep 20; done
done
wait
for r in "${ROWS[@]}"; do grep -q FAILED $E/logs/rm_${r}_thr.out && echo "RepeatMasker FAILED: row $r"; done

echo "[$(date)] harmonised lib-test + family level"
for SP in $SPECIES; do
    VAR=$(ls $T/benchlibs/$SP/thr_*.fa 2>/dev/null | xargs -n1 basename 2>/dev/null | sed 's/\.fa$//' | tr '\n' ' ')
    $PY $T/libtest_harmonised_thr.py $SP abl_wo_${ROUTE[$SP]} $VAR > $T/harm_$SP.log 2>&1 &
done
wait
SLURM_CPUS_PER_TASK=32 bash $E/family_level.sh $SPECIES > $T/famlevel_newdev.log 2>&1

echo "[$(date)] pre-existing family-level rows"
$PY - $SPECIES <<'PY'
import csv, shutil, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
T = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/thresholds"
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    old = list(csv.DictReader(open(f"{T}/backup_famlevel/{nice}_family_level.csv")))
    cur_path = f"{E}/results/{nice}_family_level.csv"
    new = list(csv.DictReader(open(cur_path)))
    key = lambda r: (r["Method"], r["class"])
    nd = {key(r): r for r in new}
    changed = [key(r) for r in old if nd.get(key(r)) != r]
    added = sorted({r["Method"] for r in new} - {r["Method"] for r in old})
    print(f"{nice}: {len(old)} old rows, {len(changed)} changed or missing, new methods {added}")
    if changed:
        # keep the published rows exactly; add only the new thr_* rows
        shutil.copy(cur_path, cur_path.replace(".csv", ".regenerated.csv"))
        with open(cur_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(old[0].keys()), lineterminator="\n")
            w.writeheader()
            w.writerows(old)
            w.writerows(r for r in new if r["Method"].startswith("thr_") and r["Method"] not in {o["Method"] for o in old})
        print(f"  restored the {len(old)} old rows and appended the thr_* rows; regenerated file kept as *.regenerated.csv")
        for k in changed[:10]:
            print("  changed:", k)
PY

echo "[$(date)] summary"
$PY - $SPECIES <<'PY'
import csv, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
for sp in sys.argv[1:]:
    nice = sp[0].upper() + sp[1:]
    harm = list(csv.DictReader(open(f"{E}/results/{nice}_bench_harmonised_thr.csv")))
    fam = list(csv.DictReader(open(f"{E}/results/{nice}_family_level.csv")))
    print(f"\n{nice}")
    for m in sorted({x["Method"] for x in harm}):
        h = {x["TE type"]: x for x in harm if x["Method"] == m}
        fa = [x for x in fam if x["Method"] == m and x["class"] == "all"]
        dc = h.get("any_dontcare")
        print(f"  {m:28s} any_dc F1 {float(dc['F1']):.3f} (s{float(dc['sens']):.2f}/p{float(dc['prec']):.2f})"
              + (f" | perfect {fa[0]['perfect']} good {fa[0]['good']} lib {fa[0]['n_library']}" if fa else ""))
PY
echo "[$(date)] all done"
