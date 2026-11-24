#!/bin/bash
# score_mm_rice.sh [arm...] - score the MMseqs2-dedup rice libraries (combos_mm, method
# mm_<arm>_edta+heliano+hite) exactly like abl_bench.sh: rm_one.sh -> harmonised lib-test
# (+ any_dontcare) -> family level, next to the CD-HIT rows already scored.
ARMS="${*:-wo}"; SP=oryza_sativa; CB=edta+heliano+hite
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; A=$B/ablation
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
cd $E
for ARM in $ARMS; do
    M=mm_${ARM}_$CB; L=$A/combos_mm/$SP/$CB/library_$ARM.fa
    [[ -s $L ]] || { echo "no $L"; continue; }
    mkdir -p $SP/$M
    if ! { cmp -s $L $SP/$M/library.fa && [[ -s $SP/$M/$SP.fna.out ]]; }; then
        cp $L $SP/$M/library.fa
        n=$(grep -nF -- "$SP	$M	" specs.tsv | head -1 | cut -d: -f1)
        if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$SP" "$M" "$L" >> specs.tsv; n=$(grep -c . specs.tsv); fi
        echo "[$(date)] RepeatMasker $M (row $(( n - 1 )))"
        SLURM_ARRAY_TASK_ID=$(( n - 1 )) bash $E/rm_one.sh > $E/logs/rm_$(( n - 1 ))_mm.out 2>&1
    fi
done
echo "[$(date)] harmonised lib-test + family level"
M=$(ls -d $E/$SP/fresh_* $E/$SP/abl_wo_* $E/$SP/abl_with_* $E/$SP/mm_* 2>/dev/null | xargs -n1 basename | tr '\n' ' ')
$PY $A/libtest_harmonised_abl.py $SP $M > $A/results/harm_${SP}_mm.log 2>&1
SLURM_CPUS_PER_TASK=16 bash $E/family_level.sh $SP > $A/results/famlevel_mm.log 2>&1
$PY - <<'PY'
import csv
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench/results"
harm = list(csv.DictReader(open(f"{E}/Oryza_sativa_bench_harmonised_abl.csv")))
fam = list(csv.DictReader(open(f"{E}/Oryza_sativa_family_level.csv")))
for m in sorted({x['Method'] for x in harm}):
    if 'edta+heliano+hite' not in m and not m.startswith('fresh_hite'):
        continue
    h = {x['TE type']: x for x in harm if x['Method'] == m}
    fa = [x for x in fam if x['Method'] == m and x['class'] == 'all']
    g = lambda t, k='F1': float(h[t][k]) if t in h else float('nan')
    print(f"{m:32s} any_dc {g('any_dontcare'):.4f} (sens {g('any_dontcare','sens'):.3f} prec {g('any_dontcare','prec'):.3f})"
          f" ltr {g('ltr'):.3f} nonltr {g('nonltr'):.3f} tir {g('tir'):.3f} hel {g('helitron'):.3f}"
          f" | perfect {fa[0]['perfect'] if fa else '-'} good {fa[0]['good'] if fa else '-'} n_lib {fa[0]['n_library'] if fa else '-'}")
PY
echo "[$(date)] done"
