#!/bin/bash
# refined_bench.sh <species>... - add the refined-dedup teannotate arms (with / w/o
# TEtrimmer) to specs.tsv, RepeatMasker them, then lib-test + CSV (all methods).
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
cd $E
FIRST=$(grep -c . specs.tsv)
DONE=()
for SP in "$@"; do
    ok=0
    for arm in teannotate_ref:pipeline_refined teannotate_wo_ref:pipeline_wo_tetrimmer_refined; do
        M=${arm%%:*}; L=$B/libraries/$SP/${arm#*:}/pipeline_library.fa
        [[ -s $L ]] || { echo "SKIP $SP/$M: no library"; continue; }
        ok=1
        grep -qP "^$SP\t$M\t" specs.tsv && { echo "$SP/$M already in specs.tsv"; continue; }
        mkdir -p $SP/$M; cp "$L" $SP/$M/library.fa
        printf '%s\t%s\t%s\n' "$SP" "$M" "$L" >> specs.tsv
        echo "$SP/$M: $(grep -c '>' $SP/$M/library.fa) seqs"
    done
    (( ok )) && DONE+=($SP)
done
LAST=$(( $(grep -c . specs.tsv) - 1 ))
(( ${#DONE[@]} )) || { echo "nothing to benchmark"; exit 1; }
DEP=""
if (( LAST >= FIRST )); then
    A=$(sbatch --parsable -p htc --array=$FIRST-$LAST -c 16 --mem=32G -t 12:00:00 --job-name=edta_rm_ref \
          -o "$E/logs/rm_%a.out" -e "$E/logs/rm_%a.err" "$E/rm_one.sh")
    DEP="--dependency=afterany:$A"; echo "RM array $A (rows $FIRST-$LAST)"
fi
# also wait for other benchmark RepeatMasker / lib-test jobs still in the queue
OTHER=$(squeue -h -u leonid.sidorov -n edta_rm,edta_rm_fungi,edta_rm_wo,edta_rm_ref,edta_libtest_fungi,edta_libtest_wo -o %A | sort -u | paste -sd: -)
if [[ -n $OTHER ]]; then DEP="${DEP:---dependency=afterany}:$OTHER"; fi
L=$(sbatch --parsable -p htc -c 4 --mem=32G -t 06:00:00 --job-name=edta_libtest_ref $DEP \
      -o "$E/logs/libtest_ref_%j.out" -e "$E/logs/libtest_ref_%j.err" "$E/libtest_species.sh" "${DONE[@]}")
echo "lib-test $L for ${DONE[*]}"
