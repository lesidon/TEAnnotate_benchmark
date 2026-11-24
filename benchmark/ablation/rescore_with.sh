#!/bin/bash
# rescore_with.sh - 12 Nov 2026: re-run the harmonised lib-test for the with-TEtrimmer ablation arm
# (abl_with_<route>, fresh_*, teannotate_ref10/_wo_ref10) of the four first dev genomes. The original
# with-arm scores in <Sp>_bench_harmonised_abl.csv were overwritten by the later wo-arm run (same file
# name). Uses the existing RepeatMasker outputs; writes <Sp>_bench_harmonised_abl_with.csv only.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; A=$B/ablation
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
cd $E
for SP in septoria_linicola zymoseptoria_tritici drosophila_melanogaster arabidopsis_thaliana; do
    M=$(ls -d $E/$SP/fresh_* $E/$SP/abl_with_* 2>/dev/null | xargs -n1 basename | tr '\n' ' ')
    $PY $A/libtest_harmonised_abl_with.py $SP teannotate_ref10 teannotate_wo_ref10 $M > $A/results/harm_${SP}_with_rescore.log 2>&1 &
done
wait
grep -h 'any_dontcare\|rows' $A/results/harm_*_with_rescore.log
