#!/bin/bash
# rescore_wo_cel.sh - 12 Nov 2026: re-run the harmonised lib-test for the without-TEtrimmer ablation arm
# (abl_wo_*, fresh_*) of C. elegans. Its scores in Caenorhabditis_elegans_bench_harmonised_abl.csv were
# overwritten by the later with-arm run (same file name). Uses the existing RepeatMasker outputs; writes
# Caenorhabditis_elegans_bench_harmonised_abl_wo.csv only.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; A=$B/ablation
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
cd $E
SP=caenorhabditis_elegans
M=$(ls -d $E/$SP/fresh_* $E/$SP/abl_wo_* 2>/dev/null | xargs -n1 basename | tr '\n' ' ')
$PY $A/libtest_harmonised_abl_wo.py $SP $M > $A/results/harm_${SP}_wo_rescore.log 2>&1
cat $A/results/harm_${SP}_wo_rescore.log
