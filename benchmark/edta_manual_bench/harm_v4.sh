#!/bin/bash
E=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
for sp in arabidopsis_thaliana drosophila_melanogaster; do
  $PY $E/libtest_harmonised.py $sp teannotate_wo_ref4 hite earlgrey edta2 > $E/logs/harm_v4_$sp.log 2>&1 &
done
wait
# the CSV written above is partial (four methods); the v5 benchmark job rewrites it in full
