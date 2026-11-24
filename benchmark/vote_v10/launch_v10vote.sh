#!/bin/bash
# Replay the v10 two-step vote on the frozen v9 runs of the 4 development genomes + C. briggsae:
# pooled Unknown-at-vote consensi -> all five classifiers (per organism_type) -> new vote ->
# lib-test of the relabelled annotation.
set -e
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
V=$B/vote_v10; PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
SPS="arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola caenorhabditis_briggsae"
cd $V; mkdir -p logs
chmod +x run_cls.sh
for f in vote_v10_common.py prep_unknown.py apply_vote.py bench_eval.py; do $PY -m py_compile $f; done
echo "syntax ok"
$PY prep_unknown.py $SPS
DEPS=""
for org in plant animal fungus; do
    j=$(sbatch --parsable -p htc -c 16 --mem=60G -t 12:00:00 -J v10cls_$org -o $V/logs/cls_${org}_%j.out $V/run_cls.sh $org)
    DEPS="$DEPS:$j"; echo "$org classification: $j"
done
j=$(sbatch --parsable -p htc -c 16 --mem=48G -t 06:00:00 --dependency=afterok$DEPS -J v10bench \
    -o $V/logs/apply_bench_%j.out --wrap "cd $V && $PY apply_vote.py $SPS && $PY bench_eval.py $SPS")
echo "apply + bench: $j (after$DEPS)"
squeue -u leonid.sidorov -o '%i %P %j %T %R' | grep -E 'v10'

