#!/bin/bash
# v10-vote replay on the held-out genomes (evaluation only): classify their Unknown-at-vote
# consensi now; apply the vote + score after the held-out bench (8484269) wrote the v9 .out.
set -e
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
V=$B/vote_v10; PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
SPS="caenorhabditis_elegans blumeria_hordei ectocarpus_sp7"
cd $V
for f in prep_unknown.py apply_vote.py; do $PY -m py_compile $f; done
TAG=_ho $PY prep_unknown.py $SPS
DEPS=""
for org in animal fungus; do
    [[ -s unknown/${org}_ho.fasta ]] || continue
    j=$(sbatch --parsable -p htc -c 16 --mem=60G -t 12:00:00 -J v10cls_${org}_ho -o $V/logs/cls_${org}_ho_%j.out $V/run_cls.sh $org _ho)
    DEPS="$DEPS:$j"; echo "$org classification: $j"
done
j=$(sbatch --parsable -p htc -c 16 --mem=48G -t 06:00:00 --dependency=afterok$DEPS:8484269 -J v10bench_ho \
    -o $V/logs/apply_bench_ho_%j.out \
    --wrap "cd $V && TAG=_ho $PY apply_vote.py $SPS && TAG=_ho $PY bench_eval.py $SPS && $PY diag_relabel.py $SPS && $PY agree_check.py $SPS")
echo "apply + bench + diagnostics: $j (after$DEPS:8484269)"

