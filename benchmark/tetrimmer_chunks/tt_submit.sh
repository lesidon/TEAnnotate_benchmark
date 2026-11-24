#!/bin/bash
# tt_submit.sh <task_dir> <chunk_root> <job name> [sbatch options for one chunk...] - run every
# chunk of <chunk_root> that is not complete yet, each as its own Slurm job <job name>_cXX, fed by
# tt_feed.sh: at most TT_PARALLEL chunk jobs queued at a time, inside the per-user submit limit of
# the partition's QOS (htc allows 16 queued jobs and counts every array task, so one array over
# all chunks can be refused), and a chunk that stopped early (1-day htc limit) is submitted again
# and continues with --continue_analysis. One feeding round runs here, so chunks queue at once;
# then a feeder job (mem partition, 1 CPU, up to 6 days) keeps feeding. Prints only the feeder's
# job id: it ends OK only when every chunk is complete, so "sbatch -d afterok:<id>" waits for all.
# Chunk job default: htc, 36 CPUs, 150 G, 1 day. TT_PARALLEL (default 5), TT_THREADS (default
# 1.5 x CPUs, tt_run_chunk.sh), TT_MAX_TRIES (5) and TT_SUBMIT_SPARE (2) reach the feeder.
set -eo pipefail
X=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/tetrimmer_chunks
T=$1; R=$2; J=$3; shift 3
OPTS=("$@"); (( ${#OPTS[@]} )) || OPTS=(-p htc -c 36 --mem=150G -t 1-00:00:00)
for v in $(compgen -e | grep -E '^SLURM_'); do unset $v; done     # submit as from a login shell
mkdir -p $R/logs $R/feed
ENV=(TT_PARALLEL=${TT_PARALLEL:-5} TT_MAX_TRIES=${TT_MAX_TRIES:-5} TT_SUBMIT_SPARE=${TT_SUBMIT_SPARE:-2})
[[ -n $TT_THREADS ]] && ENV+=(TT_THREADS=$TT_THREADS)
env "${ENV[@]}" TT_FEED_ONCE=1 bash $X/tt_feed.sh $T $R $J "${OPTS[@]}" >> $R/logs/${J}_feed_first.out 2>&1
SB=(sbatch); [[ -n $TT_FEED_DRY ]] && SB=(echo sbatch)
"${SB[@]}" --parsable -p ${TT_FEED_PARTITION:-mem} -c 1 --mem=1G -t ${TT_FEED_TIME:-6-00:00:00} \
    -J ${J}_feed -o $R/logs/${J}_feed_%j.out \
    --wrap "${ENV[*]} bash $X/tt_feed.sh $T $R $J ${OPTS[*]}"
