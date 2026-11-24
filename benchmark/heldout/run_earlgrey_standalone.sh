#!/bin/bash
# run_earlgrey_standalone.sh <species> - EarlGrey baseline (04_run_tools.sh, mirrors the
# pipeline's EarlGrey module) with the pipeline's cached EarlGrey conda env.
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
export EARLGREY_ENV=$B/work/conda_cache/env-e8cd3caef48f718c-8aca2d793016020c698f928e1801d622
export THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "[$(date)] EarlGrey standalone $1 on $(hostname), $THREADS threads"
bash $B/04_run_tools.sh "$1" earlgrey
L=$B/libraries/$1/earlgrey/earlgrey_library.fa
[[ -s $L ]] || { echo "no EarlGrey library for $1" >&2; exit 1; }
echo "[$(date)] EarlGrey library: $(grep -c '>' $L) seqs"
