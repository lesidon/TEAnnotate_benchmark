#!/bin/bash
set -euo pipefail
R=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/maize_mask/parallel_scoring_20261120
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(date -Is) scoring stage $1, task ${SLURM_ARRAY_TASK_ID:-none} on $(hostname)"
if [[ $1 == base || $1 == family ]]; then
    exec "$PY" "$R/parallel_scoring.py" "$1" --task "${SLURM_ARRAY_TASK_ID:?array task required}"
else
    exec "$PY" "$R/parallel_scoring.py" "$1"
fi
