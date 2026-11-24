#!/bin/bash
#SBATCH -o logs/score_%x_%j.out
set -euo pipefail
B=$(readlink -f /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark)
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
echo "$(date -Is) maize scoring stage $1 on $(hostname)"
"$PY" "$B/maize_mask/maize_scoring_20261117.py" "$1"
echo "$(date -Is) completed maize scoring stage $1"
