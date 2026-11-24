#!/bin/bash
#SBATCH --job-name=cbri_ensemble_pair
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=ct[01-10],gn[19,26]
set -euo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
R=$B/edta_manual_bench/ensemble_ablation_20261123
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python "$R/compare.py" score
