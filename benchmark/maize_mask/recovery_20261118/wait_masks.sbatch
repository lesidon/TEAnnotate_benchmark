#!/bin/bash
#SBATCH --job-name=maize_masks_ready
#SBATCH --partition=mem
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=6-00:00:00
#SBATCH --exclude=ct01
set -euo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
exec /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python "$B/maize_mask/recovery_20261118/wait_masks.py"
