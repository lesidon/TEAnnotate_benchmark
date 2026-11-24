#!/bin/bash
#SBATCH --job-name=ho_stub2
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/logs/ho_stub2_%j.out
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
for SP in caenorhabditis_elegans blumeria_hordei; do
    echo "===== $SP"
    bash $B/heldout/run_heldout_v9.sh $SP stub || echo "STUB FAILED: $SP"
done
