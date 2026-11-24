#!/bin/bash
#SBATCH --job-name=ho_stub
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/logs/ho_stub_%j.out
# Wiring test of the frozen v9 main.nf on one genome per route (nextflow -stub).
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
for SP in caenorhabditis_elegans blumeria_hordei chlamydomonas_reinhardtii; do
    echo "===== $SP"
    bash $B/heldout/run_heldout_v9.sh $SP stub || echo "STUB FAILED: $SP"
done

