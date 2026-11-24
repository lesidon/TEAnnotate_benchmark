#!/bin/bash
#SBATCH --job-name=ens_leak
#SBATCH --partition=mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v2/logs/leakage_%j.out
# Leakage control for the MycoMobilome candidates: same training corpora and thresholds as
# the PanTEon held-out sets (RepBase 29.10 + Dfam curated 3.8; >=70% identity over >=50%
# coverage -> leaked).
set -eo pipefail
X=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v2
MMSEQS=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/mmseqs
P=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/ensemble_prediction
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
cd $X
mkdir -p leakage/tmp
echo "[$(date)] start on $(hostname)"
$PY scripts/leakage_filter.py --query data/mycomobilome_pe.fasta \
    --training repbase=$P/training_fasta/repbase.fasta dfam=$P/training_fasta/dfam.fasta \
    --identity 0.7 --coverage 0.5 --threads ${SLURM_CPUS_PER_TASK:-32} --mmseqs $MMSEQS \
    --tmpdir $X/leakage/tmp --outdir leakage --prefix mycomobilome_pe --write-leaked
cat leakage/mycomobilome_pe.leakage_summary.tsv
echo "[$(date)] done"

