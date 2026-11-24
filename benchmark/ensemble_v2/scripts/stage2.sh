#!/bin/bash
#SBATCH --job-name=ens_pools
#SBATCH --partition=htc
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v2/logs/stage2_%j.out
# pools (PanTEon held-out + MycoMobilome held-out) -> MMseqs2 clusters -> cluster-aware split
set -eo pipefail
X=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v2
P=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/ensemble_prediction
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
MMSEQS=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/mmseqs
cd $X
$PY scripts/make_pools.py --panteon-leakage-dir $P/leakage_out --panteon-tsv-dir $P/panteon_partitions \
    --mm-fasta data/mycomobilome_pe.fasta --mm-tsv data/mycomobilome_pe.tsv \
    --mm-heldout-ids leakage/mycomobilome_pe.heldout.ids --label-map scripts/label_map.json --outdir pools
for k in plant animal fungi; do
    rm -rf split/tmp_$k
    $MMSEQS easy-cluster pools/$k.fasta split/$k split/tmp_$k --min-seq-id 0.5 -c 0.5 --cov-mode 0 \
        --threads ${SLURM_CPUS_PER_TASK:-32} > logs/mmseqs_cluster_$k.log 2>&1
    $PY scripts/split.py --truth pools/$k.truth.tsv --clusters split/${k}_cluster.tsv --out split/$k
    rm -rf split/tmp_$k split/${k}_all_seqs.fasta
done
echo "[$(date)] stage 2 done"

