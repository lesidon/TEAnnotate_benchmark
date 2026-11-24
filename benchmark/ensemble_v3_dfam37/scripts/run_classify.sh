#!/bin/bash
# run_classify.sh <plant|animal|fungi> - all five classifiers (frozen v9 modules, v9 defaults
# incl. the shared 0.7 confidence threshold) on one kingdom pool -> preds/<k>/voting_table.tsv
set -eo pipefail
K=$1
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v9; X=$B/ensemble_v2
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
case $K in plant) ORG=plant; CS=arabidopsis_thaliana;; animal) ORG=animal; CS=drosophila_melanogaster;;
           fungi) ORG=fungus; CS=drosophila_melanogaster;; *) echo "unknown kingdom $K"; exit 1;; esac
R=$WORK_DIR/ensemble_v2/$K; mkdir -p $R; cd $R
C=$WORK_DIR/pipeline/$CS
OPTS=(-profile singularity -work-dir $R/work --input $X/pools/$K.fasta --organism_type $ORG --id ens_$K
      --outdir $X/preds/$K --ensemble_kingdom_aware false --ensemble_method all
      --max_cpus ${SLURM_CPUS_PER_TASK:-32} --max_memory 144GB
      --singularity_cache_dir $WORK_DIR/singularity_heldout --conda_cachedir $WORK_DIR/conda_cache
      -c $B/heldout/heldout_hite16.config)
for c in teclass2 deepte create neuralte; do [[ -d $C/${c}_cache ]] && OPTS+=(--${c}_cache $C/${c}_cache); done
echo "[$(date)] $K ($ORG): $(grep -c '>' $X/pools/$K.fasta) sequences on $(hostname); pipeline $(cat $H/COMMIT)"
nextflow run $H/run_ensemble.nf "${OPTS[@]}" -resume
echo "[$(date)] done: $(($(wc -l < $X/preds/$K/voting_table.tsv) - 1)) rows in preds/$K/voting_table.tsv"

