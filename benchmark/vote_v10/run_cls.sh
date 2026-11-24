#!/bin/bash
# run_cls.sh <plant|animal|fungus> [tag] - all five classifiers (frozen v9 modules) on the
# pooled Unknown-at-vote consensi -> vote_v10/preds/<org>[tag]/ (per-model *.preds.tsv).
set -eo pipefail
ORG=$1; TAG=${2:-}
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v9; V=$B/vote_v10
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
case $ORG in plant) CS=arabidopsis_thaliana;; animal|fungus) CS=drosophila_melanogaster;;
             *) echo "unknown organism_type $ORG"; exit 1;; esac
IN=$V/unknown/$ORG$TAG.fasta
R=$WORK_DIR/vote_v10/$ORG$TAG; mkdir -p $R; cd $R
C=$WORK_DIR/pipeline/$CS
OPTS=(-profile singularity -work-dir $R/work --input $IN --organism_type $ORG --id u_$ORG
      --outdir $V/preds/$ORG$TAG --ensemble_kingdom_aware false --ensemble_method all
      --max_cpus ${SLURM_CPUS_PER_TASK:-16} --max_memory 60GB
      --singularity_cache_dir $WORK_DIR/singularity_heldout --conda_cachedir $WORK_DIR/conda_cache
      -c $B/heldout/heldout_hite16.config)
for c in teclass2 deepte create neuralte; do [[ -d $C/${c}_cache ]] && OPTS+=(--${c}_cache $C/${c}_cache); done
echo "[$(date)] $ORG$TAG: $(grep -c '>' $IN) sequences on $(hostname); pipeline $(cat $H/COMMIT)"
nextflow run $H/run_ensemble.nf "${OPTS[@]}" -resume
find $V/preds/$ORG$TAG -name '*.preds.tsv' | sort
echo "[$(date)] done"

