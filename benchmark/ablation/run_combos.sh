#!/bin/bash
# run_combos.sh <species> [arm] - the pipeline's library stages (run_reclassify.nf --dedup of
# pipeline_head_v10: CD-HIT + DEDUP_REFINE, evidence gate with the combination's provenance,
# curation, classification) on every tool combination of a species. arm = wo (default,
# --skip_tetrimmer, the screen) or with (TEtrimmer; for the shortlisted combinations given
# in $COMBOS). Final library -> combos/<sp>/<combo>/library_<arm>.fa, trace -> trace_<arm>.tsv.
set -eo pipefail
SP="$1"; ARM="${2:-wo}"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
A=$B/ablation; H=$B/pipeline_head_v10
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$GENOME_DIR/$KING/$SP/$SP.fna
C=$WORK_DIR/pipeline/drosophila_melanogaster
[[ -d $A/combos/$SP ]] || $PY $A/make_combos.py $A $SP
LIST=${COMBOS:-$(ls $A/combos/$SP)}
echo "[$(date)] $SP ($ORG) arm=$ARM on $(hostname); combos: $(echo $LIST | wc -w); pipeline $(cat $H/COMMIT)"
for CB in $LIST; do
    D=$A/combos/$SP/$CB
    [[ -s $D/library_$ARM.fa ]] && { echo "$CB: done"; continue; }
    LD=$A/reclass/$SP/$CB; mkdir -p $LD; cd $LD
    OPTS=(-profile singularity -work-dir $LD/work
          --input $D/raw.fa --provenance $D/provenance.tsv --dedup --genome $G --id $SP --organism_type $ORG
          --max_cpus ${SLURM_CPUS_PER_TASK:-16}
          --max_memory "$(( ${SLURM_MEM_PER_NODE:-32768} * 9 / 10 / 1024 ))GB"
          --conda_cachedir $WORK_DIR/conda_cache --singularity_cache_dir $WORK_DIR/singularity_heldout
          --pfam_db $C/pfam_cache/pfam_db
          --teclass2_repo $C/teclass2_cache/repo/TEclass2 --teclass2_model $C/teclass2_cache/model/teclass2_model
          --create_repo $C/create_cache/repo/CREATE --create_model $C/create_cache/model/create_model
          --neuralte_repo $C/neuralte_cache/repo/NeuralTE
          -with-trace $LD/trace_$ARM.tsv)
    [[ $ARM == wo ]] && OPTS+=(--skip_tetrimmer)
    if nextflow run $H/run_reclassify.nf "${OPTS[@]}" -c $A/trace_raw.config --outdir $LD/out_$ARM -resume > $LD/nf_$ARM.log 2>&1; then
        cp $LD/out_$ARM/$SP/05_classification/$SP.te_library.fa $D/library_$ARM.fa
        cp $LD/trace_$ARM.tsv $D/trace_$ARM.tsv
        echo "[$(date)] $CB: $(grep -c '^>' $D/raw.fa) raw -> $(grep -c '^>' $D/library_$ARM.fa) final"
    else
        echo "[$(date)] $CB: FAILED (see $LD/nf_$ARM.log)"; tail -5 $LD/nf_$ARM.log
    fi
done
echo "[$(date)] $SP done"
