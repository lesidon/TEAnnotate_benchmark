#!/bin/bash
# run_refined_v10.sh <species> - final development rerun: the library stages
# (run_reclassify.nf --dedup: CD-HIT + DEDUP_REFINE, evidence gate, curation,
# classification) of pipeline_head_v10 on the species' raw de novo library (the same input
# as v9), twice from one launch dir (the second run reuses the dedup and the gate):
#   with TEtrimmer  -> libraries/<sp>/pipeline_refined_v10/pipeline_library.fa
#   w/o  TEtrimmer  -> libraries/<sp>/pipeline_wo_tetrimmer_refined_v10/pipeline_library.fa
# Assets are passed explicitly (Drosophila's caches: the pinned commits, Pfam 38.2), so
# nothing is downloaded and every species uses the same ones.
set -eo pipefail
SP="$1"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v10
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$GENOME_DIR/$KING/$SP/$SP.fna
R=$WORK_DIR/pipeline/$SP
RAW=$R/out/$SP/02_denovo_te/$SP.raw_library.fa; [[ -s $WORK_DIR/edta_fixed/$SP/raw_library.fa ]] && RAW=$WORK_DIR/edta_fixed/$SP/raw_library.fa
[[ -s $RAW && -s $G ]] || { echo "missing raw library or genome for $SP" >&2; exit 1; }
C=$WORK_DIR/pipeline/drosophila_melanogaster
LD=$WORK_DIR/refined_dedup_v10/$SP
mkdir -p "$LD"; cd "$LD"
OPTS=(-profile singularity -work-dir "$LD/work"
      --input "$RAW" --dedup --genome "$G" --id "$SP" --organism_type "$ORG"
      --max_cpus "${SLURM_CPUS_PER_TASK:-16}"
      --max_memory "$(( ${SLURM_MEM_PER_NODE:-65536} * 9 / 10 / 1024 ))GB"
      --conda_cachedir "$WORK_DIR/conda_cache" --singularity_cache_dir "$WORK_DIR/singularity_heldout"
      --pfam_db "$C/pfam_cache/pfam_db"
      --teclass2_repo "$C/teclass2_cache/repo/TEclass2" --teclass2_model "$C/teclass2_cache/model/teclass2_model"
      --create_repo "$C/create_cache/repo/CREATE" --create_model "$C/create_cache/model/create_model"
      --neuralte_repo "$C/neuralte_cache/repo/NeuralTE")
echo "[$(date)] $SP ($ORG) on $(hostname): $(grep -c '^>' $RAW) raw seqs ($RAW); pipeline $(cat $H/COMMIT)"
nextflow run "$H/run_reclassify.nf" "${OPTS[@]}" --outdir "$LD/out_with" -resume
nextflow run "$H/run_reclassify.nf" "${OPTS[@]}" --skip_tetrimmer --outdir "$LD/out_wo" -resume
for arm in with:pipeline_refined_v10 wo:pipeline_wo_tetrimmer_refined_v10; do
    LIB="$LD/out_${arm%%:*}/$SP/05_classification/$SP.te_library.fa"
    [[ -s $LIB ]] || { echo "no final library in out_${arm%%:*}" >&2; exit 1; }
    mkdir -p "$LIB_DIR/$SP/${arm#*:}"
    cp "$LIB" "$LIB_DIR/$SP/${arm#*:}/pipeline_library.fa"
    echo "[$(date)] ${arm%%:*} TEtrimmer: $(grep -c '^>' "$LIB") seqs -> libraries/$SP/${arm#*:}/"
done
for t in $(find "$LD/out_with" "$LD/out_wo" -name "*.dedup_refine.tsv" -o -name "*.tetrimmer_postprocess.tsv" | sort -u); do
    echo "== ${t#$LD/}"; awk -F"\t" "NR>1{print \$2}" $t | sort | uniq -c
done
for e in $(find "$LD/out_with" -name "*.evidence.tsv"); do
    echo "== ${e#$LD/} (tier action)"; awk -F"\t" "NR>1{print \$9, \$10}" $e | sort | uniq -c
done

