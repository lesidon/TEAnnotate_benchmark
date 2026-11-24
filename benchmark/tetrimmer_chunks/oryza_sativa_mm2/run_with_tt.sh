#!/bin/bash
# run_with_tt.sh (copy of ablation/run_mm2_rice.sh + the merged chunked TEtrimmer output) [arm] - validation of the both-strand MMseqs2 search + greedy dedup path (pipeline_head_v13) on rice:
# the library stages of the plant route combination (edta+heliano+hite, the same raw library
# and provenance as the ablation) with --dedup_mmseqs_min_seqs 1, so every CD-HIT step
# clusters with MMseqs2 instead. Compared with abl_<arm>_edta+heliano+hite (CD-HIT, v10).
# Final library -> ablation/combos_mm2/oryza_sativa/edta+heliano+hite/library_<arm>.fa
set -eo pipefail
ARM="${1:-wo}"; SP=oryza_sativa; CB=edta+heliano+hite
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
A=$B/ablation; H=$B/pipeline_head_v13
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$GENOME_DIR/$KING/$SP/$SP.fna
C=$WORK_DIR/pipeline/drosophila_melanogaster
D=$A/combos/$SP/$CB
OUT=$A/combos_mm2/$SP/$CB; mkdir -p $OUT
LD=$A/reclass_mm2/$SP/$CB; mkdir -p $LD; cd $LD
echo "[$(date)] $SP $CB arm=$ARM on $(hostname); pipeline $(cat $H/COMMIT); MMseqs2 forced"
OPTS=(-profile singularity -work-dir $LD/work
      --input $D/raw.fa --provenance $D/provenance.tsv --dedup --genome $G --id $SP --organism_type $ORG
      --dedup_mmseqs_min_seqs 1
      --max_cpus ${SLURM_CPUS_PER_TASK:-16}
      --max_memory "$(( ${SLURM_MEM_PER_NODE:-32768} * 9 / 10 / 1024 ))GB"
      --conda_cachedir $WORK_DIR/conda_cache --singularity_cache_dir $WORK_DIR/singularity_heldout
      --pfam_db $C/pfam_cache/pfam_db
      --teclass2_repo $C/teclass2_cache/repo/TEclass2 --teclass2_model $C/teclass2_cache/model/teclass2_model
      --create_repo $C/create_cache/repo/CREATE --create_model $C/create_cache/model/create_model
      --neuralte_repo $C/neuralte_cache/repo/NeuralTE
      -with-trace $LD/trace_$ARM.tsv)
[[ $ARM == wo ]] && OPTS+=(--skip_tetrimmer)
nextflow run $H/run_reclassify.nf "${OPTS[@]}" -c $A/trace_raw.config -c /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/tetrimmer_chunks/oryza_sativa_mm2/precomputed.config --outdir $LD/out_$ARM -resume
cp $LD/out_$ARM/$SP/05_classification/$SP.te_library.fa $OUT/library_$ARM.fa
cp $LD/trace_$ARM.tsv $OUT/trace_$ARM.tsv
echo "[$(date)] $(grep -c '^>' $D/raw.fa) raw -> $(grep -c '^>' $OUT/library_$ARM.fa) final"
grep -h 'MMseqs2 instead' $LD/work/*/*/.command.err 2>/dev/null || true
