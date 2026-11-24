#!/bin/bash
# run_heldout_v12.sh <species> [run root] [raw library] [provenance]
# As run_heldout_v10.sh, but with pipeline_head_v12 (= v10 library stages + MMseqs2
# clustering for libraries of >= 20,000 sequences, which only the maize raw library reaches).
# Own launch dir (heldout_v12/<sp>): never -resume from another checkout path.
#   with -> libraries/<sp>/pipeline_v12/pipeline_library.fa
#   wo   -> libraries/<sp>/pipeline_wo_tetrimmer_v12/pipeline_library.fa
set -eo pipefail
SP="$1"; ROOT="${2:-heldout_v9}"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v12
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
D=$WORK_DIR/$ROOT/$SP/out_with/$SP/02_denovo_te
RAW="${3:-$D/$SP.raw_library.fa}"; PROV="${4:-$D/$SP.provenance.tsv}"
[[ -s $RAW && -s $PROV && -s $G ]] || { echo "missing raw library, provenance or genome for $SP" >&2; exit 1; }
C=$WORK_DIR/pipeline/drosophila_melanogaster
LD=$WORK_DIR/heldout_v12/$SP
mkdir -p $LD; cd $LD
OPTS=(-profile singularity -work-dir $LD/work
      --input $RAW --provenance $PROV --dedup --genome $G --id $SP --organism_type $ORG
      --max_cpus ${SLURM_CPUS_PER_TASK:-16}
      --max_memory $(( ${SLURM_MEM_PER_NODE:-65536} * 9 / 10 / 1024 ))GB
      --conda_cachedir $WORK_DIR/conda_cache --singularity_cache_dir $WORK_DIR/singularity_heldout
      --pfam_db $C/pfam_cache/pfam_db
      --teclass2_repo $C/teclass2_cache/repo/TEclass2 --teclass2_model $C/teclass2_cache/model/teclass2_model
      --create_repo $C/create_cache/repo/CREATE --create_model $C/create_cache/model/create_model
      --neuralte_repo $C/neuralte_cache/repo/NeuralTE)
echo "[$(date)] $SP ($ORG) on $(hostname): $(grep -c '^>' $RAW) raw seqs ($RAW); pipeline $(head -n 1 $H/COMMIT)"
nextflow run $H/run_reclassify.nf "${OPTS[@]}" --outdir $LD/out_with -resume
nextflow run $H/run_reclassify.nf "${OPTS[@]}" --skip_tetrimmer --outdir $LD/out_wo -resume
for arm in with:pipeline_v12 wo:pipeline_wo_tetrimmer_v12; do
    LIB=$LD/out_${arm%%:*}/$SP/05_classification/$SP.te_library.fa
    [[ -s $LIB ]] || { echo "no final library in out_${arm%%:*}" >&2; exit 1; }
    mkdir -p $LIB_DIR/$SP/${arm#*:}
    cp $LIB $LIB_DIR/$SP/${arm#*:}/pipeline_library.fa
    echo "[$(date)] ${arm%%:*} TEtrimmer: $(grep -c '^>' $LIB) seqs -> libraries/$SP/${arm#*:}/"
done
grep -h 'MMseqs2 instead' $LD/work/*/*/.command.err 2>/dev/null || true
echo "[$(date)] done"
