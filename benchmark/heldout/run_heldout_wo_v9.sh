#!/bin/bash
# run_heldout_wo_v9.sh <species> - finish a held-out genome whose with-TEtrimmer arm of
# run_heldout_v9b*.sh is done: (1) the copies that runner makes after both arms (with-arm
# library and the route's tool libraries), (2) the w/o-TEtrimmer arm as the library stages
# only (pipeline_head_v9/run_reclassify.nf --dedup --skip_tetrimmer) on the with-arm's raw
# library and provenance, so both arms share one raw library. Needed because the second
# main.nf run re-ran Earl Grey instead of resuming it (zebrafish, 2026-10-21): a new,
# stochastic raw library would make the arms incomparable, and cost days.
#   with -> libraries/<sp>/pipeline_v9/pipeline_library.fa
#   wo   -> libraries/<sp>/pipeline_wo_tetrimmer_v9/pipeline_library.fa
set -eo pipefail
SP="$1"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v9
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
R=$WORK_DIR/heldout_v9/$SP
D=$R/out_with/$SP/02_denovo_te
RAW=$D/$SP.raw_library.fa; PROV=$D/$SP.provenance.tsv
WITH=$R/out_with/$SP/05_classification/$SP.classified_library.fa
[[ -s $RAW && -s $PROV && -s $WITH && -s $G ]] || { echo "missing with-arm outputs or genome for $SP" >&2; exit 1; }

# (1) with-arm copies, as run_heldout_v9b*.sh does after both arms
mkdir -p $LIB_DIR/$SP/pipeline_v9
cp $WITH $LIB_DIR/$SP/pipeline_v9/pipeline_library.fa
echo "[$(date)] with TEtrimmer: $(grep -c '^>' $WITH) seqs -> libraries/$SP/pipeline_v9/"
for t in hite earlgrey edta heliano; do
    L=$(find $D -name "$SP.${t}_library.fa" -print -quit 2>/dev/null || true)
    if [[ -n $L && -s $L ]]; then
        mkdir -p $LIB_DIR/$SP/$t; cp $L $LIB_DIR/$SP/$t/${t}_library.fa
        echo "[$(date)] $t: $(grep -c '^>' $L) seqs -> libraries/$SP/$t/"
    fi
done

# (2) w/o-TEtrimmer arm: library stages on the same raw library
C=$WORK_DIR/pipeline/drosophila_melanogaster
LD=$WORK_DIR/heldout_v9_wo/$SP
mkdir -p $LD; cd $LD
OPTS=(-profile singularity -work-dir $LD/work
      --input $RAW --dedup --provenance $PROV --genome $G --id $SP --organism_type $ORG
      --skip_tetrimmer --outdir $LD/out_wo
      --max_cpus ${SLURM_CPUS_PER_TASK:-16}
      --max_memory $(( ${SLURM_MEM_PER_NODE:-65536} * 9 / 10 / 1024 ))GB
      --conda_cachedir $WORK_DIR/conda_cache --singularity_cache_dir $WORK_DIR/singularity_heldout)
for c in pfam teclass2 deepte create neuralte; do
    [[ -d $C/${c}_cache ]] && OPTS+=(--${c}_cache $C/${c}_cache)
done
echo "[$(date)] $SP ($ORG) w/o TEtrimmer on $(hostname): $(grep -c '^>' $RAW) raw seqs; pipeline $(cat $H/COMMIT | head -n 1)"
nextflow run $H/run_reclassify.nf "${OPTS[@]}" -resume
LIB=$(find $LD/out_wo -type f \( -name "$SP.classified_library.fa" -o -name "$SP.te_library.fa" \) -print -quit)
[[ -s $LIB ]] || { echo "no final library in $LD/out_wo" >&2; exit 1; }
mkdir -p $LIB_DIR/$SP/pipeline_wo_tetrimmer_v9
cp $LIB $LIB_DIR/$SP/pipeline_wo_tetrimmer_v9/pipeline_library.fa
echo "[$(date)] w/o TEtrimmer: $(grep -c '^>' $LIB) seqs -> libraries/$SP/pipeline_wo_tetrimmer_v9/"
for t in $(find $R/out_with $LD/out_wo -name "*.dedup_refine.tsv" -o -name "*.tetrimmer_postprocess.tsv" -o -name "*.evidence.tsv" | sort -u); do
    echo "== $t"; awk -F"\t" 'NR>1{print $2}' $t | sort | uniq -c | head
done
echo "[$(date)] done"
