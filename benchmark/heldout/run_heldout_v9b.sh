#!/bin/bash
# run_heldout_v9b.sh <species> (v9 + HiTE at 16 CPUs) [stub] - the frozen v9 pipeline (pipeline_head_v9/main.nf,
# full run: de novo -> dedup -> evidence gate -> curation -> classification -> masking) on a
# held-out genome, twice from one launch dir so the second run reuses the de novo stage:
#   with TEtrimmer  -> libraries/<sp>/pipeline_v9/pipeline_library.fa
#   --skip_tetrimmer -> libraries/<sp>/pipeline_wo_tetrimmer_v9/pipeline_library.fa
# The route's own tool libraries (HiTE + EarlGrey, or HiTE + EDTA + HELIANO on the plant
# route) are copied to libraries/<sp>/<tool>/<tool>_library.fa as the standalone baselines
# (same modules, same CLI). Re-entrant: rerun the same command to resume after a timeout.
# 'stub' as 2nd argument runs nextflow -stub into a throwaway launch dir (wiring test).
set -eo pipefail
SP="$1"; MODE="${2:-}"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v9
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
[[ -n $ROW ]] || { echo "no genomes.tsv row for $SP" >&2; exit 1; }
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
[[ -s $G ]] || { echo "missing genome $G" >&2; exit 1; }
R=$WORK_DIR/heldout_v9/$SP
[[ $MODE == stub ]] && R=$WORK_DIR/heldout_v9_stub/$SP
mkdir -p $R; cd $R
printf 'sample,genome,organism_type,gff\n%s,%s,%s,\n' "$SP" "$G" "$ORG" > samplesheet.csv
C=$WORK_DIR/pipeline/drosophila_melanogaster          # species-independent model/DB caches
CPUS=${SLURM_CPUS_PER_TASK:-32}
MEM=$(( ${SLURM_MEM_PER_NODE:-163840} * 9 / 10 / 1024 ))GB
OPTS=(-profile singularity -work-dir $R/work --input $R/samplesheet.csv
      --max_cpus $CPUS --max_memory $MEM
      --hite_ltrdetector_dir $B/hite_clone/bin/FiLTR-main/bin/LtrDetector
      --singularity_cache_dir $WORK_DIR/singularity_heldout
      --conda_cachedir $WORK_DIR/conda_cache -c $B/heldout/heldout_hite16.config)
for c in pfam teclass2 deepte create neuralte; do
    [[ -d $C/${c}_cache ]] && OPTS+=(--${c}_cache $C/${c}_cache)
done
[[ $MODE == stub ]] && OPTS+=(-stub)
echo "[$(date)] $SP ($ORG, $(grep -c '>' $G) seqs) on $(hostname), $CPUS cpus, $MEM; pipeline $(cat $H/COMMIT)"
nextflow run $H/main.nf "${OPTS[@]}" --outdir $R/out_with -resume
nextflow run $H/main.nf "${OPTS[@]}" --skip_tetrimmer --outdir $R/out_wo -resume
[[ $MODE == stub ]] && { echo "[$(date)] stub run OK"; find $R/out_with $R/out_wo -name "*.classified_library.fa"; exit 0; }
for arm in with:pipeline_v9 wo:pipeline_wo_tetrimmer_v9; do
    LIB=$(find $R/out_${arm%%:*} -type f -name "$SP.classified_library.fa" -print -quit)
    [[ -s $LIB ]] || { echo "no classified library in out_${arm%%:*}" >&2; exit 1; }
    mkdir -p $LIB_DIR/$SP/${arm#*:}
    cp $LIB $LIB_DIR/$SP/${arm#*:}/pipeline_library.fa
    echo "[$(date)] ${arm%%:*} TEtrimmer: $(grep -c '^>' $LIB) seqs -> libraries/$SP/${arm#*:}/"
done
for t in hite earlgrey edta heliano; do
    L=$(find $R/out_with/$SP/02_denovo_te -name "$SP.${t}_library.fa" -print -quit 2>/dev/null || true)
    if [[ -n $L && -s $L ]]; then
        mkdir -p $LIB_DIR/$SP/$t; cp $L $LIB_DIR/$SP/$t/${t}_library.fa
        echo "[$(date)] $t: $(grep -c '^>' $L) seqs -> libraries/$SP/$t/"
    fi
done
for t in $(find $R/out_with $R/out_wo -name "*.dedup_refine.tsv" -o -name "*.tetrimmer_postprocess.tsv" -o -name "*.evidence.tsv" | sort -u); do
    echo "== ${t#$R/}"; awk -F"\t" 'NR>1{print $2}' $t | sort | uniq -c | head
done
echo "[$(date)] done"

