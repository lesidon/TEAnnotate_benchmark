#!/bin/bash
# run_heldout_v13_maize.sh prep|finish - maize on pipeline_head_v13 (as run_heldout_v13.sh), with
# TEtrimmer finished in chunks outside the pipeline (tetrimmer_chunks/, benchmark only).
#   prep    the w/o-TEtrimmer arm first (-> libraries/zea_mays/pipeline_wo_tetrimmer_v13/), then the
#           with arm until its TETRIMMER task has started (genome BLAST db built, first sequence
#           done); stop it, deal the remaining sequences into chunks, submit the chunk array
#           (htc) and "finish" after it.
#   finish  merge the chunks, resume the with arm (the TETRIMMER task finishes the merged output
#           with --continue_analysis) -> libraries/zea_mays/pipeline_v13/.
# TE_EVIDENCE runs on 64 CPUs here (launch-dir config): BLAST, TEsorter and hmmsearch give the
# same results on any thread number. Same launch dir and checkout for both phases (-resume).
set -eo pipefail
PHASE=$1; SP=zea_mays
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H=$B/pipeline_head_v13; X=$B/tetrimmer_chunks; R=$X/zea_mays_v13
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
RAWD=$WORK_DIR/heldout_v9e/$SP/raw_v10
RAW=$RAWD/$SP.raw_library.fa; PROV=$RAWD/$SP.provenance.tsv
[[ -s $RAW && -s $PROV && -s $G ]] || { echo "missing raw library, provenance or genome" >&2; exit 1; }
C=$WORK_DIR/pipeline/drosophila_melanogaster
LD=$WORK_DIR/heldout_v13/$SP
mkdir -p $LD $R/logs; cd $LD
cat > $LD/nextflow.config <<'EOF'
// launch-dir overrides: maize raw library ~ the largest held-out input; long limits for the
// library stages, the local executor sees the whole allocation (2026-10-26), and TE_EVIDENCE on
// 64 CPUs (2026-11-05; its BLAST/TEsorter/hmmsearch results do not depend on the thread number)
process {
    withName: 'CD_HIT_EST'   { time = 5.d }
    withName: 'DEDUP_REFINE' { time = 5.d }
    withName: 'TE_EVIDENCE'  { time = 5.d; cpus = 64 }
}
executor { cpus = (System.getenv('SLURM_CPUS_PER_TASK') ?: '16') as int }
EOF
OPTS=(-profile singularity -work-dir $LD/work
      --input $RAW --provenance $PROV --dedup --genome $G --id $SP --organism_type $ORG
      --max_cpus ${SLURM_CPUS_PER_TASK:-16}
      --max_memory $(( ${SLURM_MEM_PER_NODE:-65536} * 9 / 10 / 1024 ))GB
      --conda_cachedir $WORK_DIR/conda_cache --singularity_cache_dir $WORK_DIR/singularity_heldout
      --pfam_db $C/pfam_cache/pfam_db
      --teclass2_repo $C/teclass2_cache/repo/TEclass2 --teclass2_model $C/teclass2_cache/model/teclass2_model
      --create_repo $C/create_cache/repo/CREATE --create_model $C/create_cache/model/create_model
      --neuralte_repo $C/neuralte_cache/repo/NeuralTE)
publish() {  # publish <arm> <library dir name>
    local lib=$LD/out_$1/$SP/05_classification/$SP.te_library.fa
    [[ -s $lib ]] || { echo "no final library in out_$1" >&2; exit 1; }
    mkdir -p $LIB_DIR/$SP/$2; cp $lib $LIB_DIR/$SP/$2/pipeline_library.fa
    echo "[$(date)] $1 TEtrimmer: $(grep -c '^>' $lib) seqs -> libraries/$SP/$2/"
}
tt_task() {  # the with arm's TETRIMMER task dir (the only task writing <id>_tetrimmer)
    { grep -lE -- "--output_dir +${SP}_tetrimmer" $LD/work/*/*/.command.sh 2>/dev/null || true; } | head -n 1 | xargs -r dirname
}
echo "[$(date)] $PHASE on $(hostname), $SLURM_CPUS_PER_TASK CPUs; pipeline $(head -n 1 $H/COMMIT)"

if [[ $PHASE == prep ]]; then
    echo "[$(date)] $(grep -c '^>' $RAW) raw seqs ($RAW)"
    nextflow run $H/run_reclassify.nf "${OPTS[@]}" --skip_tetrimmer --outdir $LD/out_wo -resume
    publish wo pipeline_wo_tetrimmer_v13
    grep -h 'MMseqs2 instead' $LD/work/*/*/.command.err 2>/dev/null || true
    nextflow run $H/run_reclassify.nf "${OPTS[@]}" --outdir $LD/out_with -resume > $LD/with_until_tetrimmer.log 2>&1 &
    NF=$!; T=
    while kill -0 $NF 2>/dev/null; do
        T=$(tt_task)
        if [[ -n $T && $(tail -n +2 $T/${SP}_tetrimmer/summary.txt 2>/dev/null | wc -l) -gt 0 ]]; then break; fi
        sleep 120
    done
    kill -0 $NF 2>/dev/null || { echo "the with run ended before TEtrimmer started" >&2; tail -20 $LD/with_until_tetrimmer.log >&2; exit 1; }
    echo "[$(date)] TETRIMMER started in $T: stopping the with run"
    kill -TERM $NF || true; wait $NF || true; sleep 60
    pkill -u $USER -f "$(readlink -f $T)" || true
    $PY $X/tt_chunks.py prepare $T $R ${TT_CHUNKS:-24}
    A=$(TT_PARALLEL=${TT_PARALLEL:-6} TT_THREADS=32 bash $X/tt_submit.sh $T $R tt_maize -p htc -c 24 --mem=200G -t 1-00:00:00)
    echo "[$(date)] chunk array $A"
    F=$(sbatch --parsable -p mem -c 32 --mem=250G -t 3-00:00:00 --exclude='ct[01,05-10]' -d afterok:$A \
        -J tt_maize_finish -o $R/logs/finish_%j.out --wrap "bash $B/heldout/run_heldout_v13_maize.sh finish")
    echo "[$(date)] finish job $F (after $A)"
elif [[ $PHASE == finish ]]; then
    T=$(tt_task)
    $PY $X/tt_chunks.py merge $T $R
    # containers run with --no-home and bind only /beegfs/...: pass the resolved path (16 Nov)
    echo "params.tetrimmer_precomputed = '$(readlink -f $R)/merged/${SP}_tetrimmer'" > $R/precomputed.config
    nextflow run $H/run_reclassify.nf "${OPTS[@]}" --outdir $LD/out_with -resume -c $R/precomputed.config
    publish with pipeline_v13
fi
echo "[$(date)] $PHASE done"
