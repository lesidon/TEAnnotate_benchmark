#!/bin/bash
# run_tool.sh <species> <tool> [stub] - one discovery tool on one development genome
# (discover_one.nf, v11 modules = v10 + the EDTA fix 16bcdbb). Library ->
# ablation/libs/<sp>/<tool>.fa; the Nextflow trace (realtime, cpus, %cpu, peak_rss) ->
# ablation/libs/<sp>/<tool>.trace.tsv. A tool that ran and found nothing (e.g. HELIANO
# on Drosophila: 0 Helitron blocks) leaves an empty library, which is a valid result.
set -eo pipefail
SP="$1"; TOOL="$2"; MODE="${3:-}"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
A=$B/ablation; H=$B/pipeline_head_v11
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
ROW=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s')
KING=$(cut -f1 <<< "$ROW"); ORG=$(cut -f5 <<< "$ROW")
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
[[ -s $G ]] || { echo "missing genome for $SP" >&2; exit 1; }
R=$A/runs/$SP/$TOOL; [[ $MODE == stub ]] && R=$A/stub/$SP/$TOOL
mkdir -p $R; cd $R
CPUS=${SLURM_CPUS_PER_TASK:-16}
MEM=$(( ${SLURM_MEM_PER_NODE:-65536} * 9 / 10 / 1024 ))GB
OPTS=(-c $H/nextflow.config -c $A/ablation.config -profile singularity -work-dir $R/work
      --tool $TOOL --genome $G --id $SP --organism_type $ORG --outdir $R/out
      --max_cpus $CPUS --max_memory $MEM
      --hite_ltrdetector_dir $B/hite_clone/bin/FiLTR-main/bin/LtrDetector
      --singularity_cache_dir $WORK_DIR/singularity_heldout)
[[ $MODE == stub ]] && OPTS+=(-stub)
echo "[$(date)] $SP ($ORG) $TOOL on $(hostname): $CPUS cpus, $MEM; modules of $(cat $H/COMMIT)"
nextflow run $A/discover_one.nf "${OPTS[@]}" -with-trace $R/trace.tsv -resume
LIB=$(find $R/work -name "$SP.${TOOL}_library.fa" -print -quit)
[[ -n $LIB ]] || { echo "no $TOOL library" >&2; exit 1; }
[[ -s $LIB ]] || echo "WARN: $TOOL found no TE in $SP (empty library; a valid result)" >&2
[[ $MODE == stub ]] && { echo "stub OK: $LIB"; exit 0; }
mkdir -p $A/libs/$SP
cp $LIB $A/libs/$SP/$TOOL.fa
cp $R/trace.tsv $A/libs/$SP/$TOOL.trace.tsv
echo "[$(date)] $TOOL: $(grep -c '^>' $A/libs/$SP/$TOOL.fa) seqs"; cat $R/trace.tsv
