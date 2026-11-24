#!/bin/bash
#SBATCH --job-name=maize_finish_resume
#SBATCH --partition=mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --exclude=ct01
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
REC=$B/maize_mask/recovery_20261120/low_cpu
source "$B/env.sh"
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
set -u
export NXF_VER=24.10.6
H=$B/pipeline_head_v13
LD=$(readlink -f "$B/work/heldout_v13/zea_mays")
RAW=$(readlink -f "$B/work/heldout_v9e/zea_mays/raw_v10/zea_mays.raw_library.fa")
PROV=$(readlink -f "$B/work/heldout_v9e/zea_mays/raw_v10/zea_mays.provenance.tsv")
G=$(readlink -f "$B/genomes/plant/zea_mays/zea_mays.fna")
C=$(readlink -f "$B/work/pipeline/drosophila_melanogaster")
CONDA_CACHE=$(readlink -f "$B/work/conda_cache")
SINGULARITY_CACHE=$(readlink -f "$B/work/singularity_heldout")
PRE=$B/tetrimmer_chunks/zea_mays_v13/precomputed.config
SESSION=a7588e4e-ebdd-4249-ba73-325a0b46f57b
test "${SLURM_CPUS_PER_TASK:-0}" -eq 16
test "${SLURM_MEM_PER_NODE:-0}" -ge 131072
for F in "$RAW" "$PROV" "$G" "$PRE"; do test -s "$F"; done
test -d "$LD/.nextflow/cache/$SESSION"
cd "$LD"
echo "$(date -Is) resume cached maize session $SESSION; 16-CPU allocation; cached 64-CPU evidence; postprocess limit48h, allocation72h"
# Do not merge TEtrimmer chunks again: the successful TETRIMMER task is cached.
nextflow run "$H/run_reclassify.nf" -ansi-log false -profile singularity -work-dir "$LD/work" \
    --input "$RAW" --provenance "$PROV" --dedup --genome "$G" --id zea_mays --organism_type plant \
    --max_cpus 32 --max_memory 115GB \
    --conda_cachedir "$CONDA_CACHE" --singularity_cache_dir "$SINGULARITY_CACHE" \
    --pfam_db "$C/pfam_cache/pfam_db" \
    --teclass2_repo "$C/teclass2_cache/repo/TEclass2" --teclass2_model "$C/teclass2_cache/model/teclass2_model" \
    --create_repo "$C/create_cache/repo/CREATE" --create_model "$C/create_cache/model/create_model" \
    --neuralte_repo "$C/neuralte_cache/repo/NeuralTE" \
    --outdir "$LD/out_with" -resume "$SESSION" -c "$PRE" -c "$REC/maize_recovery.config"
LIB=$LD/out_with/zea_mays/05_classification/zea_mays.te_library.fa
test -s "$LIB"
DEST=$B/libraries/zea_mays/pipeline_v13
mkdir -p "$DEST"
cp "$LIB" "$DEST/pipeline_library.fa.tmp.$SLURM_JOB_ID"
mv "$DEST/pipeline_library.fa.tmp.$SLURM_JOB_ID" "$DEST/pipeline_library.fa"
echo "$(date -Is) final maize library published: $(grep -c '^>' "$DEST/pipeline_library.fa") sequences"
