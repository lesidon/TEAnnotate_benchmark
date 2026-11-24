#!/bin/bash
# run_fungus.sh <species> - fungal pipeline run with the original checkout + config/rerun_fungus.<sp>.config.
# Same invocation and shared caches (env.sh) as the July runs, so zymoseptoria_tritici
# reuses its July EDTA/EarlGrey tasks; septoria_linicola runs from scratch.
set -eo pipefail
SP="$1"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
PIPE=/trinity/home/leonid.sidorov/tmp/nf-core-te-annotation   # original checkout, read-only
CFG=$B/config/rerun_fungus.$SP.config
source $B/env.sh
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
conda activate nf-env
RUN=$WORK_DIR/pipeline/$SP
mkdir -p "$RUN"
KING=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
G=$GENOME_DIR/$KING/$SP/$SP.fna
[[ -s $RUN/samplesheet.csv ]] || printf 'sample,genome,organism_type,gff\n%s,%s,fungus,\n' "$SP" "$(readlink -f "$G")" > "$RUN/samplesheet.csv"
if [[ -d $RUN/out && ! -e $RUN/out.pre_hitefix_2026-09-28 ]]; then mv "$RUN/out" "$RUN/out.pre_hitefix_2026-09-28"; fi
echo "[$(date)] $SP on $(hostname): $PIPE + $CFG"
cd "$RUN"
nextflow run "$PIPE/main.nf" -profile singularity \
    -work-dir "$RUN/work" \
    --input "$RUN/samplesheet.csv" \
    --outdir "$RUN/out" \
    --max_cpus 16 --max_memory 256GB \
    -c "$CFG" \
    -resume
LIB=$(find "$RUN/out" -type f -name '*.classified_library.fa' -print -quit)
[[ -s $LIB ]] || { echo "no classified library under $RUN/out" >&2; exit 1; }
mkdir -p "$LIB_DIR/$SP/pipeline"
cp "$LIB" "$LIB_DIR/$SP/pipeline/pipeline_library.fa"
echo "[$(date)] captured $(grep -c '^>' "$LIB") sequences -> libraries/$SP/pipeline/pipeline_library.fa"
H=$(find "$RUN/out" -name "$SP.hite_library.fa" -print -quit)
echo "HiTE library: $(grep -c '^>' "$H") families, LTR $(grep '^>' "$H" | grep -ci ltr)"
