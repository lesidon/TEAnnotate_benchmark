#!/bin/bash
# run_edta_v9.sh <species> - standalone EDTA baseline for a held-out genome whose pipeline
# route has no EDTA (animal/fungus): EDTA 2.3.0 run with the v9 EDTA module script
# (edta_module_16d03c3.sh: patched TIR-Learner, race-free cleanup_nested.pl, rice-fallback
# guard; the module is unchanged from 31e2e6f to 3ffb6ae). Output:
#   work/edta_fixed/<sp>/<sp>.edta_library.fa -> libraries/<sp>/edta2/edta2_library.fa
set -eo pipefail
SP="$1"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
IMG=$WORK_DIR/singularity_cache/quay.io-biocontainers-edta-2.3.0--hdfd78af_0.img
KING=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
[[ -s $G ]] || { echo "missing genome for $SP" >&2; exit 1; }
D=$WORK_DIR/edta_v9e/$SP; mkdir -p $D; cd $D
ln -sf $G $SP.fna
CPUS=${SLURM_CPUS_PER_TASK:-16}
sed -e "s/__SP__/$SP/g" -e "s/__CPUS__/$CPUS/g" $B/edta_module_16d03c3.sh > edta_task.sh
echo "[$(date)] $SP: EDTA (v9e module script, EDTA seq-ID fix) on $(hostname), $CPUS cpus"
singularity exec --no-home --bind $B --bind $(readlink -f $B) $IMG bash -eo pipefail edta_task.sh > edta.log 2>&1 \
    || { echo "EDTA failed; tail of edta.log:"; tail -20 edta.log; exit 1; }
grep -q '^lock(%seq)' edta_src/bin/cleanup_nested.pl || { echo "race patch missing in edta_src" >&2; exit 1; }
mkdir -p $LIB_DIR/$SP/edta2
cp $SP.edta_library.fa $LIB_DIR/$SP/edta2/edta2_library.fa
echo "[$(date)] $SP: EDTA library $(grep -c '>' $SP.edta_library.fa) families -> libraries/$SP/edta2/"

