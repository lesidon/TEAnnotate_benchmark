#!/bin/bash
#SBATCH --job-name=cbri_ensemble_off
#SBATCH --partition=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --exclude=ct[01-10],gn[19,26]
set -euo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
R=$B/edta_manual_bench/ensemble_ablation_20261123
SP=caenorhabditis_briggsae
OUT=$R/$SP/ensemble_off
G=$B/genomes/animal/$SP/$SP.fna
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
test -s "$OUT/library.fa"
test ! -e "$OUT/$SP.fna.out"
singularity exec --no-home --writable-tmpfs --bind "$B" "$IMG" bash -euc '
    out=$1; genome=$2
    export HOME="$out"
    cd "$out"
    mkdir rmlib
    for f in /opt/RepeatMasker/Libraries/*; do ln -s "$f" rmlib/; done
    RepeatMasker -e ncbi -pa 4 -q -no_is -norna -nolow -div 40 -libdir "$out/rmlib" \
        -lib "$out/library.fa" -cutoff 225 -dir "$out" "$genome"' _ "$OUT" "$G" > "$OUT/rm.log" 2>&1
test -s "$OUT/$SP.fna.out"
grep -qm1 '^[[:space:]]*[0-9]' "$OUT/$SP.fna.out"
