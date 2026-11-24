#!/bin/bash
# family_level.sh <species>... - RepeatModeler2-style family evaluation (Flynn et al. 2020):
# mask every method's library with the curated reference library, then classify each
# reference family as perfect / good / present / not found (family_level.py).
# Skips libraries whose masking output is newer than the library snapshot.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
CPUS=${SLURM_CPUS_PER_TASK:-4}
for SP in "$@"; do
    REF=$E/$SP/reference/library.fa
    [[ -s $REF ]] || { echo "no reference for $SP"; continue; }
    for D in $E/$SP/*/; do
        M=$(basename $D); [[ $M == reference ]] && continue
        LIB=$D/library.fa; OUT=$D/famlevel
        [[ -s $LIB ]] || continue
        [[ -s $OUT/library.fa.out && $OUT/library.fa.out -nt $LIB ]] && { echo "$SP/$M: up to date"; continue; }
        rm -rf $OUT; mkdir -p $OUT
        # -libdir: RepeatMasker's general libraries must be buildable on fuse-overlay nodes
        singularity exec --no-home --writable-tmpfs --bind $B --bind $(readlink -f $B) $IMG bash -c "
            cd $OUT && mkdir rmlib && for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
            awk '/^>/{print \">q\" ++i; next} 1' $LIB > library.fa   # RepeatMasker: IDs <= 50 chars
            RepeatMasker -e ncbi -pa $(( CPUS / 4 > 0 ? CPUS / 4 : 1 )) -nolow -no_is -libdir $OUT/rmlib -lib $REF -dir $OUT library.fa" \
            > $OUT/rm.log 2>&1
        echo "$SP/$M: $(tail -n +4 $OUT/library.fa.out 2>/dev/null | wc -l) hits"
    done
    NICE="$(tr '[:lower:]' '[:upper:]' <<< ${SP:0:1})${SP:1}"
    $PY $E/family_level.py $SP $E/$SP $E/results/${NICE}_family_level.csv
done
