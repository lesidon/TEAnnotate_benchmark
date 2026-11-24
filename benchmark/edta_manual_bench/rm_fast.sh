#!/bin/bash
# rm_fast.sh <species> <method> <wait-job> - rm_one.sh's RepeatMasker run with -pa = CPUs/4, in a
# separate folder (<method>_fast), for a held-out masking that its 24-CPU, six-at-once bench job
# will not finish in time. When done, wait for <wait-job> to leave the queue; then, if the
# original folder still has no .out, install this one there and, when every masking of the
# genome is done, run the scoring (heldout_bench_v10.sh masks nothing more in that case).
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
SP=$1; M=$2; WAITJOB=$3
PA=$(( ${SLURM_CPUS_PER_TASK:-16} / 4 ))
KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
G=$B/genomes/$KING/$SP/$SP.fna
SRC=$E/$SP/$M; OUT=$E/$SP/${M}_fast
mkdir -p $OUT; cp $SRC/library.fa $OUT/library.fa
echo "[$(date)] $SP/$M on $(hostname): -pa $PA, $(grep -c '>' $OUT/library.fa) seqs, genome=$G"
ok=0
for attempt in 1 2; do
    rm -rf "$OUT"/$SP.fna.{out,tbl,masked,cat.gz} "$OUT/rmlib"
    singularity exec --no-home --writable-tmpfs --bind "$B" "$IMG" bash -c "
        export HOME='$OUT'; cd '$OUT'
        mkdir rmlib; for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        RepeatMasker -e ncbi -pa $PA -q -no_is -norna -nolow -div 40 -libdir '$OUT/rmlib' \
            -lib '$OUT/library.fa' -cutoff 225 -dir '$OUT' '$G'" > "$OUT/rm.log" 2>&1 || true
    if [[ -s $OUT/$SP.fna.out ]] && [[ $(tail -n +4 "$OUT/$SP.fna.out" | wc -l) -gt 0 ]]; then
        echo "[$(date)] OK attempt $attempt: $(tail -n +4 $OUT/$SP.fna.out | wc -l) hits"; ok=1; break
    fi
    echo "[$(date)] attempt $attempt produced no .out; RepeatMasker log tail:"; tail -5 "$OUT/rm.log"
done
(( ok )) || { echo "FAILED: $SP/$M"; exit 1; }

# sacct, not squeue: squeue -j fails once a finished job is purged, and a failed query must
# never read as "ended" while the bench job may still be writing into $SRC
while :; do
    st=$(sacct -X -n -j $WAITJOB -o State 2>/dev/null | head -n 1 | tr -d ' ') || st=
    [[ -n $st && ! $st =~ ^(RUNNING|PENDING|REQUEUED|SUSPENDED|RESIZING|COMPLETING) ]] && break
    sleep 600
done
echo "[$(date)] job $WAITJOB has ended ($st)"
if [[ -s $SRC/$SP.fna.out ]]; then
    echo "$SRC already has its own .out (the bench job finished it): nothing to install"; exit 0
fi
cp $OUT/$SP.fna.{out,tbl,masked,cat.gz} $SRC/
echo "[$(date)] installed $OUT/$SP.fna.out into $SRC"
missing=$(for d in $E/$SP/*/; do [[ $d == *_fast/ ]] && continue; [[ -s $d/$SP.fna.out ]] || basename $d; done)
if [[ -n $missing ]]; then echo "not scoring yet, maskings still missing: $missing"; exit 0; fi
echo "[$(date)] every masking of $SP is done: scoring"
bash $E/heldout_bench_v10.sh $SP
