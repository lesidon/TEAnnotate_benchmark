#!/bin/bash
# mask_worker.sh - one worker of the chunked maize RepeatMasker run (benchmark only; see README).
# Claims <method, chunk> items (queue.py) while an item's estimated run time fits in the job's
# remaining wall time, and masks the chunk exactly as edta_manual_bench/rm_one.sh masks a genome
# (same RepeatMasker options and image; -pa = CPUs / 4). Keeps only the .out/.tbl of a chunk.
# The worker that finishes the last chunk of a method merges the method (merge_out.py).
set -uo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
M=$B/maize_mask
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
CPUS=${SLURM_CPUS_PER_TASK:-4}; PA=$(( CPUS / 4 > 0 ? CPUS / 4 : 1 ))
JOB=${SLURM_JOB_ID:?run under Slurm}
log() { echo "[$(date '+%F %T')] $*"; }
left() {  # seconds of wall time left in this job
    local L d=0 h=0 m=0 s=0
    L=$(squeue -h -j "$JOB" -o %L 2>/dev/null | tr -d ' ')
    [[ -z $L || $L == UNLIMITED || $L == INVALID ]] && { echo 0; return; }
    [[ $L == *-* ]] && { d=${L%%-*}; L=${L#*-}; }
    IFS=: read -r a b c <<< "$L"
    if [[ -n ${c:-} ]]; then h=$a; m=$b; s=$c; elif [[ -n ${b:-} ]]; then m=$a; s=$b; else s=$a; fi
    echo $(( 10#$d * 86400 + 10#$h * 3600 + 10#$m * 60 + 10#$s ))
}
merge_ready() {
    for m in $($PY $M/queue.py merge_ready); do
        mkdir "$M/merging/$m" 2>/dev/null || continue
        if $PY $M/merge_out.py "$m" > "$M/logs/merge_$m.log" 2>&1; then log "merged $m: $(tail -n 1 $M/logs/merge_$m.log)"
        else log "merge of $m FAILED (logs/merge_$m.log)"; rmdir "$M/merging/$m"; fi
    done
}
mkdir -p $M/merging $M/logs
log "worker $JOB on $(hostname), $CPUS CPUs (-pa $PA)"
while :; do
    MTH= CH= LIB= FA= EST=
    read -r MTH CH LIB FA EST < <($PY $M/queue.py claim "$JOB" "$CPUS" "$(left)")
    [[ -n $MTH ]] || { log "nothing claimable that fits the $(left) s left"; break; }
    D=$M/runs/$MTH/$CH
    log "start $MTH/$CH (estimate $(( EST / 60 )) min)"
    t0=$(date +%s)
    rm -rf "$D"; mkdir -p "$D"
    singularity exec --no-home --writable-tmpfs --bind "$B" --bind "$(readlink -f "$B")" "$IMG" bash -c "
        export HOME='$D'; cd '$D'
        mkdir rmlib; for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        RepeatMasker -e ncbi -pa $PA -q -no_is -norna -nolow -div 40 -libdir '$D/rmlib' \
            -lib '$LIB' -cutoff 225 -dir '$D' '$FA'" > "$D/rm.log" 2>&1
    RC=$?
    dt=$(( $(date +%s) - t0 ))
    O=$D/$CH.fa.out
    if [[ $RC -eq 0 && -s $O ]] && { [[ $(tail -n +4 "$O" | wc -l) -gt 0 ]] || grep -q 'no repetitive sequences' "$O"; }; then
        rm -rf "$D"/rmlib "$D"/*.masked "$D"/*.cat.gz "$D"/*.cat "$D"/RM_*
        $PY $M/queue.py done "$MTH" "$CH" "$dt" "$CPUS"
        log "done $MTH/$CH in $(( dt / 60 )) min: $(tail -n +4 "$O" | wc -l) hits"
        merge_ready
    else
        $PY $M/queue.py fail "$MTH" "$CH"
        log "FAILED $MTH/$CH after $(( dt / 60 )) min; RepeatMasker log tail:"; tail -n 5 "$D/rm.log"
    fi
done
merge_ready
log "worker $JOB exits"
