#!/bin/bash
# io_test.sh - is TEtrimmer slow here because of BeeGFS? The same 16 rice sequences (from the
# 40-sequence chunk test, all curated by the task) run twice at the same time on one node, 16
# threads each, set up exactly like a chunk (tt_run_chunk.sh): once with the working/output dir on
# BeeGFS, once on the node's local /tmp (the job's TMPDIR; removed at the end). Inputs (genome,
# BLAST db, Pfam) stay on BeeGFS in both. Records wall time, CPU use (vmstat) and file counts.
set -o pipefail
B=/beegfs/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
X=$B/tetrimmer_chunks; D=$X/io_test
T=$B/ablation/reclass_mm2/oryza_sativa/edta+heliano+hite/work/76/e027dcb0af4f6b89b9f31ecfce5775
LIB=oryza_sativa.gated_library.fa
N=${N:-16}
L=${TMPDIR:-/tmp}/tt_io_$SLURM_JOB_ID
trap 'rm -rf $L' EXIT
mkdir -p $D/beegfs $L/local
awk -v n=$N '/^>/{k++} k<=n' $X/oryza_sativa_mm2/test/$LIB > $D/beegfs/$LIB
cp $D/beegfs/$LIB $L/local/$LIB
echo "[$(date)] $(hostname), $SLURM_CPUS_PER_TASK CPUs, $(grep -c '>' $D/beegfs/$LIB) sequences; local dir $L ($(df -h ${TMPDIR:-/tmp} | tail -1 | awk '{print $1, $2}'))"
vmstat -t 30 > $D/vmstat.txt & VM=$!
run() {  # run <dir> <label>
    local s=$(date +%s)
    /usr/bin/time -o $D/time_$2.txt -f '%e s wall, %U s user, %S s sys (CPU busy = (user+sys)/wall)' \
        env -u SLURM_MEM_PER_NODE TT_THREADS=16 bash $X/tt_run_chunk.sh $T $1 > $D/run_$2.log 2>&1
    local rc=$?
    echo "$2: rc=$rc wall=$(( $(date +%s) - s )) s, $(cat $D/time_$2.txt), files=$(find $1 -type f | wc -l), $(tail -n 1 $D/run_$2.log)" >> $D/times.txt
}
run $D/beegfs beegfs & P1=$!
run $L/local local & P2=$!
wait $P1 $P2
kill $VM
cp $L/local/tetrimmer_run.log $D/local_tetrimmer_run.log 2>/dev/null
cp $L/local/oryza_sativa_tetrimmer/summary.txt $D/local_summary.txt 2>/dev/null
cat $D/times.txt
echo "[$(date)] done"
