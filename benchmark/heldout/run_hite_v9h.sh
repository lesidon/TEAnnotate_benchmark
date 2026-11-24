#!/bin/bash
# run_hite_v9h.sh <species> - rerun of a held-out genome's HiTE task with the v9 module
# script and container options exactly as the pipeline ran them, plus --work_dir in the
# task dir. HiTE's scratch defaults to /tmp, which on Zhores is a RAM tmpfs (378 GB): the
# human run filled it (coarse_boundary BLAST ~370 GB) and HiTE exited 0 without a library.
# Output: work/hite_v9h/<sp>/<sp>.hite_library.fa -> libraries/<sp>/hite_v9h/hite_v9h_library.fa
set -eo pipefail
SP="$1"; ORIG="$2"     # ORIG: the pipeline's HiTE task dir (its .command.sh is reused)
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
IMG=$WORK_DIR/singularity_heldout/kanghu-hite-3.3.3.img
D=$WORK_DIR/hite_v9h/$SP; mkdir -p $D; cd $D
ln -sf "$(readlink -f $ORIG/$SP.fna)" $SP.fna
grep -q -- '--annotate 0' $ORIG/.command.sh || { echo "unexpected task script in $ORIG" >&2; exit 1; }
sed 's#--out_dir "\$PWD/\$OUT" \\#--out_dir "$PWD/$OUT" \\\n    --work_dir "$PWD/hite_work" \\#' $ORIG/.command.sh > hite_task.sh
grep -q -- '--work_dir "\$PWD/hite_work"' hite_task.sh || { echo "work_dir patch failed" >&2; exit 1; }
mkdir -p hite_work
echo "[$(date)] $SP: HiTE (v9 module script + --work_dir) on $(hostname)"
singularity exec --no-home --pid -B $(readlink -f $B) -B $B \
    --bind /usr/bin/ps:/usr/bin/ps --bind /lib64/libprocps.so.4:/usr/lib/x86_64-linux-gnu/libprocps.so.4 \
    --writable-tmpfs --bind $B/hite_clone/bin/FiLTR-main/bin/LtrDetector:/HiTE/bin/FiLTR-main/bin/LtrDetector \
    $IMG /bin/bash -c "cd $D; bash -euo pipefail hite_task.sh" > hite.log 2>&1 \
    || { echo "HiTE failed; tail of hite.log:"; tail -20 hite.log; exit 1; }
if grep -aq 'No space left on device' hite.log; then echo "HiTE ran out of disk space" >&2; exit 1; fi
[[ -s $SP.hite_library.fa ]] || { echo "HiTE produced no library" >&2; tail -20 hite.log; exit 1; }
rm -rf hite_work
mkdir -p $LIB_DIR/$SP/hite_v9h
cp $SP.hite_library.fa $LIB_DIR/$SP/hite_v9h/hite_v9h_library.fa
echo "[$(date)] $SP: HiTE library $(grep -c '>' $SP.hite_library.fa) families -> libraries/$SP/hite_v9h/"

