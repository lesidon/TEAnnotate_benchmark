#!/bin/bash
# watch_cancel_wo.sh <species> <pipeline job id> - when the run_heldout_v9b*.sh job of
# <species> starts its w/o-TEtrimmer main.nf run, cancel the job: its with-arm is complete
# (masking included) and both final arms come from run_heldout_v10.sh. Exits when the job
# ends for any other reason.
SP=$1; JOB=$2
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
L=$WORK_DIR/heldout_v9/$SP/.nextflow.log
echo "[$(date)] watching $SP (job $JOB)"
while squeue -h -j $JOB 2>/dev/null | grep -q .; do
    if head -n 1 $L 2>/dev/null | grep -q -- '--skip_tetrimmer'; then
        scancel $JOB
        echo "[$(date)] w/o run started; job $JOB cancelled (with arm complete)"; exit 0
    fi
    sleep 60
done
echo "[$(date)] job $JOB ended before its w/o run started"
