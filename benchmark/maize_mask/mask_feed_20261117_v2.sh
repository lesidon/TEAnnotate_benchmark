#!/bin/bash
# mask_feed.sh - keeps the worker pool of the chunked maize masking topped up (benchmark only).
# Every 10 min: if no item is left (all done or failed 3 times) print the status and exit;
# otherwise submit workers up to HTC_N on htc and MEM_N on mem, never more pending workers than
# claimable items, and leave submit slots free for the user's other jobs (htc QOS: 16 submitted
# jobs per user; mem: 64). Workers for a method whose library is not written yet (the TEAnnotate
# with-TEtrimmer library of the finish job) start once it appears.
set -uo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
M=$B/maize_mask
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
# defaults; feed.conf (re-read every round) may override HTC_N MEM_N HTC_CAP MEM_CAP
HTC_N=5; MEM_N=0
HTC_CAP=12; MEM_CAP=60      # max jobs (array tasks counted) of the user in that partition after a submit
log() { echo "[$(date '+%F %T')] $*"; }
submit() {  # submit <partition>
    if [[ $1 == htc ]]; then
        sbatch --parsable -p htc -c 16 --mem=64G -t 1-00:00:00 --exclude='ct[01-10],gn[19,26]' -J mm_w_htc \
            -o $M/logs/worker_%j.out $M/mask_worker_20261117.sh
    else  # ct02-04 left free for the maize finish job (the only nodes it may use)
        sbatch --parsable -p mem -c 24 --mem=48G -t 3-00:00:00 --exclude='ct[01-04,08]' -J mm_w_mem \
            -o $M/logs/worker_%j.out $M/mask_worker_20261117.sh
    fi
}
log "feeder $SLURM_JOB_ID: htc $HTC_N x 48 CPUs, mem $MEM_N x 24 CPUs (feed.conf may change this)"
while :; do
    [[ -f $M/feed.conf ]] && source $M/feed.conf
    REM=$($PY $M/queue.py remaining)
    if [[ $REM == 0 ]]; then log "no item left"; $PY $M/queue.py status; exit 0; fi
    CL=$($PY $M/queue.py claimable)
    PEND=$(squeue -h -u $USER -t PENDING -n mm_w_htc,mm_w_mem | wc -l)
    for P in htc mem; do
        if [[ $P == htc ]]; then WANT=$HTC_N; CAP=$HTC_CAP; else WANT=$MEM_N; CAP=$MEM_CAP; fi
        HAVE=$(squeue -h -u $USER -n mm_w_$P | wc -l)
        TOT=$(squeue -r -h -u $USER -p $P | wc -l)      # -r: array tasks count toward the QOS limit
        while (( HAVE < WANT && PEND < CL && TOT < CAP )); do
            J=$(submit $P) || { log "sbatch $P failed"; break; }
            log "submitted $P worker $J (claimable $CL, pending $PEND)"
            HAVE=$(( HAVE + 1 )); PEND=$(( PEND + 1 )); TOT=$(( TOT + 1 ))
        done
    done
    sleep 600
done
