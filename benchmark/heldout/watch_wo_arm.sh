#!/bin/bash
# watch_wo_arm.sh <species> <pipeline job id> - wait until the run_heldout_v9b*.sh job of
# <species> starts its w/o-TEtrimmer main.nf run. If that run re-submits EARLGREY instead
# of resuming it (the zebrafish failure), cancel the job and submit run_heldout_wo_v9.sh
# (library stages on the with-arm's raw library). If Earl Grey is resumed from the
# cache, leave the job alone. Exits when the job ends for any other reason.
SP=$1; JOB=$2
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
R=$WORK_DIR/heldout_v9/$SP; L=$R/.nextflow.log
BAD=$(cat $B/heldout/exclude_nodes.txt)
echo "[$(date)] watching $SP (job $JOB)"
while squeue -h -j $JOB 2>/dev/null | grep -q .; do
    if head -n 1 $L 2>/dev/null | grep -q -- '--skip_tetrimmer'; then
        if grep -q 'Cached process > .*EARLGREY' $L; then
            echo "[$(date)] w/o arm resumed Earl Grey from the cache; job left running"; exit 0
        fi
        if grep -q 'Submitted process > .*EARLGREY' $L; then
            scancel $JOB
            echo "[$(date)] w/o arm re-submitted Earl Grey; job $JOB cancelled"
            sleep 30
            j=$(sbatch --parsable -p mem -c 32 --mem=96G -t 2-00:00:00 --exclude=$BAD,ct02 -J ho_wo_${SP:0:3} \
                -o $B/logs/ho_wo_${SP}_%j.out $B/heldout/run_heldout_wo_v9.sh $SP)
            echo "[$(date)] run_heldout_wo_v9.sh $SP submitted as $j"; exit 0
        fi
    fi
    sleep 60
done
echo "[$(date)] job $JOB ended before its w/o arm submitted Earl Grey; nothing done"
