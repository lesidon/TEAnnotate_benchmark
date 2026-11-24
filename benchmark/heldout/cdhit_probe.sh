#!/bin/bash
# Probe: how cd-hit-est's first (serial) batch size depends on -T and -M on the maize
# combined library. Each setting runs 90 s; we only read the batch boundaries it prints.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
L=$WORK_DIR/heldout_v10/zea_mays/work/f7/4aee90f4b33b38618345898face9e3/zea_mays.combined_library.fa
T=$B/heldout/cdhit_probe; mkdir -p $T; cd $T
CDH=$WORK_DIR/conda_cache/env-53d8c082650edb5cdad8cdc7a79b1e38/bin/cd-hit-est
[ -n "$CDH" ] || CDH=$(command -v cd-hit-est)
echo "cd-hit-est: $CDH"; $CDH -h 2>&1 | head -n 1
for set in "16 32768" "16 8000" "16 4000" "16 3000" "32 3000"; do
    set -- $set
    timeout 120 stdbuf -o0 $CDH -i $L -o probe_T$1_M$2.fa -c 0.80 -aS 0.80 -G 0 -g 1 -b 500 -n 5 -d 0 -p 1 -T $1 -M $2 > probe_T$1_M$2.log 2>&1
    echo "T=$1 M=$2: $(grep -a 'comparing' probe_T$1_M$2.log | head -n 3 | tr -s ' ' | tr '\n' ' ') $(grep -a 'Max number of word' probe_T$1_M$2.log)"
done
rm -f probe_*.fa probe_*.fa.clstr
