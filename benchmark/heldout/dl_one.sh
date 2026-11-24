#!/bin/bash
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
H=$B/heldout
log() { echo "[$(date '+%F %T')] $*"; }
source $H/dl_genome.inc
dl_genome "$@"
