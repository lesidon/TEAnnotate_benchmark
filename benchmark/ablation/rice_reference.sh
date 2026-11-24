#!/bin/bash
# rice_reference.sh - reference annotation of rice for the benchmark (rice moved from the
# held-out to the development set): the curated EDTA-benchmark library masks the genome with
# the same RepeatMasker settings as every test library (rm_one.sh), and is the family-level
# reference (edta_manual_bench/oryza_sativa/reference/library.fa).
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; SP=oryza_sativa
L=$B/curated/$SP/$SP.curated.lib.fa
mkdir -p $E/$SP/reference $E/logs
cp $L $E/$SP/reference/library.fa
cd $E
n=$(grep -nP "^$SP\treference\t" specs.tsv | head -1 | cut -d: -f1 || true)   # no row yet: grep exits 1
if [[ -z $n ]]; then printf '%s\t%s\t%s\n' $SP reference $L >> specs.tsv; n=$(grep -c . specs.tsv); fi
echo "[$(date)] $SP reference = specs row $(( n - 1 )) ($(grep -c '>' $L) families)"
SLURM_ARRAY_TASK_ID=$(( n - 1 )) bash $E/rm_one.sh
