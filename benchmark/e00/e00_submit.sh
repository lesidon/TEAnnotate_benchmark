#!/bin/bash
# e00_submit.sh [genome...] - E00 job chain. Per genome: A = stage libraries + family-level
# masking + families (oracle libraries); G = genome masking of raw_union, refined,
# oracle_wo, oracle_with (array, after A). Then base level + report after every G.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
X=$B/e00; PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
SPECIES="${*:-septoria_linicola zymoseptoria_tritici drosophila_melanogaster caenorhabditis_elegans arabidopsis_thaliana oryza_sativa}"
mkdir -p $X/logs
DEPS=()
for SP in $SPECIES; do
    case $SP in
        oryza_sativa)  CA=16; MA=32G; TA=12:00:00; CG=32; MG=64G; TG=1-00:00:00 ;;
        arabidopsis_thaliana|caenorhabditis_elegans|drosophila_melanogaster)
                       CA=8;  MA=16G; TA=4:00:00;  CG=16; MG=32G; TG=12:00:00 ;;
        *)             CA=8;  MA=16G; TA=2:00:00;  CG=16; MG=32G; TG=6:00:00 ;;
    esac
    JA=$(sbatch --parsable -p htc -c $CA --mem=$MA -t $TA -J e00a_${SP:0:8} -o $X/logs/a_${SP}_%j.out \
        --wrap "set -e; $PY $X/e00.py collect $SP; bash $X/e00_run.sh famlevel $SP; $PY $X/e00.py families $SP")
    JG=$(sbatch --parsable -p htc -c $CG --mem=$MG -t $TG --array=0-3 -d afterok:$JA -J e00g_${SP:0:8} \
        -o $X/logs/g_${SP}_%A_%a.out \
        --wrap "L=(raw_union refined oracle_wo oracle_with); bash $X/e00_run.sh genome $SP \${L[\$SLURM_ARRAY_TASK_ID]}")
    DEPS+=($JG)
    echo "$SP: A $JA, genome maskings $JG"
done
JB=$(sbatch --parsable -p htc -c 2 --mem=24G -t 4:00:00 -d afterany:$(IFS=:; echo "${DEPS[*]}") -J e00_base \
    -o $X/logs/base_%j.out --wrap "for SP in $SPECIES; do $PY $X/e00.py baselevel \$SP; done; $PY $X/e00.py report")
echo "base level + report: $JB"
