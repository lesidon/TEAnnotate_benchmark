#!/bin/bash
# e00_submit_rest.sh <base-dependency job ids, ':'-separated> <genome>...
# The htc QOS takes 16 queued jobs per user, so here each genome's four genome maskings
# run inside ONE job, side by side (CPUs split four ways), after that genome's A job
# when it is still queued or running. Then the base level + report job.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
X=$B/e00; PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
PREV=$1; shift
DEPS=(${PREV//:/ })
for SP in "$@"; do
    case $SP in
        oryza_sativa) C=64; M=96G; T=1-00:00:00 ;;
        *)            C=32; M=48G; T=12:00:00 ;;
    esac
    JA=$(squeue -u leonid.sidorov -h -n e00a_${SP:0:8} -o %i | head -1)
    D=${JA:+-d afterok:$JA}
    JG=$(sbatch --parsable -p htc -c $C --mem=$M -t $T $D -J e00g_${SP:0:8} -o $X/logs/g_${SP}_%j.out --wrap "
        for L in raw_union refined oracle_wo oracle_with; do
            SLURM_CPUS_PER_TASK=$(( C / 4 )) bash $X/e00_run.sh genome $SP \$L > $X/logs/g_${SP}_\$L.out 2>&1 &
        done
        wait
        grep -q FAILED $X/logs/g_${SP}_{raw_union,refined,oracle_wo,oracle_with}.out && exit 1 || exit 0")
    [[ -n $JG ]] || { echo "submission failed for $SP"; exit 1; }
    DEPS+=($JG)
    echo "$SP: A ${JA:-done}, genome maskings $JG"
done
ALL="septoria_linicola zymoseptoria_tritici drosophila_melanogaster caenorhabditis_elegans arabidopsis_thaliana oryza_sativa"
JB=$(sbatch --parsable -p htc -c 2 --mem=24G -t 4:00:00 -d afterany:$(IFS=:; echo "${DEPS[*]}") -J e00_base \
    -o $X/logs/base_%j.out --wrap "for SP in $ALL; do $PY $X/e00.py baselevel \$SP; done; $PY $X/e00.py report")
echo "base level + report: $JB"
