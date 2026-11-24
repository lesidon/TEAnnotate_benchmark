#!/bin/bash
# edta_rerun_bench.sh - race-free EDTA (pipeline 31e2e6f) for the four benchmark species:
# EDTA filter + final steps, then re-mask, lib-test and family-level for the edta2 rows.
# Racy results stay in work/edta_fixed_racy/<sp>/ and edta_manual_bench/racy_edta2/<sp>/.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
SPECIES="zymoseptoria_tritici septoria_linicola drosophila_melanogaster arabidopsis_thaliana"
declare -A ROW=([zymoseptoria_tritici]=40 [septoria_linicola]=43 [drosophila_melanogaster]=46 [arabidopsis_thaliana]=49)
for sp in $SPECIES; do   # guard: the specs rows must be this species' edta2 row
    read -r s m l < <(sed -n "$(( ${ROW[$sp]} + 1 ))p" $E/specs.tsv)
    [[ $s == $sp && $m == edta2 ]] || { echo "specs row ${ROW[$sp]} is $s/$m, not $sp/edta2"; exit 1; }
done
echo "[$(date)] EDTA filter+final on $(hostname)"
declare -A PID
for sp in $SPECIES; do
    CPUS=8 bash $B/rerun_edta_filter.sh $sp > $E/logs/edta_rerun_$sp.log 2>&1 &
    PID[$sp]=$!
done
OK=""
for sp in $SPECIES; do
    if wait ${PID[$sp]}; then OK="$OK $sp"; tail -n 3 $E/logs/edta_rerun_$sp.log
    else echo "$sp: EDTA rerun FAILED"; tail -n 8 $E/logs/edta_rerun_$sp.log; fi
done
echo "[$(date)] re-mask edta2 for:$OK"
for sp in $OK; do
    [[ -d $E/racy_edta2/$sp ]] || { mkdir -p $E/racy_edta2; cp -a $E/$sp/edta2 $E/racy_edta2/$sp; }
    cp $B/work/edta_fixed/$sp/$sp.edta_library.fa $E/$sp/edta2/library.fa
    SLURM_ARRAY_TASK_ID=${ROW[$sp]} bash $E/rm_one.sh > $E/logs/rm_${ROW[$sp]}_edta_rerun.out 2>&1 &
done
wait
for sp in $OK; do tail -n 1 $E/logs/rm_${ROW[$sp]}_edta_rerun.out; done
echo "[$(date)] lib-test + family level"
[[ -n $OK ]] && bash $E/libtest_species.sh $OK
[[ -n $OK ]] && bash $E/family_level.sh $OK
echo "[$(date)] all done"
