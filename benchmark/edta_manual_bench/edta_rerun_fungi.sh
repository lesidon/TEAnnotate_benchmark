#!/bin/bash
# edta_rerun_fungi.sh - fungal EDTA baselines without the rice fallback: EDTA found no
# SINE (Zymoseptoria) / no SINE+LINE (Septoria) and --force 1 copied rice7.0.0.liban.SINE/LINE
# into the raw set. Empty those raw files (EDTA allows empty SINE/LINE, as with --force 0;
# originals kept in work/edta_fixed_racy/<sp>/raw_rice/), redo filter+final with the
# race-free module script (31e2e6f), then re-mask, lib-test and family level.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
IMG=$B/work/singularity_cache/quay.io-biocontainers-edta-2.3.0--hdfd78af_0.img
declare -A ROW=([zymoseptoria_tritici]=40 [septoria_linicola]=43)
DB=$(mktemp -d); singularity exec --no-home $IMG bash -c "cp /usr/local/share/EDTA/database/rice7.0.0.liban.* $DB/" < /dev/null
declare -A PID
for sp in zymoseptoria_tritici septoria_linicola; do
    RAWD=$B/work/edta_fixed/$sp/$sp.fna.mod.EDTA.raw
    mkdir -p $B/work/edta_fixed_racy/$sp/raw_rice
    for t in LTR:LTR.raw.fa LINE:LINE.raw.fa SINE:SINE.raw.fa; do
        f=$RAWD/$sp.fna.mod.${t#*:}
        if [[ -s $f ]] && cmp -s $f $DB/rice7.0.0.liban.${t%%:*}; then
            [[ ${t%%:*} == LTR ]] && { echo "$sp: rice LTR - not handled"; exit 1; }
            cp $f $B/work/edta_fixed_racy/$sp/raw_rice/ && : > $f
            echo "$sp: emptied rice ${t%%:*} raw file"
        fi
    done
    CPUS=8 bash $B/rerun_edta_filter.sh $sp > $E/logs/edta_rerun_$sp.log 2>&1 &
    PID[$sp]=$!
done
rm -rf $DB
OK=""
for sp in zymoseptoria_tritici septoria_linicola; do
    if wait ${PID[$sp]}; then OK="$OK $sp"; tail -n 2 $E/logs/edta_rerun_$sp.log
    else echo "$sp: EDTA rerun FAILED"; tail -n 8 $E/logs/edta_rerun_$sp.log; tail -n 20 $B/work/edta_fixed/$sp/edta_filter.log; fi
done
for sp in $OK; do
    L=$B/work/edta_fixed/$sp/$sp.edta_library.fa
    echo "$sp: $(grep -c '>' $L) seqs; SINE $(grep -c '#SINE' $L); LINE $(grep -c '#LINE' $L); Os00 names $(grep -c '>Os00' $L)"
    [[ -d $E/racy_edta2/$sp ]] || { mkdir -p $E/racy_edta2; cp -a $E/$sp/edta2 $E/racy_edta2/$sp; }
    cp $L $E/$sp/edta2/library.fa
    SLURM_ARRAY_TASK_ID=${ROW[$sp]} bash $E/rm_one.sh > $E/logs/rm_${ROW[$sp]}_edta_rerun.out 2>&1 &
done
wait
for sp in $OK; do tail -n 1 $E/logs/rm_${ROW[$sp]}_edta_rerun.out; done
[[ -n $OK ]] && bash $E/libtest_species.sh $OK
[[ -n $OK ]] && bash $E/family_level.sh $OK
echo "[$(date)] all done"
