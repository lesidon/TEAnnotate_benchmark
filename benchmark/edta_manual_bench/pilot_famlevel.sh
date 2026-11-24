#!/bin/bash
# Family-level (RepeatModeler2-style) scores for the pilot variants, same procedure as family_level.sh.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; P=$E/pilot
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
run() {
    local sp=$1 var=$2
    local REF=$E/$sp/reference/library.fa LIB=$P/$sp/$var/library.fa OUT=$P/$sp/$var/famlevel
    rm -rf $OUT; mkdir -p $OUT
    singularity exec --no-home --writable-tmpfs --bind $B --bind $(readlink -f $B) $IMG bash -c "
        cd $OUT && mkdir rmlib && for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        awk '/^>/{print \">q\" ++i; next} 1' $LIB > library.fa
        RepeatMasker -e ncbi -pa 2 -nolow -no_is -libdir $OUT/rmlib -lib $REF -dir $OUT library.fa" > $OUT/rm.log 2>&1
}
for sp in arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola; do
  for var in evid_lite evid_strict tt_evid tt_hybrid; do run $sp $var & done
done
wait
for sp in arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola; do
  $PY $E/family_level.py $sp $P/$sp $P/${sp}_family_level.csv > /dev/null 2>&1
  awk -F, 'NR==1 || $2=="all"' $P/${sp}_family_level.csv
done

