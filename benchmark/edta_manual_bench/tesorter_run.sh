#!/bin/bash
# TEsorter (REXdb protein domains) on the four teannotate w/o-TEtrimmer libraries (v4).
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
W=$E/junk_analysis/tesorter; mkdir -p $W; cd $W
IMG=$B/work/singularity_cache/quay.io-biocontainers-edta-2.3.0--hdfd78af_0.img
for sp in arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola; do
  mkdir -p $sp; cd $sp
  # TEsorter keys on the ID before '#'; keep IDs as in the library
  cp $E/$sp/teannotate_wo_ref4/library.fa lib.fa
  singularity exec --no-home --bind $B --bind $(readlink -f $B) $IMG \
      TEsorter lib.fa -db rexdb -p ${SLURM_CPUS_PER_TASK:-8} -pre $sp -tmp ./tmp > tesorter.log 2>&1 < /dev/null
  echo "$sp: $(ls) ; classified $(grep -vc '^#' $sp.cls.tsv 2>/dev/null)"
  cd $W
done

