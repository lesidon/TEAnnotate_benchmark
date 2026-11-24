#!/bin/bash
# test: run_hite_v9h.sh on a 6 Mb piece of C. elegans; scratch must land in hite_work, not /tmp
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
T=$B/tmp/hite_workdir_test; rm -rf $T $WORK_DIR/hite_v9h/celtest; mkdir -p $T/orig
ORIG=$(dirname $(ls $WORK_DIR/heldout_v9/caenorhabditis_elegans/work/54/199962b04780b54909ac49d6805f*/.command.sh | head -n 1))
( set +o pipefail; awk '/^>/{n++} n==1' $(readlink -f $ORIG/caenorhabditis_elegans.fna) | head -c 6000000 ) > $T/orig/celtest.fna
sed 's/caenorhabditis_elegans/celtest/g' $ORIG/.command.sh > $T/orig/.command.sh
before=$(ls -d /tmp/coarse_boundary_* 2>/dev/null | wc -l || true)
LIB_DIR=$T/lib bash $B/heldout/run_hite_v9h.sh celtest $T/orig
D=$WORK_DIR/hite_v9h/celtest
echo "tmp dirs before/after: $before / $(ls -d /tmp/coarse_boundary_* 2>/dev/null | wc -l || true)"
grep -a 'Temporary work dir' $D/hite.log
grep -ac 'coarse_boundary.py' $D/hite.log
grep -aoE -- '-w [^ ]+' $D/hite.log | sort -u | head -n 3
ls $D; grep -c '>' $D/celtest.hite_library.fa

