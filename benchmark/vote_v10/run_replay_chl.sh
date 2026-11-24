#!/bin/bash
set -eo pipefail
cd /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/vote_v10
TAG=_ho2 /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python prep_unknown.py chlamydomonas_reinhardtii
bash run_cls.sh plant _ho2
TAG=_ho2 /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python apply_vote.py chlamydomonas_reinhardtii
TAG=_ho2 /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python bench_eval.py chlamydomonas_reinhardtii
/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python diag_relabel.py chlamydomonas_reinhardtii
/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python agree_check.py chlamydomonas_reinhardtii
TAG=_ho2 /beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python variants_eval.py chlamydomonas_reinhardtii
