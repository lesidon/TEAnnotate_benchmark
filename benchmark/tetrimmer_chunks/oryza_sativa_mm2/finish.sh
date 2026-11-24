#!/bin/bash
# finish.sh - after every chunk: merge, resume the rice MMseqs2 validation run (with TEtrimmer;
# the TETRIMMER task finishes the merged output with --continue_analysis), score both arms.
set -eo pipefail
/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/tetrimmer_chunks/tt_chunks.py merge /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ablation/reclass_mm2/oryza_sativa/edta+heliano+hite/work/76/e027dcb0af4f6b89b9f31ecfce5775 /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/tetrimmer_chunks/oryza_sativa_mm2
bash /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/tetrimmer_chunks/oryza_sativa_mm2/run_with_tt.sh with
bash /trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ablation/score_mm2_rice.sh wo with
