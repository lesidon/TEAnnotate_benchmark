#!/bin/bash
# (1) release the held-out vertebrate pipeline jobs; replace their standalone-EDTA jobs by the
#     EDTA seq-ID fix (all three genomes have 14-character NW_ IDs);
# (2) score the four re-classification variants on the 8 genomes with replay data.
set -e
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
HO=$B/heldout; V=$B/vote_v10; PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
BAD=$(cat $HO/exclude_nodes.txt)
scontrol release 8482970 8482972 8482974 && echo "released pipeline jobs 8482970 (gga), 8482972 (tgu), 8482974 (mmu)"
scancel 8482971 8482973 8482975 && echo "cancelled old-script EDTA jobs 8482971 8482973 8482975"
for sp in gallus_gallus:gga taeniopygia_guttata:tgu mus_musculus:mmu; do
    j=$(sbatch --parsable -p mem -c 32 --mem=250G -t 6-00:00:00 --exclude=$BAD -J ho_edta9e_${sp#*:} \
        -o $B/logs/ho_edta9e_${sp#*:}_%j.out $HO/run_edta_v9e.sh ${sp%%:*})
    echo "EDTA (fixed) ${sp%%:*}: $j"
done
cd $V
$PY -m py_compile variants_eval.py bin/ensemble_vote.py && echo "syntax ok"
j=$(sbatch --parsable -p htc -c 16 --mem=64G -t 08:00:00 -J v10variants -o $V/logs/variants_%j.out \
    --wrap "cd $V && $PY variants_eval.py arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola caenorhabditis_briggsae caenorhabditis_elegans blumeria_hordei ectocarpus_sp7")
echo "variant scoring: $j"
squeue -u leonid.sidorov -o '%i %P %j %T %R' | grep -v JobHeldUser

