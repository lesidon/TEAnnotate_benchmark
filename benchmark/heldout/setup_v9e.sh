#!/bin/bash
# v9e = frozen v9 + the EDTA sequence-ID fix (16d03c3) only. Reruns: Chlamydomonas (full
# pipeline, plant route) and Blumeria's standalone EDTA baseline; then bench + v10-vote replay.
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
H9=$B/pipeline_head_v9; H=$B/pipeline_head_v9e; HO=$B/heldout; V=$B/vote_v10
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python

# 1. snapshot
if [[ ! -d $H ]]; then
    cp -a $H9 $H
    cp $HO/v9e_edta_main.nf $H/modules/local/edta/main.nf
    echo "3ffb6ae + modules/local/edta/main.nf from 16d03c3 (EDTA seq-ID fix)" > $H/COMMIT
fi
echo "== snapshot differs from v9 in:"; diff -rq $H9 $H | grep -v '/work\|\.nextflow' || true

# 2. runners
sed -e "s#^H=\$B/pipeline_head_v9\$#H=\$B/pipeline_head_v9e#" \
    -e "s#^R=\$WORK_DIR/heldout_v9/\$SP#R=\$WORK_DIR/heldout_v9e/\$SP#" \
    -e 's#^for t in $(find $R/out_with $R/out_wo -name#set +o pipefail\nfor t in $(find $R/out_with $R/out_wo -name#' \
    $HO/run_heldout_v9b.sh > $HO/run_heldout_v9e.sh
sed -e "s#edta_module_31e2e6f.sh#edta_module_16d03c3.sh#" -e "s#^D=\$WORK_DIR/edta_fixed/\$SP#D=\$WORK_DIR/edta_v9e/\$SP#" \
    -e "s#(v9 module script)#(v9e module script, EDTA seq-ID fix)#" \
    $HO/run_edta_v9.sh > $HO/run_edta_v9e.sh
chmod +x $HO/run_heldout_v9e.sh $HO/run_edta_v9e.sh
bash -n $HO/run_heldout_v9e.sh && bash -n $HO/run_edta_v9e.sh && bash -n $B/edta_module_16d03c3.sh && echo "syntax ok"
echo "== runner changes"; diff $HO/run_heldout_v9b.sh $HO/run_heldout_v9e.sh || true; diff $HO/run_edta_v9.sh $HO/run_edta_v9e.sh || true
grep -c 'edtaid' $H/modules/local/edta/main.nf

# 3. jobs (kernel-5.4.265 nodes only: ct02-04)
BAD=$(cat $HO/exclude_nodes.txt)
j1=$(sbatch --parsable -p mem -c 32 --mem=160G -t 5-00:00:00 --exclude=$BAD -J ho_v9e_chl \
     -o $B/logs/ho_v9e_chl_%j.out $HO/run_heldout_v9e.sh chlamydomonas_reinhardtii)
j2=$(sbatch --parsable -p mem -c 32 --mem=120G -t 3-00:00:00 --exclude=$BAD -J ho_edta9e_bho \
     -o $B/logs/ho_edta9e_bho_%j.out $HO/run_edta_v9e.sh blumeria_hordei)
j3=$(sbatch --parsable -p htc -c 16 --mem=64G -t 1-00:00:00 --dependency=afterok:$j1:$j2 -J ho_bench_v9e \
     -o $B/logs/ho_bench_v9e_%j.out $B/edta_manual_bench/heldout_bench_v9.sh chlamydomonas_reinhardtii blumeria_hordei)
echo "chlamydomonas v9e pipeline: $j1; blumeria EDTA (fixed): $j2; bench: $j3"
cat > $V/run_replay_chl.sh <<EOF
#!/bin/bash
set -eo pipefail
cd $V
TAG=_ho2 $PY prep_unknown.py chlamydomonas_reinhardtii
bash run_cls.sh plant _ho2
TAG=_ho2 $PY apply_vote.py chlamydomonas_reinhardtii
TAG=_ho2 $PY bench_eval.py chlamydomonas_reinhardtii
$PY diag_relabel.py chlamydomonas_reinhardtii
$PY agree_check.py chlamydomonas_reinhardtii
EOF
j4=$(sbatch --parsable -p htc -c 16 --mem=60G -t 12:00:00 --dependency=afterok:$j3 -J v10replay_chl \
     -o $V/logs/replay_chl_%j.out $V/run_replay_chl.sh)
echo "chlamydomonas new-vote replay: $j4"
squeue -u leonid.sidorov -o '%i %P %j %T %R' | grep -E 'v9e|edta9e|replay'

