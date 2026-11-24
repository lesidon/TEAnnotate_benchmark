#!/bin/bash
# Prototype of the MMseqs2 dedup path on finished CD_HIT_EST inputs (rice, Danio) and on
# maize (timing only). Writes CD-HIT-format files and compares with CD-HIT's clustering.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
E=$B/heldout/mmseqs_env/bin
H=$B/heldout
OUT=$H/mmseqs_test; mkdir -p $OUT
RICE=$B/ablation/reclass/oryza_sativa/edta+heliano+hite/work/13/c63f41e16c217cae250616eb4d608a
DANIO=$WORK_DIR/heldout_v10/danio_rerio/work/56/db347b1f235e1364f3f149c8c47aaa
MAIZE=$WORK_DIR/heldout_v10/zea_mays/work/f7/4aee90f4b33b38618345898face9e3

dedup() {   # fasta outdir threads cluster-extra-args...
    local fa=$1 d=$2 t=$3; shift 3
    mkdir -p $d; cd $d; rm -rf db* clu* rep* res* tmp*
    local t0=$(date +%s)
    $E/mmseqs createdb $fa db --dbtype 2 -v 1
    $E/mmseqs cluster db clu tmp --min-seq-id 0.80 -c 0.80 --cov-mode 1 --cluster-mode 2 \
        --threads $t -v 1 "$@"
    $E/mmseqs createtsv db db clu clusters.tsv -v 1
    local t1=$(date +%s)
    # member coordinates on its representative: strand-aware nucleotide search of every
    # sequence against the representatives, kept only for the member's own representative
    $E/mmseqs createsubdb clu db rep -v 1
    $E/mmseqs search db rep res tmp --search-type 3 --strand 2 -a --max-seqs 1000 -e 10 \
        --threads $t -v 1
    $E/mmseqs convertalis db rep res aln.tsv --format-output query,target,fident,qstart,qend,tstart,tend -v 1
    local t2=$(date +%s)
    /usr/bin/python3 $H/mmseqs_to_clstr.py --fasta $fa --clusters clusters.tsv --alignments aln.tsv --out dedup.fa
    echo "  time: cluster $((t1 - t0)) s, coordinates $((t2 - t1)) s"
    rm -rf tmp
}

for sp in danio rice; do
    case $sp in danio) T=$DANIO ;; rice) T=$RICE ;; esac
    fa=$(ls $T/*.combined_library.fa)
    for v in default sens; do
        echo "== $sp $v"
        case $v in
            default) dedup $fa $OUT/$sp/$v 16 ;;
            sens)    dedup $fa $OUT/$sp/$v 16 -s 7.5 --max-seqs 1000 ;;
        esac
        /usr/bin/python3 $H/clstr_compare.py $T/*.dedup_library.fa.clstr $OUT/$sp/$v/dedup.fa.clstr
    done
done
echo "== maize sens (timing)"
dedup $MAIZE/zea_mays.combined_library.fa $OUT/maize/sens 16 -s 7.5 --max-seqs 1000
grep -c '>' $OUT/maize/sens/dedup.fa
echo "[$(date)] done"
