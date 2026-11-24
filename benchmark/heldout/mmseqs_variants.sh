#!/bin/bash
# MMseqs2 cluster settings vs CD-HIT on Danio and rice: merged members and pair agreement.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
E=$B/heldout/mmseqs_env/bin
H=$B/heldout
OUT=$H/mmseqs_test/variants; mkdir -p $OUT
RICE=$B/ablation/reclass/oryza_sativa/edta+heliano+hite/work/13/c63f41e16c217cae250616eb4d608a
DANIO=$WORK_DIR/heldout_v10/danio_rerio/work/56/db347b1f235e1364f3f149c8c47aaa

for sp in danio rice; do
    case $sp in danio) T=$DANIO ;; rice) T=$RICE ;; esac
    fa=$(ls $T/*.combined_library.fa)
    for v in k13 k11 k9 k11_reassign k11_alnmode3; do
        case $v in
            k13) args="-k 13" ;;
            k11) args="-k 11" ;;
            k9)  args="-k 9" ;;
            k11_reassign) args="-k 11 --cluster-reassign 1" ;;
            k11_alnmode3) args="-k 11 --alignment-mode 3" ;;
        esac
        d=$OUT/$sp/$v; rm -rf $d; mkdir -p $d; cd $d
        t0=$(date +%s)
        $E/mmseqs createdb $fa db --dbtype 2 -v 1
        $E/mmseqs cluster db clu tmp --min-seq-id 0.80 -c 0.80 --cov-mode 1 --cluster-mode 2 \
            --threads 16 -v 1 $args || { echo "== $sp $v FAILED"; continue; }
        $E/mmseqs createtsv db db clu clusters.tsv -v 1
        t1=$(date +%s)
        # pair comparison needs only membership: write a coordinate-free .clstr
        : > aln.tsv
        /usr/bin/python3 $H/mmseqs_to_clstr.py --fasta $fa --clusters clusters.tsv --alignments aln.tsv --out dedup.fa 2>/dev/null
        echo "== $sp $v ($((t1 - t0)) s): $(/usr/bin/python3 $H/clstr_compare.py $T/*.dedup_library.fa.clstr dedup.fa.clstr)"
        rm -rf tmp db* clu*
    done
done
echo "[$(date)] done"
