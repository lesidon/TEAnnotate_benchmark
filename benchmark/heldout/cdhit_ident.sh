#!/bin/bash
# Does a low -M (small serial first batch) change cd-hit-est's clustering?
# Re-runs two finished CD_HIT_EST tasks (rice raw 11,168 seqs, Danio 4,261) with a low -M
# and compares representatives and cluster membership with the pipeline's own output.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
OUT=$B/heldout/cdhit_ident; mkdir -p $OUT
RICE=$B/ablation/reclass/oryza_sativa/edta+heliano+hite/work/13/c63f41e16c217cae250616eb4d608a
DANIO=$WORK_DIR/heldout_v10/danio_rerio/work/56/db347b1f235e1364f3f149c8c47aaa

run_one() {   # name taskdir
    local name=$1 t=$2 env cdh m req
    env=$(grep -oE 'conda_cache/env-[^/ ]+' $t/.command.run | head -n 1)
    cdh=$WORK_DIR/$env/bin/cd-hit-est
    cd $OUT; mkdir -p $name; cd $name
    m=3000
    for try in 1 2 3; do
        /usr/bin/time -f "%e s wall, %M KB" stdbuf -o0 $cdh -i $t/*.combined_library.fa -o low.fa \
            -c 0.80 -aS 0.80 -G 0 -g 1 -b 500 -n 5 -d 0 -p 1 -T 16 -M $m > low.log 2>&1
        req=$(grep -oE 'greater than [0-9]+' low.log | grep -oE '[0-9]+')
        [ -z "$req" ] && break
        m=$((req + 800))
    done
    echo "[$name] env=$env -M $m first batch: $(grep -a 'comparing' low.log | head -n 1 | tr -s ' ') | $(tail -n 1 low.log)"
    # representatives
    grep '>' $t/*.dedup_library.fa | cut -c2- | awk '{print $1}' | sort > ref.ids
    grep '>' low.fa | cut -c2- | awk '{print $1}' | sort > low.ids
    echo "[$name] reps: pipeline $(wc -l < ref.ids), low-M $(wc -l < low.ids), differing $(comm -3 ref.ids low.ids | wc -l)"
    # membership: one line per cluster = sorted member IDs
    awk '/^>Cluster/{if(s)print s; s=""; next}{x=$3; sub(/^>/,"",x); sub(/\.\.\.$/,"",x); s=s" "x}END{print s}' $t/*.dedup_library.fa.clstr \
        | awk '{n=split($0,a," "); asort(a); o=a[1]; for(i=2;i<=n;i++)o=o" "a[i]; print o}' | sort > ref.clusters
    awk '/^>Cluster/{if(s)print s; s=""; next}{x=$3; sub(/^>/,"",x); sub(/\.\.\.$/,"",x); s=s" "x}END{print s}' low.fa.clstr \
        | awk '{n=split($0,a," "); asort(a); o=a[1]; for(i=2;i<=n;i++)o=o" "a[i]; print o}' | sort > low.clusters
    echo "[$name] clusters: pipeline $(wc -l < ref.clusters), low-M $(wc -l < low.clusters), differing $(comm -3 ref.clusters low.clusters | wc -l)"
    echo "[$name] clstr byte-identical: $(cmp -s $t/*.dedup_library.fa.clstr low.fa.clstr && echo yes || echo no)"
}

run_one danio $DANIO &
run_one rice  $RICE &
wait
echo "[$(date)] done"
