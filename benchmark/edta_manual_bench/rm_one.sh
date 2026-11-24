#!/bin/bash
# One RepeatMasker run, exactly as in the EDTA manual (library-comparison section).
# -libdir -> task-local symlinks of RepeatMasker's Libraries: on nodes whose
# --writable-tmpfs is fuse-overlayfs, makeblastdb cannot build general.working in /opt.
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
read -r SP M LIB < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$E/specs.tsv")
KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
G=$B/genomes/$KING/$SP/$SP.fna
OUT=$E/$SP/$M
echo "[$(date)] $SP/$M on $(hostname): lib=$OUT/library.fa genome=$G"
for attempt in 1 2; do
    rm -rf "$OUT"/$SP.fna.{out,tbl,masked,cat.gz} "$OUT/rmlib"
    singularity exec --no-home --writable-tmpfs --bind "$B" "$IMG" bash -c "
        export HOME='$OUT'; cd '$OUT'
        mkdir rmlib; for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        RepeatMasker -e ncbi -pa 4 -q -no_is -norna -nolow -div 40 -libdir '$OUT/rmlib' \
            -lib '$OUT/library.fa' -cutoff 225 -dir '$OUT' '$G'" > "$OUT/rm.log" 2>&1 || true
    if [[ -s $OUT/$SP.fna.out ]] && [[ $(tail -n +4 "$OUT/$SP.fna.out" | wc -l) -gt 0 ]]; then
        echo "[$(date)] OK attempt $attempt: $(tail -n +4 $OUT/$SP.fna.out | wc -l) hits"; exit 0
    fi
    echo "[$(date)] attempt $attempt produced no .out; RepeatMasker log tail:"; tail -5 "$OUT/rm.log"
done
echo "FAILED: $SP/$M"; exit 1
