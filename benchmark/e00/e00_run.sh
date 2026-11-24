#!/bin/bash
# e00_run.sh famlevel <sp> | genome <sp> <library>
# The RepeatMasker runs of E00 with the Table 2 settings: famlevel = every stage library
# masked with the curated reference (family_level.sh); genome = the genome masked with one
# E00 library (rm_one.sh). -libdir: task-local RepeatMasker libraries (fuse-overlay nodes).
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; X=$B/e00
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
CPUS=${SLURM_CPUS_PER_TASK:-4}; PA=$(( CPUS / 4 > 0 ? CPUS / 4 : 1 ))
MODE=$1; SP=$2

rm_in() {  # <dir> <library> <target> <RepeatMasker options...>
    local D=$1 LIB=$2 T=$3; shift 3
    singularity exec --no-home --writable-tmpfs --bind "$B" --bind "$(readlink -f $B)" "$IMG" bash -c "
        export HOME='$D'; cd '$D'
        rm -rf rmlib; mkdir rmlib; for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        RepeatMasker -e ncbi -pa $PA $* -libdir '$D/rmlib' -lib '$LIB' -dir '$D' '$T'" > "$D/rm.log" 2>&1 || true
}

if [[ $MODE == famlevel ]]; then
    REF=$E/$SP/reference/library.fa
    for FA in $X/$SP/stages/*.fa; do
        S=$(basename $FA .fa); O=$X/$SP/famlevel/$S
        [[ -s $O/library.fa.out && $O/library.fa.out -nt $FA ]] && { echo "$SP/$S: up to date"; continue; }
        rm -rf $O; mkdir -p $O; cp $FA $O/library.fa
        rm_in $O $REF library.fa -nolow -no_is
        [[ -s $O/library.fa.out ]] || { echo "FAILED: $SP/$S famlevel"; tail -5 $O/rm.log; exit 1; }
        echo "[$(date)] $SP/$S: $(tail -n +4 $O/library.fa.out | wc -l) hits"
    done
elif [[ $MODE == genome ]]; then
    L=$3; O=$X/$SP/genome_rm/$L
    KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
    G=$B/genomes/$KING/$SP/$SP.fna
    echo "[$(date)] $SP/$L on $(hostname), $(grep -c '>' $X/$SP/genome_libs/$L.fa) seqs, -pa $PA"
    for attempt in 1 2; do
        rm -rf $O; mkdir -p $O; cp $X/$SP/genome_libs/$L.fa $O/library.fa
        rm_in $O $O/library.fa $G -q -no_is -norna -nolow -div 40 -cutoff 225
        if [[ $(tail -n +4 $O/$SP.fna.out 2>/dev/null | wc -l) -gt 0 ]]; then
            echo "[$(date)] OK attempt $attempt: $(tail -n +4 $O/$SP.fna.out | wc -l) hits"; exit 0
        fi
        echo "[$(date)] attempt $attempt produced no .out:"; tail -5 $O/rm.log
    done
    echo "FAILED: $SP/$L"; exit 1
fi
