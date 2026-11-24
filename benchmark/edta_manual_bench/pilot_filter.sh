#!/bin/bash
# pilot_filter.sh - evidence-filtered variants of the v4 w/o-TEtrimmer libraries, re-masked and
# scored with lib-test exactly like the benchmark rows (separate pilot/ tree; benchmark untouched).
#   evid_lite   : drop non-TE labels (Satellite/rRNA/tRNA/...) and consensi with no evidence
#   evid_strict : additionally drop consensi whose only evidence is >=2 full-length copies
# evidence = TE protein domain (TEsorter/REXdb) | same element in another tool's library
#            (>=50% at >=80% id) | >=2 full-length genomic copies
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; P=$E/pilot; mkdir -p $P
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
IMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
SPECIES="arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola"
$PY - $E $P $SPECIES <<'PY'
import sys, csv
E, P, species = sys.argv[1], sys.argv[2], sys.argv[3:]
NONTE = ("Satellite", "rRNA", "tRNA", "snRNA", "Simple", "Low_complexity", "Cent")
import os
for sp in species:
    dom = set(l.split("\t")[0].split("#")[0] for l in open(f"{E}/junk_analysis/tesorter/{sp}/{sp}.cls.tsv") if not l.startswith("#"))
    feats = {r["name"]: r for r in csv.DictReader(open(f"{E}/junk_analysis/{sp}.features.tsv"), delimiter="\t")}
    def tier(name, label):
        if label.startswith(NONTE): return 0
        if name in dom: return 1
        f = feats.get(name)
        if f is None: return 4          # no genomic hit at all
        if int(f["n_support"]) >= 1: return 2
        if int(f["n_full"]) >= 2: return 3
        return 4
    seqs, cur = [], None
    for l in open(f"{E}/{sp}/teannotate_wo_ref4/library.fa"):
        if l.startswith(">"):
            name, _, label = l[1:].split()[0].partition("#")
            cur = [name, label, l, []]; seqs.append(cur)
        else:
            cur[3].append(l)
    for var, keep in (("evid_lite", {1, 2, 3}), ("evid_strict", {1, 2})):
        os.makedirs(f"{P}/{sp}/{var}", exist_ok=True)
        n = 0
        with open(f"{P}/{sp}/{var}/library.fa", "w") as o:
            for name, label, h, body in seqs:
                if tier(name, label) in keep:
                    o.write(h); o.writelines(body); n += 1
        print(f"{sp}/{var}: kept {n} of {len(seqs)}")
PY
mask() {   # same RepeatMasker call as rm_one.sh
    local sp=$1 var=$2
    local KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$sp" '$2==s{print $1}')
    local G=$B/genomes/$KING/$sp/$sp.fna OUT=$P/$sp/$var
    singularity exec --no-home --writable-tmpfs --bind "$B" "$IMG" bash -c "
        export HOME='$OUT'; cd '$OUT'
        rm -rf rmlib; mkdir rmlib; for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        RepeatMasker -e ncbi -pa 4 -q -no_is -norna -nolow -div 40 -libdir '$OUT/rmlib' \
            -lib '$OUT/library.fa' -cutoff 225 -dir '$OUT' '$G'" > "$OUT/rm.log" 2>&1
    ( cd $OUT && $PY $B/tools/EDTA/lib-test.py --genome "$G" --reference "$E/$sp/reference/$sp.fna.out" \
          --test "$OUT/$sp.fna.out" --extended_report ) > $OUT/libtest.log 2>&1
    echo "[$(date)] $sp/$var: $(tail -n +4 $OUT/$sp.fna.out | wc -l) hits"
}
for sp in $SPECIES; do for var in evid_lite evid_strict; do mask $sp $var & done; done
wait
$PY - $P $SPECIES <<'PY'
import sys, re, glob
P, species = sys.argv[1], sys.argv[2:]
with open(f"{P}/pilot_results.csv", "w") as o:
    o.write("species,variant,TE type,sens,spec,accu,prec,FDR,F1,TP,TN,FP,FN\n")
    for sp in species:
        for rep in sorted(glob.glob(f"{P}/{sp}/*/{sp}.fna.out.repeatmasker.lib.report")):
            var = rep.split("/")[-2]
            for line in open(rep):
                f = line.rstrip("\n").split("\t")
                m = re.search(r"\.(\w+)\.lib\.report$", f[0])
                if m and len(f) >= 11:
                    o.write(",".join([sp, var, m.group(1)] + f[1:11]) + "\n")
print(open(f"{P}/pilot_results.csv").read())
PY

