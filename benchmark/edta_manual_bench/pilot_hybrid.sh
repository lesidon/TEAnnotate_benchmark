#!/bin/bash
# pilot_hybrid.sh - TEtrimmer-arm variants (v4), re-masked and scored like the benchmark rows:
#   tt_evid   : TEtrimmer-arm library, keep only consensi with evidence (domain | cross-tool |
#               >=2 full copies) and drop non-TE labels  (features computed on this library)
#   tt_hybrid : tt_evid + inputs TEtrimmer dropped (skipped, not written back) that carry a
#               TE domain or cross-tool support, taken from the w/o-TEtrimmer library
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; P=$E/pilot; W=$E/junk_analysis; mkdir -p $P
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
EIMG=$B/work/singularity_cache/quay.io-biocontainers-edta-2.3.0--hdfd78af_0.img
RIMG=$B/work/pipeline/arabidopsis_thaliana/work/singularity/dfam-tetools-1.89.2.img
SX="singularity exec --no-home --bind $B --bind $(readlink -f $B) $EIMG"
SPECIES="arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola"
for sp in $SPECIES; do   # features + domains for the TEtrimmer-arm library
  Q=$E/$sp/teannotate_ref4/library.fa
  TOOLS="hite:$E/$sp/hite/library.fa earlgrey:$E/$sp/earlgrey/library.fa heliano:$E/$sp/heliano/library.fa edta:$B/work/edta_fixed_racy/$sp/$sp.edta_library.fa"
  : > $W/$sp.ref4.hits.tsv
  for t in $TOOLS; do
    tool=${t%%:*}; L=${t#*:}; [[ -s $L ]] || continue
    $SX blastn -query $Q -subject $L -outfmt '6 qseqid sseqid pident length qstart qend qlen' -evalue 1e-10 -perc_identity 80 -dust no 2>/dev/null < /dev/null \
      | awk -v t=$tool 'BEGIN{OFS="\t"} {print t, $0}' >> $W/$sp.ref4.hits.tsv
  done
  sed -n '/^  \$PY - \$sp.hits.tsv/,/^PY$/p' $E/junk_run.sh | sed '1d;$d' > /tmp/support_$$.py   # reuse the support code
  $PY /tmp/support_$$.py $W/$sp.ref4.hits.tsv $W/$sp.ref4.support.tsv 2>/dev/null || \
    $PY - $W/$sp.ref4.hits.tsv $W/$sp.ref4.support.tsv <<'PY'
import sys, re, collections
def source(name):
    n = name.split('#')[0].replace('_', '-')
    if re.match(r'TE-\d+', n): return 'edta'
    if re.match(r'(LTR|TIR|Helitron|Homology-Non-LTR|Denovo-Non-LTR)-\d+', n): return 'hite'
    if n.startswith('insertion-HLE'): return 'heliano'
    if n.startswith('rnd-'): return 'earlgrey'
    return 'other'
iv = collections.defaultdict(list); qlen = {}
for l in open(sys.argv[1]):
    tool, q, s, pid, ln, qs, qe, ql = l.split()
    if tool == source(q): continue
    iv[(q.split('#')[0], tool)].append((int(qs), int(qe))); qlen[q.split('#')[0]] = int(ql)
sup = collections.Counter()
for (q, tool), v in iv.items():
    v.sort(); cov = 0; s, e = v[0]
    for a, b in v[1:]:
        if a > e + 1: cov += e - s + 1; s, e = a, b
        else: e = max(e, b)
    cov += e - s + 1
    if cov >= 0.5 * qlen[q]: sup[q] += 1
with open(sys.argv[2], 'w') as o:
    for q, n in sup.items(): o.write(f'{q}\t{n}\n')
PY
  rm -f /tmp/support_$$.py
  $PY $E/junk_features.py $sp teannotate_ref4 $W/$sp.ref4.support.tsv $W/$sp.ref4.features.tsv
  mkdir -p $W/tesorter/$sp.ref4 && cd $W/tesorter/$sp.ref4 && cp $Q lib.fa
  $SX TEsorter lib.fa -db rexdb -p 8 -pre $sp -tmp ./tmp > tesorter.log 2>&1 < /dev/null
  cd $E
done
$PY - $B $E $P $SPECIES <<'PY'
import sys, os, csv, glob
B, E, P, species = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
NONTE = ("Satellite", "rRNA", "tRNA", "snRNA", "Simple", "Low_complexity", "Cent")
def tiers(feat_tsv, cls_tsv):
    dom = set(l.split("\t")[0].split("#")[0] for l in open(cls_tsv) if not l.startswith("#"))
    feats = {r["name"]: r for r in csv.DictReader(open(feat_tsv), delimiter="\t")}
    def t(name, label):
        if label.startswith(NONTE): return 0
        if name in dom: return 1
        f = feats.get(name)
        if f is None: return 4
        if int(f["n_support"]) >= 1: return 2
        if int(f["n_full"]) >= 2: return 3
        return 4
    return t
def read(path):
    seqs, cur = [], None
    for l in open(path):
        if l.startswith(">"):
            name, _, label = l[1:].split()[0].partition("#"); cur = [name, label, l, []]; seqs.append(cur)
        else: cur[3].append(l)
    return seqs
for sp in species:
    W = f"{E}/junk_analysis"
    t4 = tiers(f"{W}/{sp}.ref4.features.tsv", f"{W}/tesorter/{sp}.ref4/{sp}.cls.tsv")
    tw = tiers(f"{W}/{sp}.features.tsv", f"{W}/tesorter/{sp}/{sp}.cls.tsv")
    ref4 = read(f"{E}/{sp}/teannotate_ref4/library.fa")
    wo = read(f"{E}/{sp}/teannotate_wo_ref4/library.fa")
    # inputs TEtrimmer dropped
    dropped = set()
    for cmd in glob.glob(f"{B}/work/refined_dedup_v4/{sp}/work/*/*/.command.sh"):
        txt = open(cmd).read()
        if "TETRIMMER_CMD" in txt and "--input_file" in txt:
            s = glob.glob(cmd.replace(".command.sh", "*_tetrimmer/summary.txt"))
            if s:
                od = s[0].rsplit("/", 1)[0]
                written = set(l[1:].split("#")[0].strip() for l in open(f"{od}/TEtrimmer_consensus.fasta") if l.startswith(">"))
                dropped = {r["input_name"] for r in csv.DictReader(open(s[0])) if r["status"] == "skipped" and r["input_name"] not in written}
                break
    evid = [x for x in ref4 if t4(x[0], x[1]) in (1, 2, 3)]
    rescue = [x for x in wo if x[0] in dropped and tw(x[0], x[1]) in (1, 2)]
    for var, seqs in (("tt_evid", evid), ("tt_hybrid", evid + rescue)):
        os.makedirs(f"{P}/{sp}/{var}", exist_ok=True)
        with open(f"{P}/{sp}/{var}/library.fa", "w") as o:
            for name, label, h, body in seqs: o.write(h); o.writelines(body)
        print(f"{sp}/{var}: {len(seqs)} seqs (TEtrimmer arm {len(ref4)}; rescued {len(rescue) if var == 'tt_hybrid' else 0})")
PY
mask() {
    local sp=$1 var=$2
    local KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$sp" '$2==s{print $1}')
    local G=$B/genomes/$KING/$sp/$sp.fna OUT=$P/$sp/$var
    singularity exec --no-home --writable-tmpfs --bind "$B" "$RIMG" bash -c "
        export HOME='$OUT'; cd '$OUT'
        rm -rf rmlib; mkdir rmlib; for f in /opt/RepeatMasker/Libraries/*; do ln -s \$f rmlib/; done
        RepeatMasker -e ncbi -pa 4 -q -no_is -norna -nolow -div 40 -libdir '$OUT/rmlib' \
            -lib '$OUT/library.fa' -cutoff 225 -dir '$OUT' '$G'" > "$OUT/rm.log" 2>&1
    ( cd $OUT && $PY $B/tools/EDTA/lib-test.py --genome "$G" --reference "$E/$sp/reference/$sp.fna.out" \
          --test "$OUT/$sp.fna.out" --extended_report ) > $OUT/libtest.log 2>&1
    echo "[$(date)] $sp/$var: $(tail -n +4 $OUT/$sp.fna.out | wc -l) hits"
}
for sp in $SPECIES; do for var in tt_evid tt_hybrid; do mask $sp $var & done; done
wait
$PY - $P $SPECIES <<'PY'
import sys, re, glob
P, species = sys.argv[1], sys.argv[2:]
with open(f"{P}/pilot_results_tt.csv", "w") as o:
    o.write("species,variant,TE type,sens,spec,accu,prec,FDR,F1,TP,TN,FP,FN\n")
    for sp in species:
        for var in ("tt_evid", "tt_hybrid"):
            rep = f"{P}/{sp}/{var}/{sp}.fna.out.repeatmasker.lib.report"
            for line in open(rep):
                f = line.rstrip("\n").split("\t")
                m = re.search(r"\.(\w+)\.lib\.report$", f[0])
                if m and len(f) >= 11: o.write(",".join([sp, var, m.group(1)] + f[1:11]) + "\n")
print("done")
PY

