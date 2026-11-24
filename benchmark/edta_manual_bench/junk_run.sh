#!/bin/bash
# Cross-tool support + intrinsic features for the teannotate w/o-TEtrimmer libraries (v4).
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
W=$E/junk_analysis; mkdir -p $W; cd $W
IMG=$B/work/singularity_cache/quay.io-biocontainers-edta-2.3.0--hdfd78af_0.img
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
SX="singularity exec --no-home --bind $B --bind $(readlink -f $B) $IMG"
for sp in arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola; do
  Q=$E/$sp/teannotate_wo_ref4/library.fa
  # the standalone tool libraries that fed this route (EDTA: the version the v4 run used)
  TOOLS="hite:$E/$sp/hite/library.fa earlgrey:$E/$sp/earlgrey/library.fa heliano:$E/$sp/heliano/library.fa edta:$B/work/edta_fixed_racy/$sp/$sp.edta_library.fa"
  : > $sp.hits.tsv
  for t in $TOOLS; do
    tool=${t%%:*}; L=${t#*:}
    [[ -s $L ]] || continue
    $SX blastn -query $Q -subject $L -outfmt '6 qseqid sseqid pident length qstart qend qlen' -evalue 1e-10 -perc_identity 80 -dust no 2>/dev/null < /dev/null \
      | awk -v t=$tool 'BEGIN{OFS="\t"} {print t, $0}' >> $sp.hits.tsv
  done
  # support = number of OTHER tools whose library covers >=50% of the consensus at >=80% identity
  $PY - $sp.hits.tsv $sp.support.tsv <<'PY'
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
  $PY $E/junk_features.py $sp teannotate_wo_ref4 $sp.support.tsv $sp.features.tsv
done
head -3 arabidopsis_thaliana.features.tsv

