#!/bin/bash
# Offline test of the domain-label conflict rule on the v7 libraries (no re-masking):
# relabel hits in the existing RepeatMasker .out, then re-score with lib-test.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
W=$E/domainfix; mkdir -p $W
IMG=$B/work/singularity_cache/quay.io-biocontainers-tesorter-1.5.1--pyhdfd78af_0.img
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
cat > $W/relabel.py <<'PY'
import sys, csv, re
sys.path.insert(0, sys.argv[1])                         # pipeline bin (te_evidence.tesorter_label)
from te_evidence import tesorter_label
cls_tsv, out_in, rule, out_out, report = sys.argv[2:7]
def group(lab):
    l = lab.lower(); top = l.split('/')[0]
    if 'helitron' in l or top == 'rc': return 'helitron'
    if top in ('ltr', 'dirs', 'erv') or re.search('gypsy|copia|pao|retrovirus|caulimo', l): return 'ltr'
    if top in ('line', 'sine', 'ple', 'penelope', 'retroposon'): return 'nonltr'
    if top in ('dna', 'tir', 'mite'): return 'tir'
    return 'unknown'
dom = {}
for l in open(cls_tsv):
    if l.startswith('#') or not l.strip(): continue
    f = l.rstrip('\n').split('\t')
    dom[f[0].split('#')[0]] = (f[1], f[2])
new, n_conf, n_mix = {}, 0, 0
for name, (order, sf) in dom.items():
    new[name] = ('mixture', None) if order == 'mixture' else ('class', tesorter_label(order, sf))
relab = {}
changed = 0
with open(out_in) as fi, open(out_out, 'w') as fo:
    for i, line in enumerate(fi):
        f = line.split()
        if i < 3 or len(f) < 11:
            fo.write(line); continue
        name, lab = f[9], f[10]
        kind, nl = new.get(name, (None, None))
        target = None
        if kind == 'class' and nl and group(nl) != 'unknown' and group(nl) != group(lab):
            target = nl
        elif kind == 'mixture' and rule == 'R2' and group(lab) != 'unknown':
            target = 'Unknown'
        if target:
            relab[name] = (lab, target); f[10] = target; changed += 1
        fo.write(' '.join(f) + '\n')
with open(report, 'w') as r:
    for n, (a, b) in sorted(relab.items()): r.write(f'{n}\t{a}\t{b}\n')
print(f'{rule}: {len(relab)} consensi relabelled ({changed} hits)')
PY
run_sp() {
  sp=$1
  KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$sp" '$2==s{print $1}')
  G=$B/genomes/$KING/$sp/$sp.fna
  for m in teannotate_ref7 teannotate_wo_ref7; do
    D=$W/$sp/$m; mkdir -p $D; cd $D
    singularity exec --no-home --bind $B --bind $(readlink -f $B) $IMG \
        TEsorter $E/$sp/$m/library.fa -db rexdb -p 8 -pre ts -tmp ./tmp > ts.log 2>&1 < /dev/null
    for rule in R0 R1 R2; do
      mkdir -p $D/$rule; cd $D/$rule
      if [[ $rule == R0 ]]; then cp $E/$sp/$m/$sp.fna.out .; else
        $PY $W/relabel.py $B/pipeline_head_v7/bin $D/ts.cls.tsv $E/$sp/$m/$sp.fna.out $rule $sp.fna.out relabelled.tsv; fi
      $PY $B/tools/EDTA/lib-test.py --genome $G --reference $E/$sp/reference/$sp.fna.out --test $sp.fna.out > lt.log 2>&1
      mkdir -p h; cd h
      $PY -c "import sys; sys.path.insert(0,'$E'); import libtest_harmonised as h; h.harmonise('../$sp.fna.out', '$sp.fna.out')"
      $PY $B/tools/EDTA/lib-test.py --genome $G --reference $E/$sp/reference/harmonised/$sp.fna.out --test $sp.fna.out > lt.log 2>&1
      cd $D
    done
  done
}
for sp in arabidopsis_thaliana drosophila_melanogaster zymoseptoria_tritici septoria_linicola; do run_sp $sp & done
wait
$PY - <<'PY'
import re, glob, os
W = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench/domainfix"
def read(rep):
    out = {}
    for line in open(rep):
        f = line.rstrip('\n').split('\t'); m = re.search(r'\.(\w+)\.lib\.report$', f[0])
        if m and len(f) >= 11: out[m.group(1)] = float(f[6])
    return out
for sp in ['arabidopsis_thaliana', 'drosophila_melanogaster', 'zymoseptoria_tritici', 'septoria_linicola']:
    print(f"\n{sp}")
    for m in ['teannotate_ref7', 'teannotate_wo_ref7']:
        for score in ('default', 'harmonised'):
            line = f"  {m:20s} {score:10s}"
            for rule in ('R0', 'R1', 'R2'):
                d = f"{W}/{sp}/{m}/{rule}" + ("/h" if score == 'harmonised' else "")
                rep = glob.glob(f"{d}/*.fna.out.repeatmasker.lib.report")
                if not rep: line += f"  {rule}: -"; continue
                r = read(rep[0])
                line += f"  {rule}: " + " ".join(f"{c}={r.get(c,0):.3f}" for c in ('ltr', 'nonltr', 'tir', 'helitron'))
            print(line)
        rel = f"{W}/{sp}/{m}/R1/relabelled.tsv"
        if os.path.exists(rel):
            import collections
            c = collections.Counter((l.split('\t')[1].split('/')[0], l.split('\t')[2].strip().split('/')[0]) for l in open(rel))
            print(f"    R1 relabels: " + ", ".join(f"{a}->{b} x{n}" for (a, b), n in c.most_common(8)))
PY

