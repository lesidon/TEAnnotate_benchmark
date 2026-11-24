"""How much of each method's lib-test 'any' FP lies on repeats the REFERENCE annotates but
lib-test drops from its truth set (labels matching no class regex: Unknown, Satellite,
Cent, DNA/unknown, DNA/P, TcMar-Mariner, ...)? Reproduces lib-test's category regexes.
  default    : truth = reference annotations in a class (lib-test default, include_unknown off)
  don't-care : as default, but test bp on reference 'others' regions count neither TP nor FP
  inclusive  : truth = ALL reference annotations (lib-test --include_unknown behaviour)
"""
import re, sys, collections
import numpy as np

B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E = f"{B}/edta_manual_bench"
CATS = {
 'ltr': ['RLG', 'RLC', 'RLB', 'RLR', 'RLE', 'LTR', 'RLX', 'Gypsy', 'Copia'],
 'nonltr': ['SINE', 'LINE', 'Penelope', 'RIT', 'RIL', 'RST', 'RIX', 'RSX', 'nonLTR', 'YR'],
 'tir': ['TIR', 'hAT', 'hAT-Ac', 'MULE', 'MLE', 'MuDR', 'Tourist', 'CACT', 'PILE', 'POLE', 'Stowaway',
         'TcMar-Stowaway', 'PIF', 'Harbinger', 'Tc1', 'En-Spm', 'EnSpm', 'CMC-EnSpm', 'PiggyBac', 'Mirage',
         'P-element', 'Transib', 'DTA', 'DTH', 'DTT', 'DTM', 'DTC', 'DTX', 'DTR', 'DTE', 'Merlin', 'DTP',
         'DTB', 'polinton'],
 'mite': ['MITE'], 'helitron': ['Helitron', 'DHH', 'DHX', 'helitron']}
RX = re.compile('|'.join(sum(CATS.values(), [])), re.I)

def sizes(sp):
    king = [l.split('\t')[0] for l in open(f"{B}/config/genomes.tsv") if l.strip() and not l.startswith('#') and l.split('\t')[1].strip() == sp][0]
    s, name = {}, None
    for l in open(f"{B}/genomes/{king}/{sp}/{sp}.fna"):
        if l.startswith('>'): name = l[1:].split()[0]; s[name] = 0
        else: s[name] += len(l.strip())
    return s

def masks(out, size, split=False):
    a = {c: np.zeros(n, bool) for c, n in size.items()}
    o = {c: np.zeros(n, bool) for c, n in size.items()} if split else None
    lab_bp = collections.Counter()
    for i, l in enumerate(open(out)):
        f = l.split()
        if i < 3 or len(f) < 11 or f[4] not in a: continue
        s, e, lab = int(f[5]) - 1, int(f[6]), f[10]
        if split and not RX.search(lab) and lab != 'Unspecified':
            o[f[4]][s:e] = True; lab_bp[lab] += e - s
        else:
            a[f[4]][s:e] = True
    return a, o, lab_bp

sp = sys.argv[1]
size = sizes(sp)
ref_in, ref_oth, oth_labels = masks(f"{E}/{sp}/reference/{sp}.fna.out", size, split=True)
for c in ref_oth: ref_oth[c] &= ~ref_in[c]
oth_bp = sum(int(v.sum()) for v in ref_oth.values())
print(f"\n{sp}: reference repeats outside lib-test's truth set: {oth_bp/1e6:.2f} Mb; top labels:",
      ", ".join(f"{k} {v/1e3:.0f}kb" for k, v in oth_labels.most_common(6)))
methods = [("hite", f"{E}/{sp}/hite"), ("earlgrey", f"{E}/{sp}/earlgrey"), ("edta", f"{E}/{sp}/edta2"),
           ("teannotate w/o TT", f"{E}/{sp}/teannotate_wo_ref4"), ("teannotate TT", f"{E}/{sp}/teannotate_ref4"),
           ("evid_lite", f"{E}/pilot/{sp}/evid_lite"), ("evid_strict", f"{E}/pilot/{sp}/evid_strict"), ("tt_hybrid", f"{E}/pilot/{sp}/tt_hybrid")]
print(f"  {'method':18s} {'FP Mb':>6s} {'on ref-others':>13s} | F1 default  don't-care  inclusive")
for name, d in methods:
    try:
        t, _, _ = masks(f"{d}/{sp}.fna.out", size)
    except FileNotFoundError:
        continue
    tp = fp = fn = fpo = tpi = fni = 0
    for c in size:
        T, R, O = t[c], ref_in[c], ref_oth[c]
        tp += int((T & R).sum()); fn += int((R & ~T).sum())
        fp += int((T & ~R).sum()); fpo += int((T & O).sum())
        tpi += int((T & (R | O)).sum()); fni += int(((R | O) & ~T).sum())
    f_def = 2 * tp / (2 * tp + fp + fn)
    f_dc = 2 * tp / (2 * tp + (fp - fpo) + fn)
    f_inc = 2 * tpi / (2 * tpi + (fp - fpo) + fni)
    print(f"  {name:18s} {fp/1e6:6.2f} {fpo/1e6:7.2f} ({fpo/max(1,fp):3.0%}) |   {f_def:.3f}      {f_dc:.3f}      {f_inc:.3f}")
