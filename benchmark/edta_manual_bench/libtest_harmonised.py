#!/usr/bin/env python3
"""Harmonised lib-test: one label vocabulary for the reference and every library, then
EDTA lib-test as usual, plus an 'any_dontcare' row.

Why: lib-test's class regexes miss common RepeatMasker labels (bare 'DNA', 'DNA/P',
'DNA/TcMar*', 'DNA/Mariner', 'DNA/unknown', ...). In default mode, reference repeats
with such labels drop out of the truth set while test hits on them count as false
positives. Harmonised labels (applied identically to reference and tests):
  DNA transposons (any 'DNA', 'TIR', 'MITE' label)  -> 'TIR/<superfamily>' ('MITE/*' kept)
  any Helitron / 'RC'                                -> 'RC/Helitron'
  LTR / DIRS / ERV / Gypsy / Copia / Pao ...         -> 'LTR/<superfamily>'
  LINE / SINE / Penelope / nonLTR                    -> 'LINE/*' / 'SINE/*'
  non-TE (Satellite, Simple, Low_complexity, rRNA,
          tRNA, snRNA, Cent, rDNA, ...)              -> 'Unspecified' (never a TE class)
  everything else                                    -> 'Unknown'
any_dontcare: as lib-test 'any', but test bases on reference 'Unknown'/'Unspecified'
regions count neither as TP nor as FP (the reference annotates a repeat there but
cannot say which).

usage: libtest_harmonised.py <species> <method_dir>... (dirs under edta_manual_bench/<sp>/)
"""
import os, re, subprocess, sys
import numpy as np

B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E = f"{B}/edta_manual_bench"
PY = "/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python"
LIBTEST = f"{B}/tools/EDTA/lib-test.py"
LABEL = {"teannotate_ref9": "teannotate with tetrimmer (v9)", "teannotate_wo_ref9": "teannotate w/o tetrimmer (v9)",
         "teannotate_ref8": "teannotate with tetrimmer (v8)", "teannotate_wo_ref8": "teannotate w/o tetrimmer (v8)",
         "teannotate_ref7": "teannotate with tetrimmer (v7)", "teannotate_wo_ref7": "teannotate w/o tetrimmer (v7)",
         "teannotate_ref6": "teannotate with tetrimmer (v6)", "teannotate_wo_ref6": "teannotate w/o tetrimmer (v6)",
         "teannotate_ref5": "teannotate with tetrimmer", "teannotate_wo_ref5": "teannotate w/o tetrimmer",
         "teannotate_ref4": "teannotate with tetrimmer (v4)", "teannotate_wo_ref4": "teannotate w/o tetrimmer (v4)",
         "edta2": "edta", "edta": "edta (rice fallback)", "hite": "hite", "earlgrey": "earlgrey", "heliano": "heliano"}

NON_TE = re.compile(r"^(satellite|simple_repeat|simple|low_complexity|rrna|trna|snrna|scrna|srprna|"
                    r"cent|centr|rdna|artefact|artifact|segmental|buffer|other)", re.I)

def canon(lab):
    if lab == "Unspecified":
        return lab
    l = lab.lower()
    top, _, sub = lab.partition("/")
    tl = top.lower()
    if NON_TE.match(l):
        return "Unspecified"
    if "helitron" in l or tl == "rc":
        return "RC/Helitron"
    if tl in ("ltr", "dirs", "erv") or re.search(r"gypsy|copia|bel-pao|retrovirus|caulimo|ngaro", l):
        return lab if (tl == "ltr" and sub) else "LTR/" + (sub or (top if tl != "ltr" else "unknown"))
    if tl in ("line", "sine"):
        return lab if sub else f"{top}/unknown"
    if tl in ("ple", "penelope"):
        return "LINE/Penelope"
    if tl == "nonltr":
        return "LINE/unknown"
    if tl == "retroposon":
        return "SINE/" + (sub or "unknown")
    if tl == "mite":
        return lab
    if tl in ("dna", "tir", "class_ii", "dnaauto", "dnanona"):
        if re.search(r"maverick|polinton", l):
            return "TIR/Maverick"
        return "TIR/" + (sub or "unknown")
    return "Unknown"

def harmonise(src, dst):
    with open(src) as fi, open(dst, "w") as fo:
        for i, line in enumerate(fi):
            f = line.split()
            if i < 3 or len(f) < 11:
                fo.write(line); continue
            f[10] = canon(f[10])
            fo.write(" ".join(f) + "\n")

def sizes(sp):
    king = [l.split("\t")[0] for l in open(f"{B}/config/genomes.tsv")
            if l.strip() and not l.startswith("#") and l.split("\t")[1].strip() == sp][0]
    g = f"{B}/genomes/{king}/{sp}/{sp}.fna"
    s, name = {}, None
    for l in open(g):
        if l.startswith(">"): name = l[1:].split()[0]; s[name] = 0
        else: s[name] += len(l.strip())
    return g, s

def masks(out, size, ref=False):
    a = {c: np.zeros(n, bool) for c, n in size.items()}
    o = {c: np.zeros(n, bool) for c, n in size.items()} if ref else None
    for i, l in enumerate(open(out)):
        f = l.split()
        if i < 3 or len(f) < 11 or f[4] not in a: continue
        s, e = int(f[5]) - 1, int(f[6])
        if ref and f[10] in ("Unknown", "Unspecified"):
            o[f[4]][s:e] = True
        else:
            a[f[4]][s:e] = True
    return a, o

def main():
    sp, methods = sys.argv[1], sys.argv[2:]
    genome, size = sizes(sp)
    rdir = f"{E}/{sp}/reference/harmonised"; os.makedirs(rdir, exist_ok=True)
    ref = f"{rdir}/{sp}.fna.out"
    harmonise(f"{E}/{sp}/reference/{sp}.fna.out", ref)
    ref_in, ref_oth = masks(ref, size, ref=True)
    for c in ref_oth: ref_oth[c] &= ~ref_in[c]
    rows = []
    for m in methods:
        src = f"{E}/{sp}/{m}/{sp}.fna.out"
        if not os.path.exists(src):
            continue
        d = f"{E}/{sp}/{m}/harmonised"; os.makedirs(d, exist_ok=True)
        tst = f"{d}/{sp}.fna.out"
        harmonise(src, tst)
        with open(f"{d}/libtest.log", "w") as log:
            subprocess.run([PY, LIBTEST, "--genome", genome, "--reference", ref, "--test", tst,
                            "--extended_report"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
        rep = f"{d}/{sp}.fna.out.repeatmasker.lib.report"
        label = LABEL.get(m, m)
        for line in open(rep):
            f = line.rstrip("\n").split("\t")
            mm = re.search(r"\.(\w+)\.lib\.report$", f[0])
            if mm and len(f) >= 11:
                rows.append([label, mm.group(1)] + f[1:11])
        # any_dontcare
        t, _ = masks(tst, size)
        tp = fp = fn = tn = 0
        for c in size:
            T, R, O = t[c], ref_in[c], ref_oth[c]
            tp += int((T & R).sum()); fn += int((R & ~T).sum())
            fp += int((T & ~R & ~O).sum()); tn += int((~T & ~R & ~O).sum())
        sens = tp / max(1, tp + fn); prec = tp / max(1, tp + fp); spec = tn / max(1, tn + fp)
        f1 = 2 * tp / max(1, 2 * tp + fp + fn)
        rows.append([label, "any_dontcare", f"{sens}", f"{spec}", f"{(tp+tn)/max(1,tp+tn+fp+fn)}",
                     f"{prec}", f"{1-prec}", f"{f1}", str(tp), str(tn), str(fp), str(fn)])
        print(f"{sp}/{m}: any_dontcare F1 {f1:.3f}", flush=True)
    nice = sp[0].upper() + sp[1:]
    out = f"{E}/results/{nice}_bench_harmonised.csv"
    with open(out, "w") as o:
        o.write("Method,TE type,sens,spec,accu,prec,FDR,F1,TP,TN,FP,FN\n")
        for r in rows:
            o.write(",".join(r) + "\n")
    print(out, len(rows), "rows")

if __name__ == "__main__":
    main()
