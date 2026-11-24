#!/usr/bin/env python3
"""variants_eval.py <sp>... - genome-level effect of re-classification variants on the frozen
v9 runs. All variants share the v9 library sequences (hence identical RepeatMasker hits and
"any" F1); they differ only in the labels of the consensi that were Unknown at the vote:
  none        no re-classification (they stay Unknown)
  v9          v9 rule (weighted, kingdom subsets of v9, min_agree 2) = the v9 library
  two_step    two-step vote, any Order (27320e2)
  restricted  two-step vote, Orders DNA,RC only (ad70d21, default) - run here with the
              pipeline voter on the same all-five predictions
-> vote_v10/variants/<sp>/<arm>/<variant>/ and vote_v10/results/variants.tsv"""
import os, subprocess, sys
from multiprocessing import Pool
from vote_v10_common import ARMS, SUBSET, V, paths, read_fasta
import bench_eval as BE

PY = "/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python"
VARIANTS = ("none", "v9", "two_step", "restricted")
CLASSES = ("any", "ltr", "nonltr", "tir", "helitron")


def labels(sp, arm):
    d = f"{V}/apply/{sp}/{arm}"
    rel = {}
    for i, l in enumerate(open(f"{d}/relabel.tsv")):
        f = l.rstrip("\n").split("\t")
        if i:
            rel[f[0]] = (f[1], f[2])
    org = paths(sp, arm)[1]
    preds = sorted(os.path.join(d, x) for x in os.listdir(d) if x.endswith(".preds.tsv"))
    r = subprocess.run([PY, f"{V}/bin/ensemble_vote.py", "--library", f"{d}/input_library.fa", "--pred", *preds,
                        "--out", f"{d}/restricted_library.fa", "--table", f"{d}/restricted_voting_table.tsv",
                        "--method", "two_step", "--kingdom", org, "--include-models", SUBSET[org],
                        "--relabel-orders", "DNA,RC"], capture_output=True, text=True)
    if r.returncode:
        sys.exit(r.stderr)
    res = {n.split()[0].partition("#")[0]: n.split()[0].partition("#")[2]
           for n, _ in read_fasta(f"{d}/restricted_library.fa")}
    return {"none": {k: "Unknown" for k in rel},
            "v9": {k: v[0] for k, v in rel.items()},
            "two_step": {k: v[1] for k, v in rel.items()},
            "restricted": {k: res[k] for k in rel}}


def main():
    jobs = []
    for sp in sys.argv[1:]:
        for arm in ARMS:
            v9out = paths(sp, arm)[4]
            if not (os.path.exists(v9out) and os.path.exists(f"{V}/apply/{sp}/{arm}/relabel.tsv")):
                continue
            lab = labels(sp, arm)
            for var in VARIANTS:
                d = f"{V}/variants/{sp}/{arm}/{var}"
                os.makedirs(d, exist_ok=True)
                with open(v9out) as fi, open(f"{d}/{sp}.fna.out", "w") as fo:
                    for i, line in enumerate(fi):
                        f = line.split()
                        if i >= 3 and len(f) >= 11 and f[9] in lab[var]:
                            f[10] = lab[var][f[9]]
                            line = " ".join(f) + "\n"
                        fo.write(line)
                jobs.append((sp, arm, var, d))
            n = {v: sum(1 for x in lab[v].values() if not x.lower().startswith("unknown")) for v in VARIANTS}
            print(f"{sp}/{arm}: Unknown-at-vote consensi {len(lab['none'])}, relabelled by "
                  + ", ".join(f"{v} {n[v]}" for v in VARIANTS[1:]), flush=True)
    with Pool(min(16, len(jobs))) as p:
        res = p.map(BE.score, jobs)
    os.makedirs(f"{V}/results", exist_ok=True)
    R = {(sp, arm, var): (dft, har, unk) for sp, arm, var, dft, har, unk in res}
    with open(f"{V}/results/variants{os.environ.get('TAG', '')}.tsv", "w") as fo:
        fo.write("species\tarm\tvariant\tmode\tclass\tF1\tprecision\tsensitivity\n")
        for (sp, arm, var), (dft, har, unk) in sorted(R.items()):
            for mode, rep in (("default", dft), ("harmonised", har)):
                for c, v in rep.items():
                    fo.write(f"{sp}\t{arm}\t{var}\t{mode}\t{c}\t{v['F1']:.4f}\t{v['prec']:.4f}\t{v['sens']:.4f}\n")
            fo.write(f"{sp}\t{arm}\t{var}\tall\tunknown_bp_share\t{unk:.4f}\t\t\n")
    print("\nharmonised F1 (with / w/o TEtrimmer): " + " | ".join(VARIANTS))
    for sp in sys.argv[1:]:
        for arm in ARMS:
            if (sp, arm, "v9") not in R:
                continue
            print(f"\n{sp} ({'with' if arm == 'w' else 'w/o'} TEtrimmer)  Unknown share "
                  + " / ".join(f"{100*R[(sp, arm, v)][2]:.1f}%" for v in VARIANTS))
            for c in CLASSES:
                print(f"  {c:9s} " + "  ".join(f"{R[(sp, arm, v)][1].get(c, {}).get('F1', float('nan')):.3f}"
                                               for v in VARIANTS))


if __name__ == "__main__":
    main()

