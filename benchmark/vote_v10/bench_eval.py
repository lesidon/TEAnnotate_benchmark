#!/usr/bin/env python3
"""bench_eval.py <sp>... - genome-level effect of the v10 vote. The v10 library has the same
sequences as v9, only some Unknown labels changed, so RepeatMasker hits are identical and
the v10 annotation is the v9 bench .out with the class column of relabelled consensi
replaced. Both are scored with EDTA lib-test, default and harmonised (same procedure and
reference as the bench) -> vote_v10/bench/<sp>/<arm>/{v9,v10}/ and vote_v10/results/."""
import os, re, shutil, subprocess, sys
from multiprocessing import Pool
from vote_v10_common import ARMS, E, V, paths

sys.path.insert(0, E)
import libtest_harmonised as LH  # noqa: E402  (harmonise(), sizes(), LIBTEST, PY)


def report(path):
    rows = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        m = re.search(r"\.(\w+)\.lib\.report$", f[0])
        if m and len(f) >= 11:
            rows[m.group(1)] = {"sens": float(f[1]), "prec": float(f[4]), "F1": float(f[6])}
    return rows


def unknown_share(out):
    tot = unk = 0
    for i, l in enumerate(open(out)):
        f = l.split()
        if i < 3 or len(f) < 11:
            continue
        L = int(f[6]) - int(f[5]) + 1
        tot += L
        unk += L if f[10].lower().startswith("unknown") else 0
    return unk / max(1, tot)


def score(job):
    sp, arm, ver, d = job
    genome, _ = LH.sizes(sp)
    out = f"{d}/{sp}.fna.out"
    with open(f"{d}/libtest.log", "w") as log:
        subprocess.run([LH.PY, LH.LIBTEST, "--genome", genome, "--reference", f"{E}/{sp}/reference/{sp}.fna.out",
                        "--test", out, "--extended_report"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    h = f"{d}/harmonised"; os.makedirs(h, exist_ok=True)
    LH.harmonise(out, f"{h}/{sp}.fna.out")
    with open(f"{h}/libtest.log", "w") as log:
        subprocess.run([LH.PY, LH.LIBTEST, "--genome", genome, "--reference",
                        f"{E}/{sp}/reference/harmonised/{sp}.fna.out", "--test", f"{h}/{sp}.fna.out",
                        "--extended_report"], cwd=h, stdout=log, stderr=subprocess.STDOUT)
    return (sp, arm, ver, report(f"{d}/{sp}.fna.out.repeatmasker.lib.report"),
            report(f"{h}/{sp}.fna.out.repeatmasker.lib.report"), unknown_share(out))


def main():
    jobs = []
    for sp in sys.argv[1:]:
        for arm in ARMS:
            v9out = paths(sp, arm)[4]
            rel = f"{V}/apply/{sp}/{arm}/relabel.tsv"
            if not (os.path.exists(v9out) and os.path.exists(rel)):
                continue
            new = {}
            for i, l in enumerate(open(rel)):
                f = l.rstrip("\n").split("\t")
                if i and f[2] and f[2] != f[1]:
                    new[f[0]] = f[2]
            base = f"{V}/bench/{sp}/{arm}"
            for ver in ("v9", "v10"):
                d = f"{base}/{ver}"; os.makedirs(d, exist_ok=True)
                if ver == "v9":
                    shutil.copy(v9out, f"{d}/{sp}.fna.out")
                    continue
                hit = 0
                with open(v9out) as fi, open(f"{d}/{sp}.fna.out", "w") as fo:
                    for i, line in enumerate(fi):
                        f = line.split()
                        if i >= 3 and len(f) >= 11 and f[9] in new:
                            f[10] = new[f[9]]; hit += 1
                            line = " ".join(f) + "\n"
                        fo.write(line)
                print(f"{sp}/{arm}: {len(new)} consensi relabelled vs v9, {hit} RepeatMasker hits changed", flush=True)
            jobs += [(sp, arm, ver, f"{base}/{ver}") for ver in ("v9", "v10")]
    with Pool(min(16, len(jobs))) as p:
        res = p.map(score, jobs)
    os.makedirs(f"{V}/results", exist_ok=True)
    R = {(sp, arm, ver): (dft, har, unk) for sp, arm, ver, dft, har, unk in res}
    with open(f"{V}/results/vote_v10_libtest{os.environ.get('TAG', '')}.tsv", "w") as fo:
        fo.write("species\tarm\tversion\tmode\tclass\tF1\tprecision\tsensitivity\n")
        for (sp, arm, ver), (dft, har, unk) in sorted(R.items()):
            for mode, rep in (("default", dft), ("harmonised", har)):
                for c, v in rep.items():
                    fo.write(f"{sp}\t{arm}\t{ver}\t{mode}\t{c}\t{v['F1']:.4f}\t{v['prec']:.4f}\t{v['sens']:.4f}\n")
            fo.write(f"{sp}\t{arm}\t{ver}\tall\tunknown_bp_share\t{unk:.4f}\t\t\n")
    for sp in sys.argv[1:]:
        for arm in ARMS:
            if (sp, arm, "v9") not in R:
                continue
            a, b = R[(sp, arm, "v9")], R[(sp, arm, "v10")]
            print(f"\n{sp} ({'with' if arm == 'w' else 'w/o'} TEtrimmer): Unknown share of masked bp "
                  f"{100*a[2]:.1f}% -> {100*b[2]:.1f}%")
            for mode, k in (("harmonised", 1), ("default", 0)):
                print(f"  {mode:10s} " + "  ".join(
                    f"{c} {a[k][c]['F1']:.3f}->{b[k].get(c, {}).get('F1', float('nan')):.3f}" for c in sorted(a[k])))


if __name__ == "__main__":
    main()

