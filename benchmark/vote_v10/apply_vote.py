#!/usr/bin/env python3
"""apply_vote.py <sp>... - replay the ensemble step of each v9 run with the new pipeline
voter (vote_v10/bin/ensemble_vote.py = commit 27320e2): the v9 library with its
Unknown-at-vote consensi reset to Unknown + all-five predictions (vote_v10/preds/<org>) ->
vote_v10/apply/<sp>/<arm>/{v10_library.fa, v10_voting_table.tsv, relabel.tsv}.
Everything upstream of the vote is identical to v9, so this is the v10 library."""
import glob, os, subprocess, sys
from collections import Counter
from vote_v10_common import ARMS, SUBSET, V, paths, read_fasta

PY = "/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python"
for sp in sys.argv[1:]:
    for arm in ARMS:
        pre, org, lib, vt, _ = paths(sp, arm)
        if not (os.path.exists(lib) and os.path.exists(vt)):
            continue
        d = f"{V}/apply/{sp}/{arm}"
        os.makedirs(d, exist_ok=True)
        voted = {l.split("\t")[0] for i, l in enumerate(open(vt)) if i}
        v9 = {}
        with open(f"{d}/input_library.fa", "w") as fh:
            for name, seq in read_fasta(lib):
                tok, _, desc = name.partition(" ")
                base, _, lab = tok.partition("#")
                v9[base] = lab
                if base in voted or lab.lower().startswith("unknown"):
                    tok = base + "#Unknown"
                fh.write(">" + tok + (" " + desc if desc else "") + "\n" + "\n".join(seq) + "\n")
        preds = []
        tag = os.environ.get("TAG", "")
        for f in sorted(glob.glob(f"{V}/preds/{org}{tag}/**/*.preds.tsv", recursive=True)):
            model = os.path.basename(f)[:-len(".preds.tsv")].rsplit(".", 1)[-1]
            out = f"{d}/{sp}.{model}.preds.tsv"
            n = 0
            with open(out, "w") as fo:
                for l in open(f):
                    if l.startswith(pre + "_"):
                        fo.write(l[len(pre) + 1:]); n += 1
            preds.append(out)
        r = subprocess.run([PY, f"{V}/bin/ensemble_vote.py", "--library", f"{d}/input_library.fa",
                            "--pred", *preds, "--out", f"{d}/v10_library.fa",
                            "--table", f"{d}/v10_voting_table.tsv", "--method", "all",
                            "--kingdom", org, "--include-models", SUBSET[org], "--min-agree", "2",
                            "--relabel-orders", "all"],     # unrestricted two-step (27320e2)
                           capture_output=True, text=True)
        if r.returncode:
            print(r.stderr); sys.exit(1)
        v10 = {n.split()[0].partition("#")[0]: n.split()[0].partition("#")[2]
               for n, _ in read_fasta(f"{d}/v10_library.fa")}
        trans = Counter()
        with open(f"{d}/relabel.tsv", "w") as fo:
            fo.write("repeat_id\tv9_label\tv10_label\n")
            for b, lab in v9.items():
                if b in voted or lab.lower().startswith("unknown") or v10.get(b) != lab:
                    fo.write(f"{b}\t{lab}\t{v10.get(b)}\n")
                    trans[(lab.split('/')[0], v10.get(b, '').split('/')[0])] += 1
        print(f"{sp}/{arm} ({org}): {r.stderr.strip().splitlines()[-1]}")
        print("   v9 -> v10 order:", ", ".join(f"{a}->{b}:{n}" for (a, b), n in trans.most_common()))

