"""Per-consensus intrinsic features (computable inside the pipeline, no reference needed)
joined with lib-test's per-sequence TP/FP (the answer key, used only to evaluate).

usage: junk_features.py <species> <method> <support.tsv> <out.tsv>
features from the method's own genome masking (<sp>.fna.out):
  n_hits, masked_bp, full_frac (hits covering >=80% of the consensus), n_full,
  med_div (median % divergence), array_frac (hits whose next hit of the same consensus
  starts within 50 bp -> tandem arrays), cons_len, label, source tool, n_support
"""
import re, sys, statistics, collections

sp, method, support_file, out = sys.argv[1:5]
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
D = f"{E}/{sp}/{method}"

def source(name):
    n = name.replace("_", "-")
    if re.match(r"TE-\d+", n): return "EDTA"
    if re.match(r"(LTR|TIR|Helitron|Homology-Non-LTR|Denovo-Non-LTR)-\d+", n): return "HiTE"
    if n.startswith("insertion-HLE"): return "HELIANO"
    if n.startswith("rnd-"): return "EarlGrey"
    return "other"

lib = {}
name = None
for l in open(f"{D}/library.fa"):
    if l.startswith(">"):
        name, _, label = l[1:].split()[0].partition("#")
        lib[name] = [label or "NA", 0]
    elif name:
        lib[name][1] += len(l.strip())

hits = collections.defaultdict(list)
for i, l in enumerate(open(f"{D}/{sp}.fna.out")):
    f = l.split()
    if i < 3 or len(f) < 14:
        continue
    chrom, qb, qe, strand, rep = f[4], int(f[5]), int(f[6]), f[8], f[9]
    if strand == "+":
        rb, re_ = int(f[11]), int(f[12])
    else:
        rb, re_ = int(f[13]), int(f[12])
    hits[rep].append((chrom, qb, qe, abs(re_ - rb) + 1, float(f[1])))

ans = {}
rep = open(f"{D}/{sp}.fna.out.sequence_assessment.report").read().splitlines()
hdr = rep[0].split("\t")
for l in rep[1:]:
    d = dict(zip(hdr, l.split("\t")))
    ans[d["sequence_id"]] = (int(d["tp_all_cats"]), int(d["fp_all_cats"]))

support = collections.defaultdict(int)
for l in open(support_file):
    q, n = l.split()
    support[q] = int(n)

with open(out, "w") as o:
    o.write("species\tmethod\tname\tsource\tlabel\tcons_len\tn_hits\tmasked_bp\tn_full\tfull_frac\tmed_div\tarray_frac\tn_support\ttp\tfp\n")
    for name, (label, ln) in lib.items():
        h = sorted(hits.get(name, []))
        if not h:
            continue
        n_full = sum(1 for x in h if x[3] >= 0.8 * ln)
        arr = 0
        for a, b in zip(h, h[1:]):
            if a[0] == b[0] and b[1] - a[2] <= 50:
                arr += 1
        masked = sum(x[2] - x[1] + 1 for x in h)
        tp, fp = ans.get(name, (0, 0))
        o.write(f"{sp}\t{method}\t{name}\t{source(name)}\t{label}\t{ln}\t{len(h)}\t{masked}\t{n_full}\t"
                f"{n_full/len(h):.3f}\t{statistics.median(x[4] for x in h):.1f}\t{arr/max(1,len(h)-1):.3f}\t"
                f"{support[name]}\t{tp}\t{fp}\n")
print(f"{sp}/{method}: {len(lib)} consensi -> {out}")
