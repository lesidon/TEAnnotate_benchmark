#!/usr/bin/env python3
"""build_dataset.py - candidate sequences for the extended ensemble benchmark.

MycoMobilome v1.1 (clustered 80-80-80 protein-evidence library): keeps consensi whose
classification is supported by TE protein hits (_PE; _DA = contradicting evidence and
_NE = none are dropped) and resolves to a TE order. Sequences from the held-out genomes
of the whole-pipeline benchmark are flagged (heldout_species=1) so the split puts them,
and every sequence clustered with them, in the ensemble TEST set only: tuned weights then
never see a genome the pipeline is evaluated on blind.

Writes <outdir>/<prefix>.fasta (plain ids MM000001..., the classifiers' input) and
<outdir>/<prefix>.tsv (truth: family_id, species, kingdom, label, label_order,
label_superfamily, source, evidence, orig_id, seq_length, heldout_species).

usage: build_dataset.py --mm-fasta <MycoMobilome ...proteinEvidence...fasta> --outdir data
"""
import argparse
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from te_label_utils import load_label_overrides, normalize_label  # noqa: E402

# held-out genomes of the whole-pipeline benchmark; a one-word entry flags the whole genus
# (Blumeria hordei / graminis, Ectocarpus spp.)
HELDOUT_SPECIES = [
    "Caenorhabditis elegans", "Caenorhabditis briggsae", "Blumeria",
    "Chlamydomonas reinhardtii", "Ectocarpus", "Gallus gallus", "Taeniopygia guttata",
    "Mus musculus",
]
NON_TE_ORDERS = {"Unknown", "Other"}
HDR = re.compile(r"^>(\S+?)#(\S+)(?:\s+@(.+))?")


def is_heldout(species):
    s = " ".join((species or "").strip().lower().replace("_", " ").split())
    for b in HELDOUT_SPECIES:
        b = b.lower()
        if s == b or s.startswith(b + " "):
            return True
    return False


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line.rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mm-fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="mycomobilome_pe")
    ap.add_argument("--label-map", default=None)
    ap.add_argument("--min-length", type=int, default=100)
    a = ap.parse_args()
    load_label_overrides(a.label_map)
    os.makedirs(a.outdir, exist_ok=True)
    stats = Counter()
    orders, species_seen, heldout_sp = Counter(), set(), Counter()
    n = 0
    with open(os.path.join(a.outdir, a.prefix + ".fasta"), "w") as fa, \
         open(os.path.join(a.outdir, a.prefix + ".tsv"), "w") as tv:
        tv.write("family_id\tspecies\tkingdom\tlabel\tlabel_order\tlabel_superfamily\t"
                 "source\tevidence\torig_id\tseq_length\theldout_species\n")
        for hdr, seq in read_fasta(a.mm_fasta):
            stats["read"] += 1
            m = HDR.match(hdr)
            if not m:
                stats["unparsed"] += 1
                continue
            orig, label, species = m.group(1), m.group(2), (m.group(3) or "").strip()
            ev = orig.rsplit("_", 1)[-1]
            if ev != "PE":
                stats["evidence_" + ev] += 1
                continue
            order = normalize_label(label, "order")
            if order in NON_TE_ORDERS:
                stats["non_te_or_unknown"] += 1
                continue
            ho = is_heldout(species)
            if ho:
                heldout_sp[species] += 1
            if len(seq) < a.min_length:
                stats["short"] += 1
                continue
            n += 1
            fid = "MM%06d" % n
            fa.write(">%s\n" % fid)
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            tv.write("\t".join([fid, species, "fungi", label, order,
                                normalize_label(label, "superfamily"), "MycoMobilome_v1.1",
                                ev, orig, str(len(seq)), "1" if ho else "0"]) + "\n")
            orders[order] += 1
            species_seen.add(species)
    print("records:", dict(stats))
    print("kept: %d consensi from %d species" % (n, len(species_seen)))
    print("orders:", orders.most_common())
    if heldout_sp:
        print("held-out species (forced to the test split):", heldout_sp.most_common())


if __name__ == "__main__":
    main()

