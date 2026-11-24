#!/usr/bin/env python3
"""Build a species-level curated TE library from the RepBase FASTA release.

Selects every RepBase entry whose taxon (third header field) is the target species or one
of its ancestral clades in the NCBI taxonomy lineage. This mirrors how RepeatMasker
applies RepBase (-species includes ancestral clades): mouse = Mus musculus + Mus + Murinae
+ ... + Mammalia + Eutheria + Vertebrata ... A species-only selection would leave the
ancestral repeats (MIR, L2, MER families) out of the truth set, so tools that find them
would be scored as false positives.

RepBase headers are '>name<TAB>class<TAB>taxon'. Each class is mapped to a RepeatMasker
'Class/Superfamily' token (name#Class/Superfamily) so lib-test.py and the harmonised scorer
can split the library by TE category. A side TSV records name, taxon, RepBase class, label.

usage: repbase_species_lib.py --repbase-dir DIR --taxid N --out lib.fa [--tsv out.tsv]
                              [--lineage "Name1;Name2;..."]   (skips the NCBI lookup)
"""
import argparse
import glob
import os
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

EXACT = {
    # LTR retrotransposons
    "gypsy": "LTR/Gypsy", "copia": "LTR/Copia", "bel": "LTR/Pao", "dirs": "LTR/DIRS",
    "ngaro": "LTR/Ngaro", "erv1": "LTR/ERV1", "erv2": "LTR/ERVK", "erv3": "LTR/ERVL",
    "erv4": "LTR/ERV4", "endogenous retrovirus": "LTR/ERV", "ltr retrotransposon": "LTR/unknown",
    "caulimoviridae": "LTR/Caulimovirus", "lentivirus": "LTR/ERV", "lokiretrovirus": "LTR/ERV",
    # non-LTR retrotransposons
    "l1": "LINE/L1", "tx1": "LINE/Tx1", "l2": "LINE/L2", "l2a": "LINE/L2", "l2b": "LINE/L2",
    "cr1": "LINE/CR1", "rte": "LINE/RTE", "rtex": "LINE/RTE-X", "jockey": "LINE/I-Jockey",
    "i": "LINE/I", "ingi": "LINE/Ingi", "r1": "LINE/R1", "r2": "LINE/R2", "r4": "LINE/R4",
    "nesl": "LINE/R2-NeSL", "tad1": "LINE/Tad1", "daphne": "LINE/Daphne", "rex1": "LINE/Rex-Babar",
    "proto1": "LINE/Proto1", "proto2": "LINE/Proto2", "nimb": "LINE/Nimb", "kiri": "LINE/Kiri",
    "crack": "LINE/Crack", "loa": "LINE/Loa", "cre": "LINE/CRE", "hero": "LINE/Hero",
    "outcast": "LINE/Outcast", "vingi": "LINE/Vingi", "randi": "LINE/Dualen",
    "non-ltr retrotransposon": "LINE/unknown",
    "penelope": "LINE/Penelope", "penelope/poseidon": "LINE/Penelope", "neptune": "LINE/Penelope",
    "athena": "LINE/Penelope", "naiad/chlamys": "LINE/Penelope",
    "sine": "SINE/unknown", "sine1/7sl": "SINE/7SL", "sine2/trna": "SINE/tRNA",
    "sine3/5s": "SINE/5S", "sineu/snrna": "SINE/U",
    # DNA transposons
    "mariner/tc1": "DNA/TcMar", "hat": "DNA/hAT", "charlie-galluhop": "DNA/hAT-Charlie",
    "mudr": "DNA/MULE-MuDR", "harbinger": "DNA/PIF-Harbinger", "enspm/cacta": "DNA/CMC-EnSpm",
    "enspm": "DNA/CMC-EnSpm", "piggybac": "DNA/PiggyBac", "p": "DNA/P", "merlin": "DNA/Merlin",
    "transib": "DNA/CMC-Transib", "kolobok": "DNA/Kolobok", "academ": "DNA/Academ",
    "sola1": "DNA/Sola", "sola2": "DNA/Sola", "sola3": "DNA/Sola", "dada": "DNA/Dada",
    "zator": "DNA/Zator", "isl2eu": "DNA/ISL2EU", "is3eu": "DNA/IS3EU", "is481eu": "DNA/IS481EU",
    "ginger1": "DNA/Ginger", "ginger2/tdd": "DNA/Ginger", "novosib": "DNA/Novosib",
    "zisupton": "DNA/Zisupton", "polinton": "DNA/Maverick", "dna transposon": "DNA/unknown",
    "transposon": "DNA/unknown", "replitron": "DNA/Replitron", "mite": "MITE/unknown",
    "crypton": "DNA/Crypton", "cryptonv": "DNA/Crypton", "cryptons": "DNA/Crypton",
    "cryptona": "DNA/Crypton", "cryptoni": "DNA/Crypton", "cryptonf": "DNA/Crypton",
    "helitron": "RC/Helitron",
    # non-TE repeats
    "satellite": "Satellite", "sat": "Satellite", "msat": "Satellite", "tandem repeat": "Satellite",
    "simple repeat": "Simple_repeat", "(ca)n related": "Simple_repeat", "trna": "tRNA",
    "rrna": "rRNA", "ribosomal intergenic spacer": "rRNA", "snrna": "snRNA",
    "multicopy gene": "Other", "integrated virus": "Other", "dna virus": "Other",
    "y chromosome": "Other",
}


def rm_label(cls):
    c = (cls or "").strip().lower()
    if c in EXACT:
        return EXACT[c]
    if c.startswith("crypton"):
        return "DNA/Crypton"
    if re.match(r"^(eutr|eutrep|mare)\d*$", c):    # eutherian / mammalian repeats of unknown type
        return "Unknown"
    if "satellite" in c:
        return "Satellite"
    if "ltr" in c and "non" not in c:
        return "LTR/unknown"
    if "sine" in c:
        return "SINE/unknown"
    if "retrotransposon" in c or "line" in c:
        return "LINE/unknown"
    if "dna" in c or "transposon" in c:
        return "DNA/unknown"
    return "Unknown"


def ncbi_lineage(taxid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml"
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                root = ET.fromstring(r.read())
            t = root.find("Taxon")
            names = [x.findtext("ScientificName") for x in t.find("LineageEx")]
            return t.findtext("ScientificName"), names
        except Exception as e:  # network hiccup: retry
            print(f"lineage lookup failed ({e}); retry {attempt + 1}", file=sys.stderr)
            time.sleep(10)
    sys.exit(f"could not fetch the NCBI lineage for taxid {taxid}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repbase-dir", required=True)
    ap.add_argument("--taxid")
    ap.add_argument("--lineage", help="';'-separated species + ancestral clade names")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tsv")
    a = ap.parse_args()
    if a.lineage:
        names = [x.strip() for x in a.lineage.split(";") if x.strip()]
        species = names[0]
    else:
        species, lineage = ncbi_lineage(a.taxid)
        names = [species] + lineage
    wanted = {n.lower() for n in names}

    seqs, meta, order = {}, {}, []
    for f in sorted(glob.glob(os.path.join(a.repbase_dir, "*.ref"))):
        name = None
        for line in open(f, errors="replace"):
            if line.startswith(">"):
                p = line[1:].rstrip("\n").split("\t")
                name = None
                if len(p) >= 3 and p[2].strip().lower() in wanted:
                    nm = re.sub(r"[^\w.\-:+()]", "_", p[0].strip())
                    if nm in seqs:            # same entry in several taxon files: keep one
                        continue
                    name = nm
                    seqs[name] = []
                    meta[name] = (p[2].strip(), p[1].strip(), rm_label(p[1]))
                    order.append(name)
            elif name:
                seqs[name].append(line.strip())
    if not order:
        sys.exit(f"no RepBase entries for {species} (lineage {names})")

    with open(a.out, "w") as o:
        for n in order:
            s = "".join(seqs[n]).upper()
            if not s:
                continue
            o.write(f">{n}#{meta[n][2]}\n")
            for i in range(0, len(s), 80):
                o.write(s[i:i + 80] + "\n")
    if a.tsv:
        with open(a.tsv, "w") as t:
            t.write("name\ttaxon\trepbase_class\trm_label\n")
            for n in order:
                t.write("\t".join((n,) + meta[n]) + "\n")
    from collections import Counter
    print(f"{species}: {len(order)} entries -> {a.out}")
    print("  taxa:", Counter(m[0] for m in meta.values()).most_common())
    print("  labels:", Counter(m[2] for m in meta.values()).most_common())
    unk = Counter(m[1] for m in meta.values() if m[2] == "Unknown")
    if unk:
        print("  RepBase classes mapped to Unknown:", unk.most_common(20))


if __name__ == "__main__":
    main()

