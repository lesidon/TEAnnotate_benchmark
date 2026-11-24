#!/usr/bin/env python3
"""Convert a RepBase FASTA (.ref) library into a RepeatMasker-ready library.

RepBase FASTA headers are '>name<TAB>class<TAB>species' (RepBase's own class
vocabulary), which RepeatMasker -lib does not read — so, like the Dfam FASTA, a
library fed straight in is annotated 'Unknown' and lib-test.py cannot split it into
TE categories. This rewrites each header to '>name#Class' with the RepBase class
mapped to a RepeatMasker-style token (LTR / LINE / SINE / DNA / RC / Satellite / …).

Tolerant: headers already containing '#' are kept; if there is no class field the
entry becomes '#Unknown'; falls back to splitting on runs of >=2 spaces when the
header is not tab-delimited.

Usage:
  python repbase_to_rmlib.py --fasta drorep.ref --out drosophila.curated.lib.fa
"""
import argparse
import gzip
import re
import sys

_MULTISPACE = re.compile(r"\s{2,}")


def load_text(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    head = raw.lstrip()[:1]
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode("utf-8", "replace")
    if head == b"<":
        raise ValueError(f"{path} looks like an HTML page, not a RepBase FASTA")
    if head == b"{":
        raise ValueError(f"{path} looks like JSON, not a RepBase FASTA")
    return raw.decode("utf-8", "replace")


def rm_class(raw):
    """Map a RepBase class phrase to a RepeatMasker-style Class token."""
    r = (raw or "").strip().lower()
    if not r:
        return "Unknown"
    # Order matters: test 'non-LTR' before 'LTR'.
    if "non-ltr" in r or "nonltr" in r or "line" in r or "l1" in r or "cr1" in r or "jockey" in r:
        return "LINE"
    if "sine" in r or "alu" in r:
        return "SINE"
    if ("ltr" in r or "endogenous retrovirus" in r or "erv" in r or "gypsy" in r
            or "copia" in r or "bel" in r or "pao" in r or "dirs" in r):
        return "LTR"
    if "helitron" in r or "rolling" in r or r == "rc":
        return "RC/Helitron"
    if "maverick" in r or "polinton" in r:
        return "DNA/Maverick"
    if ("dna transposon" in r or "terminal inverted" in r or r == "tir" or "mite" in r
            or "hat" in r or "mariner" in r or "tc1" in r or "mudr" in r or "mutator" in r
            or "harbinger" in r or "pif" in r or "piggybac" in r or "transib" in r
            or "cacta" in r or "enspm" in r or "p element" in r or "transposon" in r):
        return "DNA"
    if "satellite" in r:
        return "Satellite"
    if "simple" in r or "microsatellite" in r or "tandem" in r:
        return "Simple_repeat"
    # Unrecognised: keep the original token (spaces -> underscores) so nothing is lost.
    return re.sub(r"\s+", "_", raw.strip())


def parse_header(line):
    body = line[1:].rstrip("\n")
    parts = body.split("\t")
    if len(parts) < 2:
        parts = _MULTISPACE.split(body)
    name = parts[0].strip()
    cls = parts[1].strip() if len(parts) >= 2 else ""
    return name, cls


def convert(fasta_in, out_path):
    n = unknown = kept = 0
    text = load_text(fasta_in)
    with open(out_path, "w", encoding="utf-8") as out:
        for line in text.splitlines():
            if line.startswith(">"):
                n += 1
                if "#" in line:                      # already RepeatMasker-style
                    kept += 1
                    out.write(line + "\n")
                    continue
                name, cls_raw = parse_header(line)
                cls = rm_class(cls_raw)
                if cls == "Unknown":
                    unknown += 1
                out.write(f">{name.replace(' ', '_')}#{cls}\n")
            else:
                out.write(line + "\n")
    return n, unknown, kept


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="RepBase FASTA (.ref) input")
    ap.add_argument("--out", required=True, help="RepeatMasker-ready library FASTA")
    args = ap.parse_args(argv)

    n, unknown, kept = convert(args.fasta, args.out)
    sys.stderr.write(f"[repbase_to_rmlib] {n} entries -> {args.out}; "
                     f"{kept} already RM-style, {unknown} left Unknown\n")
    if n == 0:
        sys.stderr.write("[repbase_to_rmlib] ERROR: no FASTA records parsed "
                         "(HTML/JSON/empty input?)\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
