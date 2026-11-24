#!/usr/bin/env python3
"""Normalise a curated TE library to RepeatMasker 'name#Class/Superfamily' headers.

Keeps the curators' classification and only repairs spellings that RepeatMasker,
lib-test.py or the harmonised scorer would misread:
  'LINE/Tad1?' -> 'LINE/Tad1'        (uncertainty marks)
  'TIR/Tc1/Mariner' -> 'DNA/TcMar-Mariner'
  'DNA/Helitron' -> 'RC/Helitron'
  'Simple/repeat' -> 'Simple_repeat', 'Low/complexity' -> 'Low_complexity'
  'tandem_repeat' -> 'Satellite', 'DIRS/x' -> 'LTR/x', 'PLE/x' -> 'LINE/Penelope'
  'Crypton/x' -> 'DNA/Crypton', bare 'RNA' (rDNA / 7S rRNA entries) -> 'rRNA'
                                     (RepeatMasker's own conventions)
  bare 'LTR' / 'LINE' / 'SINE' / 'DNA' -> '<class>/unknown'
Names keep their first whitespace-free token (made unique). Prints label counts before
and after so every change is visible.

usage: normalize_lib.py <in.fa> <out.fa>
"""
import re
import sys
from collections import Counter

FIX = [
    (re.compile(r"\?+"), ""),
    (re.compile(r"^TIR/Tc1/Mariner$", re.I), "DNA/TcMar-Mariner"),
    (re.compile(r"^TIR/", re.I), "DNA/"),
    (re.compile(r"^(DNA|RC)/Helitron.*$", re.I), "RC/Helitron"),
    (re.compile(r"^Helitron$", re.I), "RC/Helitron"),
    (re.compile(r"^Simple/repeat$", re.I), "Simple_repeat"),
    (re.compile(r"^Low/complexity$", re.I), "Low_complexity"),
    (re.compile(r"^tandem_repeat.*$", re.I), "Satellite"),
    (re.compile(r"^DIRS/", re.I), "LTR/"),
    (re.compile(r"^DIRS$", re.I), "LTR/DIRS"),
    (re.compile(r"^(PLE|Penelope)(/.*)?$", re.I), "LINE/Penelope"),
    (re.compile(r"^Crypton(/.*)?$", re.I), "DNA/Crypton"),
    (re.compile(r"^RNA$", re.I), "rRNA"),
    (re.compile(r"^(LTR|LINE|SINE|DNA)$", re.I), lambda m: m.group(1).upper() + "/unknown"),
]


def fix(label):
    for rx, rep in FIX:
        label = rx.sub(rep, label)
    return label or "Unknown"


def main():
    src, dst = sys.argv[1], sys.argv[2]
    before, after, seen = Counter(), Counter(), Counter()
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                tok = line[1:].split()[0]
                name, _, lab = tok.partition("#")
                lab = lab or "Unknown"
                name = re.sub(r"[^\w.\-:+()]", "_", name)
                seen[name] += 1
                if seen[name] > 1:
                    name = f"{name}_dup{seen[name]}"
                new = fix(lab)
                before[lab] += 1
                after[new] += 1
                fo.write(f">{name}#{new}\n")
            elif line.strip():
                fo.write(line.strip().upper() + "\n")
    print(f"{src} -> {dst}: {sum(after.values())} sequences")
    changed = {k: fix(k) for k in before if fix(k) != k}
    if changed:
        print("  relabelled:", "; ".join(f"{k} -> {v} (x{before[k]})" for k, v in changed.items()))
    print("  labels:", after.most_common())


if __name__ == "__main__":
    main()

