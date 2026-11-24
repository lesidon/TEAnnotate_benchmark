#!/usr/bin/env python3
"""Extract only the 'Unknown' consensi from a FASTA library, so the
re-classification ensemble runs on — and its votes apply to — just the
unclassified fraction (e.g. `rnd_1_family_4#Unknown`). Already-classified consensi
keep their labels and are not re-classified.

A record is treated as Unknown when the label after the first '#' is empty or one
of {unknown*, unspecified*, unclassified*, na, none} (case-insensitive; the order
prefix is checked too, so 'Unknown/Unknown' counts), or when there is no '#'.
"""
import argparse
import sys

UNKNOWN_EXACT = {'', 'na', 'none'}
UNKNOWN_PREFIX = ('unknown', 'unspecified', 'unclassified')


def is_unknown(header):
    tok = header.split()[0] if header.split() else header
    if '#' not in tok:
        return True
    cls = tok.split('#', 1)[1].strip().lower()
    cls0 = cls.split('/', 1)[0]
    if cls in UNKNOWN_EXACT or cls0 in UNKNOWN_EXACT:
        return True
    return cls.startswith(UNKNOWN_PREFIX) or cls0.startswith(UNKNOWN_PREFIX)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    keep = False
    n = 0
    with open(args.inp) as fh, open(args.out, 'w') as out:
        for line in fh:
            if line.startswith('>'):
                keep = is_unknown(line[1:].rstrip('\n'))
                if keep:
                    n += 1
            if keep:
                out.write(line)
    sys.stderr.write('extract_unknown: {} Unknown records kept\n'.format(n))


if __name__ == '__main__':
    main()
