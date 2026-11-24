#!/usr/bin/env python3
"""Normalize TEclass2 predictions to a [base_id <TAB> RM_label] TSV for ensemble
voting. TEclass2 output is whitespace/tab-delimited: <id> <Wicker_superfamily>
<scores...>. We map the superfamily to a RepeatMasker Order/Superfamily and emit
one row per classified sequence (header rows and unmapped labels are skipped).
"""
import argparse
import sys

WICKER_TO_RM = {
    'Copia': 'LTR/Copia', 'Gypsy': 'LTR/Gypsy', 'ERV': 'LTR/ERV', 'Pao': 'LTR/Pao',
    'L1_L2': 'LINE/L1', 'Jockey': 'LINE/Jockey', 'RTE': 'LINE/RTE', 'SINE': 'SINE',
    'hAT': 'DNA/hAT', 'TcMar': 'DNA/TcMar', 'Merlin': 'DNA/Merlin',
    'Transib': 'DNA/CMC-Transib', 'P': 'DNA/P', 'PIF': 'DNA/PIF-Harbinger',
    'Maverick': 'DNA/Maverick', 'Crypton': 'DNA/Crypton',
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pred', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-score', type=float, default=0.0,
                    help='abstain (no vote) below this top softmax score')
    args = ap.parse_args()

    n = 0
    with open(args.pred) as fh, open(args.out, 'w') as out:
        for line in fh:
            f = line.split()
            # Skip blank/header rows: a data row's 2nd field is a known Wicker
            # superfamily (the header row's is not).
            if len(f) < 2 or f[1] not in WICKER_TO_RM:
                continue
            scores = []
            for tok in f[2:]:
                try:
                    scores.append(float(tok))
                except ValueError:
                    pass
            if scores and max(scores) < args.min_score:
                continue
            base = f[0].split('#', 1)[0]
            out.write('{}\t{}\n'.format(base, WICKER_TO_RM[f[1]]))
            n += 1
    sys.stderr.write('normalize_teclass2: {} votes\n'.format(n))


if __name__ == '__main__':
    main()
