#!/usr/bin/env python3
"""Normalize DeepTE predictions to a [base_id <TAB> RM_label] TSV for ensemble
voting. DeepTE writes opt_DeepTE.txt as <id>\\t<DeepTE_label> (e.g. ClassI_LTR_Copia,
ClassII_DNA_Mutator_DTM, ClassIII_Helitron, ClassI_nLTR_LINE_L1, unknown). We map
each label to a RepeatMasker Order/Superfamily and emit one row per classified
sequence (unmapped/unknown labels are skipped, i.e. the model abstains).
"""
import argparse
import sys

UNKNOWN = {'', 'unknown', 'unspecified', 'unclassified', 'na', 'none'}
LTR_SF  = {'copia': 'LTR/Copia', 'gypsy': 'LTR/Gypsy', 'bel': 'LTR/Pao',
           'pao': 'LTR/Pao', 'erv': 'LTR/ERV'}
LINE_SF = {'l1': 'LINE/L1', 'rte': 'LINE/RTE', 'jockey': 'LINE/Jockey',
           'cr1': 'LINE/CR1', 'i_': 'LINE/I'}
DNA_SF  = {'hat': 'DNA/hAT', 'mariner': 'DNA/TcMar', 'tc1': 'DNA/TcMar',
           'mutator': 'DNA/MULE-MuDR', 'mule': 'DNA/MULE-MuDR',
           'harbinger': 'DNA/PIF-Harbinger', 'pif': 'DNA/PIF-Harbinger',
           'cacta': 'DNA/CMC-EnSpm', 'enspm': 'DNA/CMC-EnSpm',
           'piggybac': 'DNA/PiggyBac', 'merlin': 'DNA/Merlin',
           'transib': 'DNA/CMC-Transib', 'p_': 'DNA/P'}


def deepte_to_rm(label):
    l = label.strip().lower()
    if not l or l in UNKNOWN:
        return None
    if 'helitron' in l:
        return 'RC/Helitron'
    if 'maverick' in l or 'polinton' in l:
        return 'DNA/Maverick'
    if 'sine' in l:
        return 'SINE'
    if 'line' in l or 'nltr' in l:
        for k, v in LINE_SF.items():
            if k in l:
                return v
        return 'LINE'
    if 'ltr' in l:
        for k, v in LTR_SF.items():
            if k in l:
                return v
        return 'LTR'
    if 'classii' in l or l.startswith('dna') or '_dna' in l or 'tir' in l:
        for k, v in DNA_SF.items():
            if k in l:
                return v
        return 'DNA'
    if 'classi' in l:
        return 'LTR'
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pred', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    n = 0
    with open(args.pred) as fh, open(args.out, 'w') as out:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            rm = deepte_to_rm(parts[1])
            if not rm:
                continue
            base = parts[0].split('#', 1)[0]
            out.write('{}\t{}\n'.format(base, rm))
            n += 1
    sys.stderr.write('normalize_deepte: {} votes\n'.format(n))


if __name__ == '__main__':
    main()
