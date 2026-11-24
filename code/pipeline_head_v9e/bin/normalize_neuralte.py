#!/usr/bin/env python3
"""Normalize NeuralTE predictions to a [base_id <TAB> RM_label] TSV for ensemble
voting. NeuralTE (CSU-KangHu/NeuralTE) run with --is_predict 1 --is_wicker 0 writes
`classified.info`, a 2-column CSV:

    #Seq_Name,Predict_Label
    rnd_1_family_1,LTR_Copia
    rnd_1_family_2,Unknown

With --is_wicker 0 the labels are DeepTE-style underscore strings (e.g. LTR_Copia,
DNA_nMITE_hAT, RC_Helitron). We map each to a RepeatMasker Order/Superfamily and
emit one row per confidently-classified sequence; Unknown / unmappable labels
abstain (no vote).
"""
import argparse
import sys

UNKNOWN = {'', 'unknown', 'unclassified', 'unspecified', 'na', 'none', 'other'}

# Exact NeuralTE (is_wicker 0) labels -> RepeatMasker Order/Superfamily.
EXACT = {
    'ltr_copia': 'LTR/Copia', 'ltr_gypsy': 'LTR/Gypsy', 'ltr_erv': 'LTR/ERV',
    'ltr_bel': 'LTR/Pao', 'ltr_belpao': 'LTR/Pao', 'dirs_dirs': 'DIRS',
    'rc_helitron': 'RC/Helitron',
    'dna_nmite_tc': 'DNA/TcMar', 'dna_nmite_mariner': 'DNA/TcMar',
    'dna_nmite_harbinger': 'DNA/PIF-Harbinger', 'dna_nmite_hat': 'DNA/hAT',
    'dna_nmite_cacta': 'DNA/CMC-EnSpm', 'dna_nmite_mudr': 'DNA/MULE-MuDR',
    'dna_nmite_mutator': 'DNA/MULE-MuDR',
    # MITE variants collapse to the same superfamily as their nMITE parent.
    'dna_mite_tc': 'DNA/TcMar', 'dna_mite_harbinger': 'DNA/PIF-Harbinger',
    'dna_mite_hat': 'DNA/hAT', 'dna_mite_cacta': 'DNA/CMC-EnSpm',
    'dna_mite_mudr': 'DNA/MULE-MuDR',
}

# Substring keyword fallback (most specific first) — robust to label-style drift.
KEYWORDS = [
    ('helitron', 'RC/Helitron'), ('maverick', 'DNA/Maverick'),
    ('polinton', 'DNA/Maverick'), ('crypton', 'DNA/Crypton'),
    ('copia', 'LTR/Copia'), ('gypsy', 'LTR/Gypsy'),
    ('bel', 'LTR/Pao'), ('pao', 'LTR/Pao'), ('erv', 'LTR/ERV'),
    ('dirs', 'DIRS'), ('penelope', 'PLE/Penelope'), ('ple', 'PLE'),
    ('cacta', 'DNA/CMC-EnSpm'), ('enspm', 'DNA/CMC-EnSpm'),
    ('harbinger', 'DNA/PIF-Harbinger'), ('pif', 'DNA/PIF-Harbinger'),
    ('mariner', 'DNA/TcMar'), ('tcmar', 'DNA/TcMar'), ('_tc', 'DNA/TcMar'),
    ('mudr', 'DNA/MULE-MuDR'), ('mutator', 'DNA/MULE-MuDR'), ('mule', 'DNA/MULE-MuDR'),
    ('hat', 'DNA/hAT'), ('merlin', 'DNA/Merlin'), ('piggybac', 'DNA/PiggyBac'),
    ('l1', 'LINE/L1'), ('rte', 'LINE/RTE'), ('jockey', 'LINE/Jockey'),
    ('cr1', 'LINE/CR1'), ('r2', 'LINE/R2'),
    # Orders last (only if no superfamily keyword matched). 'non-ltr'/'nltr' are
    # Non-LTR retro = LINE-ish; handled before the bare 'ltr' catch.
    ('non_ltr', 'LINE'), ('non-ltr', 'LINE'), ('nltr', 'LINE'),
    ('sine', 'SINE'), ('line', 'LINE'), ('ltr', 'LTR'),
    ('tir', 'DNA'), ('dna', 'DNA'),
]


def label_to_rm(label):
    l = (label or '').strip().lower()
    if not l or l in UNKNOWN:
        return None
    if l in EXACT:
        return EXACT[l]
    for key, rm in KEYWORDS:
        if key in l:
            return rm
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pred', required=True, help='NeuralTE classified.info')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    n = 0
    with open(args.pred) as fh, open(args.out, 'w') as out:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split(',', 1)
            if len(parts) < 2:
                continue
            base = parts[0].split('#', 1)[0].strip()
            rm = label_to_rm(parts[1])
            if not base or not rm:
                continue
            out.write('{}\t{}\n'.format(base, rm))
            n += 1
    sys.stderr.write('normalize_neuralte: {} votes\n'.format(n))


if __name__ == '__main__':
    main()
