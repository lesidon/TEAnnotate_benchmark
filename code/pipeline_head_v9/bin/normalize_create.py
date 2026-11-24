#!/usr/bin/env python3
"""Normalize CREATE predictions to a [base_id <TAB> RM_label] TSV for ensemble
voting. CREATE (yangqi-cs/CREATE) writes `relabeled_seqs.fasta` whose every header
is `>{predicted_label}\\t{original_description}`. We read the predicted label (the
text before the first TAB) and the original sequence id (first token after it),
map the label to a RepeatMasker Order/Superfamily, and emit one row per sequence.

Confirmed output header format (leaf name, a TAB, then our rewritten header):
  "> Copia\\tTE#rnd_1_family_1"
CREATE classifies in the Wicker taxonomy and, depending on --prob_thr, may stop at
an Order (e.g. "LTR", "TIR") rather than a Superfamily — so we accept full names,
Wicker 3-letter codes, and a keyword fallback. Unknown / too-broad / low-confidence
labels abstain (no vote).
"""
import argparse
import sys

UNKNOWN = {'', 'unknown', 'unclassified', 'unspecified', 'na', 'none', 'other'}

# Interior / too-broad tree nodes CREATE may emit when it stops early: no usable
# Order, so abstain rather than guess (esp. 'non-ltr', which must NOT match 'ltr').
NONINFORMATIVE = {'te', 'classi', 'class1', 'classii', 'class2',
                  'non-ltr', 'nonltr', 'sub1', 'sub2'}

# Exact Wicker superfamily / order names -> RepeatMasker label.
NAME_TO_RM = {
    'copia': 'LTR/Copia', 'gypsy': 'LTR/Gypsy', 'bel-pao': 'LTR/Pao', 'pao': 'LTR/Pao',
    'retrovirus': 'LTR/ERV', 'erv': 'LTR/ERV',
    'dirs': 'DIRS', 'ngaro': 'DIRS', 'viper': 'DIRS',
    'penelope': 'PLE/Penelope', 'ple': 'PLE',
    'l1': 'LINE/L1', 'l1_l2': 'LINE/L1', 'rte': 'LINE/RTE', 'jockey': 'LINE/Jockey',
    'r2': 'LINE/R2', 'i': 'LINE/I', 'line': 'LINE',
    'trna': 'SINE/tRNA', '7sl': 'SINE/7SL', '5s': 'SINE/5S', 'sine': 'SINE',
    'tc1-mariner': 'DNA/TcMar', 'mariner': 'DNA/TcMar', 'tcmar': 'DNA/TcMar',
    'hat': 'DNA/hAT', 'mutator': 'DNA/MULE-MuDR', 'mule': 'DNA/MULE-MuDR',
    'pif-harbinger': 'DNA/PIF-Harbinger', 'harbinger': 'DNA/PIF-Harbinger',
    'cacta': 'DNA/CMC-EnSpm', 'enspm': 'DNA/CMC-EnSpm', 'merlin': 'DNA/Merlin',
    'transib': 'DNA/CMC-Transib', 'p': 'DNA/P', 'piggybac': 'DNA/PiggyBac',
    'crypton': 'DNA/Crypton', 'helitron': 'RC/Helitron',
    'maverick': 'DNA/Maverick', 'polinton': 'DNA/Maverick',
    'tir': 'DNA', 'dna': 'DNA', 'ltr': 'LTR',
}

# Wicker 3-letter codes -> RepeatMasker label.
CODE_TO_RM = {
    'RLC': 'LTR/Copia', 'RLG': 'LTR/Gypsy', 'RLB': 'LTR/Pao', 'RLR': 'LTR/ERV',
    'RLE': 'LTR/ERV', 'RYD': 'DIRS', 'RYN': 'DIRS', 'RYV': 'DIRS',
    'RPP': 'PLE/Penelope', 'RIL': 'LINE/L1', 'RIT': 'LINE/RTE', 'RIJ': 'LINE/Jockey',
    'RIR': 'LINE/R2', 'RII': 'LINE/I', 'RST': 'SINE/tRNA', 'RSL': 'SINE/7SL',
    'RSS': 'SINE/5S', 'DTT': 'DNA/TcMar', 'DTA': 'DNA/hAT', 'DTM': 'DNA/MULE-MuDR',
    'DTP': 'DNA/PIF-Harbinger', 'DTH': 'DNA/PIF-Harbinger', 'DTC': 'DNA/CMC-EnSpm',
    'DTE': 'DNA/Merlin', 'DTR': 'DNA/CMC-Transib', 'DTB': 'DNA/PiggyBac',
    'DYC': 'DNA/Crypton', 'DHH': 'RC/Helitron', 'DMM': 'DNA/Maverick',
}
# Order from a Wicker code's first two letters (when the full code is unknown).
CODE_PREFIX_TO_RM = {
    'RL': 'LTR', 'RY': 'DIRS', 'RP': 'PLE', 'RI': 'LINE', 'RS': 'SINE',
    'DT': 'DNA', 'DY': 'DNA/Crypton', 'DH': 'RC/Helitron', 'DM': 'DNA/Maverick',
}
# Substring keyword fallback (ordered: most specific first).
KEYWORDS = [
    ('helitron', 'RC/Helitron'), ('maverick', 'DNA/Maverick'),
    ('polinton', 'DNA/Maverick'), ('crypton', 'DNA/Crypton'),
    ('copia', 'LTR/Copia'), ('gypsy', 'LTR/Gypsy'), ('pao', 'LTR/Pao'),
    ('erv', 'LTR/ERV'), ('dirs', 'DIRS'), ('penelope', 'PLE/Penelope'),
    ('mariner', 'DNA/TcMar'), ('tcmar', 'DNA/TcMar'), ('hat', 'DNA/hAT'),
    ('mutator', 'DNA/MULE-MuDR'), ('mule', 'DNA/MULE-MuDR'),
    ('harbinger', 'DNA/PIF-Harbinger'), ('cacta', 'DNA/CMC-EnSpm'),
    ('enspm', 'DNA/CMC-EnSpm'), ('merlin', 'DNA/Merlin'),
    ('sine', 'SINE'), ('line', 'LINE'), ('ltr', 'LTR'),
    ('tir', 'DNA'), ('dna', 'DNA'),
]


def label_to_rm(label):
    raw = (label or '').strip()
    l = raw.lower()
    if not l or l in UNKNOWN or l in NONINFORMATIVE:
        return None
    # Exact full-name match (handles "Copia", "LTR/Copia", "LTR", ...).
    leaf = l.split('/')[-1]
    if leaf in NAME_TO_RM:
        return NAME_TO_RM[leaf]
    if l in NAME_TO_RM:
        return NAME_TO_RM[l]
    # Wicker 3-letter code (case-sensitive uppercase).
    code = raw.split('/')[-1].upper()
    if code in CODE_TO_RM:
        return CODE_TO_RM[code]
    if len(code) >= 2 and code[:2] in CODE_PREFIX_TO_RM:
        return CODE_PREFIX_TO_RM[code[:2]]
    # Keyword fallback.
    for key, rm in KEYWORDS:
        if key in l:
            return rm
    return None


def iter_headers(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                yield line[1:].rstrip('\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pred', required=True, help='CREATE relabeled_seqs.fasta')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    n = 0
    with open(args.out, 'w') as out:
        for header in iter_headers(args.pred):
            # ">{pred_label}\t{original_description}"
            if '\t' not in header:
                continue
            pred_label, orig = header.split('\t', 1)
            toks = orig.split()
            if not toks:
                continue
            # The CREATE module rewrites input headers to ">TE#<orig_id>" (CREATE
            # reads the TRUE label from the text BEFORE '#'), so the real family id
            # is the part AFTER the first '#'. Fall back to the whole token if the
            # placeholder scheme was not applied.
            first = toks[0]
            base = first.split('#', 1)[1] if '#' in first else first
            rm = label_to_rm(pred_label)
            if not base or not rm:
                continue
            out.write('{}\t{}\n'.format(base, rm))
            n += 1
    sys.stderr.write('normalize_create: {} votes\n'.format(n))


if __name__ == '__main__':
    main()
