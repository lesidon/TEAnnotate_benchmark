#!/usr/bin/env python3
"""Normalize Terrier predictions to a [base_id <TAB> RM_label] TSV for ensemble
voting. Terrier (rbturnbull/terrier, `bio-terrier`) writes --output-csv with, per
sequence, a direct call plus every class probability. Confirmed columns:

    file,accession,prediction,probability,original_id,original_classification,
    description,LTR,SINE,DNA,...,LTR/Copia,LTR/Gypsy,...,DNA/hAT,...,LINE/L1,...

We use the `prediction` column (Terrier's own arg-max — an Order like "LTR" or a
superfamily like "LTR/Copia"), keyed by `accession`/`original_id`, and map it to a
RepeatMasker Order/Superfamily; `probability` gates it via --min-prob. If no
`prediction` column is present (older layouts), we arg-max the per-class probability
columns instead. Unknown / low-probability / unmappable predictions abstain.
"""
import argparse
import csv
import sys

UNKNOWN = {'', 'unknown', 'unclassified', 'unspecified', 'na', 'none', 'other',
           'satellite', 'structural_rna', 'simple_repeat', 'tandem'}

# Columns that are metadata, not a per-class probability (excluded from arg-max).
NAME_HINTS = ('accession', 'original_id', 'seqname', 'seq_name', 'sequence',
              'name', 'id', 'query')
PRED_HINTS = ('prediction', 'predicted', 'pred', 'classification', 'class', 'label')
PROB_HINTS = ('probability', 'prob', 'score', 'confidence')
META_COLS = set(NAME_HINTS) | set(PRED_HINTS) | set(PROB_HINTS) | {
    'file', 'description', 'original_classification'}

EXACT = {
    'ltr/copia': 'LTR/Copia', 'ltr/gypsy': 'LTR/Gypsy', 'ltr/bel-pao': 'LTR/Pao',
    'ltr/pao': 'LTR/Pao', 'ltr/erv': 'LTR/ERV', 'ltr/dirs': 'DIRS', 'dirs': 'DIRS',
    'line/l1': 'LINE/L1', 'line/rte': 'LINE/RTE', 'line/jockey': 'LINE/Jockey',
    'line/cr1': 'LINE/CR1', 'line/r2': 'LINE/R2', 'line/i': 'LINE/I',
    'line/l2': 'LINE/L2', 'line/r1': 'LINE/R1',
    'sine/trna': 'SINE/tRNA', 'sine/7sl': 'SINE/7SL', 'sine/5s': 'SINE/5S',
    'dna/tcmar': 'DNA/TcMar', 'dna/hat': 'DNA/hAT', 'dna/mule': 'DNA/MULE-MuDR',
    'dna/mule-mudr': 'DNA/MULE-MuDR', 'dna/pif-harbinger': 'DNA/PIF-Harbinger',
    'dna/harbinger': 'DNA/PIF-Harbinger', 'dna/cmc': 'DNA/CMC-EnSpm',
    'dna/cmc-enspm': 'DNA/CMC-EnSpm', 'dna/maverick': 'DNA/Maverick',
    'dna/crypton': 'DNA/Crypton', 'dna/piggybac': 'DNA/PiggyBac',
    'dna/merlin': 'DNA/Merlin', 'dna/p': 'DNA/P',
    'rc/helitron': 'RC/Helitron', 'rc': 'RC/Helitron',
    'ple/penelope': 'PLE/Penelope', 'ple': 'PLE',
}

# Order-level fallbacks and substring keywords (most specific first).
KEYWORDS = [
    ('helitron', 'RC/Helitron'), ('maverick', 'DNA/Maverick'),
    ('polinton', 'DNA/Maverick'), ('crypton', 'DNA/Crypton'),
    ('copia', 'LTR/Copia'), ('gypsy', 'LTR/Gypsy'), ('bel', 'LTR/Pao'),
    ('pao', 'LTR/Pao'), ('erv', 'LTR/ERV'), ('dirs', 'DIRS'),
    ('penelope', 'PLE/Penelope'), ('ple', 'PLE'),
    ('cacta', 'DNA/CMC-EnSpm'), ('enspm', 'DNA/CMC-EnSpm'),
    ('harbinger', 'DNA/PIF-Harbinger'), ('pif', 'DNA/PIF-Harbinger'),
    ('mariner', 'DNA/TcMar'), ('tcmar', 'DNA/TcMar'), ('tc1', 'DNA/TcMar'),
    ('mudr', 'DNA/MULE-MuDR'), ('mule', 'DNA/MULE-MuDR'), ('mutator', 'DNA/MULE-MuDR'),
    ('hat', 'DNA/hAT'), ('merlin', 'DNA/Merlin'), ('piggybac', 'DNA/PiggyBac'),
    ('l1', 'LINE/L1'), ('rte', 'LINE/RTE'), ('jockey', 'LINE/Jockey'),
    ('cr1', 'LINE/CR1'), ('r2', 'LINE/R2'),
    ('non-ltr', 'LINE'), ('nonltr', 'LINE'),
    ('sine', 'SINE'), ('line', 'LINE'), ('ltr', 'LTR'),
    ('tir', 'DNA'), ('dna', 'DNA'),
]


def class_to_rm(label):
    l = (label or '').strip().lower()
    if not l or l in UNKNOWN:
        return None
    if l in EXACT:
        return EXACT[l]
    for key, rm in KEYWORDS:
        if key in l:
            return rm
    return None


def is_number(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def find_col(header, hints):
    low = [f.strip().lower() for f in header]
    for hint in hints:                       # exact match first
        if hint in low:
            return low.index(hint)
    for hint in hints:                       # then substring
        for i, f in enumerate(low):
            if hint in f:
                return i
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pred', required=True, help='Terrier --output-csv file')
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-prob', type=float, default=0.0,
                    help='abstain if the prediction/arg-max probability is below this')
    args = ap.parse_args()

    # Counters so an empty result is never silent — we always say why.
    stats = {'rows': 0, 'votes': 0, 'no_base': 0, 'low_prob': 0,
             'unmapped': 0, 'bad_row': 0}
    unmapped_examples = []

    with open(args.pred, newline='', encoding='utf-8-sig') as fh:  # tolerate BOM
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            header = []
        low = [f.strip().lower() for f in header]
        name_i = find_col(header, NAME_HINTS)
        if name_i is None:
            name_i = 0
        pred_i = find_col(header, PRED_HINTS)
        prob_i = find_col(header, PROB_HINTS)
        class_cols = [i for i in range(len(header)) if low[i] not in META_COLS]

        if not header:
            sys.stderr.write('normalize_terrier: ERROR empty CSV (no header) '
                             '-> 0 votes\n')
        sys.stderr.write(
            'normalize_terrier: columns name=%s pred=%s prob=%s (%d class cols)\n'
            % (header[name_i] if name_i < len(header) else '?',
               header[pred_i] if pred_i is not None else 'NONE(argmax)',
               header[prob_i] if prob_i is not None else 'NONE',
               len(class_cols)))

        with open(args.out, 'w') as out:
            for row in reader:
                stats['rows'] += 1
                try:
                    if not row or len(row) <= name_i:
                        stats['bad_row'] += 1
                        continue
                    base = row[name_i].split('#', 1)[0].strip().lstrip('>')
                    if not base:
                        stats['no_base'] += 1
                        continue

                    if pred_i is not None and len(row) > pred_i and row[pred_i].strip():
                        label = row[pred_i]
                        if (prob_i is not None and len(row) > prob_i
                                and is_number(row[prob_i])
                                and float(row[prob_i]) < args.min_prob):
                            stats['low_prob'] += 1
                            continue
                    else:
                        best_p, label = -1.0, None
                        for i in class_cols:
                            if i < len(row) and is_number(row[i]):
                                p = float(row[i])
                                if p > best_p:
                                    best_p, label = p, header[i]
                        if label is None or best_p < args.min_prob:
                            stats['low_prob'] += 1
                            continue

                    rm = class_to_rm(label)
                    if not rm:
                        stats['unmapped'] += 1
                        if len(unmapped_examples) < 5:
                            unmapped_examples.append(label)
                        continue
                    out.write('{}\t{}\n'.format(base, rm))
                    stats['votes'] += 1
                except Exception as exc:                 # never abort the task
                    stats['bad_row'] += 1
                    if stats['bad_row'] <= 3:
                        sys.stderr.write(
                            'normalize_terrier: skipped a malformed row (%s)\n' % exc)

    sys.stderr.write(
        'normalize_terrier: %d votes from %d rows '
        '(skipped: %d low-prob, %d unmapped, %d no-id, %d malformed)\n'
        % (stats['votes'], stats['rows'], stats['low_prob'], stats['unmapped'],
           stats['no_base'], stats['bad_row']))
    if stats['rows'] == 0:
        sys.stderr.write(
            'normalize_terrier: WARNING the CSV had a header but NO data rows — '
            'Terrier likely crashed mid-run (e.g. GPU/driver). Terrier abstains.\n')
    elif stats['votes'] == 0 and stats['unmapped']:
        sys.stderr.write(
            'normalize_terrier: WARNING no votes; %d predictions did not map to an '
            'RM label, e.g. %s\n' % (stats['unmapped'], ', '.join(unmapped_examples)))


if __name__ == '__main__':
    main()
