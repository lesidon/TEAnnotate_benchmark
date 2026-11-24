#!/usr/bin/env python3
"""Post-process TEtrimmer output before the class-aware merge (CD-HIT + refine_dedup.py).

Reads TEtrimmer's un-merged consensus library (TEtrimmer_consensus.fasta; its own
length-based merge is replaced by the pipeline's class-aware one) and:

  1. drops Need_check consensi (inline tier tag, or a legacy need_check folder);
  2. re-checks labels. TEtrimmer keeps the input label on every consensus it derives
     from an input, even when the consensus is rebuilt from a different element that
     dominated the BLAST hits. For example, a chimeric LTR-int carrying part of a
     high-copy LINE comes back as a full-length LINE still labelled LTR/Gypsy. A label
     is changed only on positive evidence: another input of a DIFFERENT class covers
     >= --adopt-cov of the output (same family by the 80-80-80 rule) and at least
     --adopt-margin more of it than the output's own input does. The output then takes
     that input's label. Plain trimming or extension leaves the own-input coverage
     high, so those labels are untouched.
  3. with --evidence (TE_EVIDENCE's table), re-injects unchanged every input that has
     no surviving consensus but carries evidence of being a TE (tiers in
     --reinject-tiers: 1 = TE domain, 2 = sequence support, T = user library), and
     every input the gate kept unassessed because its domain search failed (action
     keep_unassessed): missing evidence is not negative evidence here either.
     TEtrimmer writes a skipped input back only when it has >= 2 full-length copies,
     so real low-copy and degraded families (old DNA transposons, fragmented
     Helitrons and LINEs) are otherwise lost; inputs whose consensi were all
     Need_check are recovered the same way.
  4. with --lineage and --out-user, the user's records (origin "user") are written to
     --out-user exactly as given, and consensi TEtrimmer rebuilt from them are dropped:
     the user's curated record is kept, not re-curated. The caller merges them back as
     user records, so the curated merge protects them as in Stage 2.
  5. writes clean headers (>name#label), without TEtrimmer's tab-separated tier fields.

Writes the library, a TSV report (one row per consensus) and, with --out-lineage, the
lineage of the output (each consensus with the input it came from as parent).
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

from refine_dedup import (base_of, blast_pairs, covered_bp, group_of, label_of,
                          read_fasta, read_lineage, write_lineage, wrap)


def read_consensus(path):
    """[(id, label, tier, seq)] from TEtrimmer's FASTA (header: name#label\\tTEtrimmer\\ttier)."""
    recs, head, buf = [], None, []

    def flush():
        if head is not None:
            fields = head.split('\t')
            tok = fields[0].split()[0]
            tier = fields[2].strip() if len(fields) > 2 else ''
            recs.append((base_of(tok), label_of(tok), tier, ''.join(buf)))

    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                flush()
                head, buf = line[1:], []
            else:
                buf.append(line.strip())
    flush()
    return recs


def legacy_need_check_ids(tdir):
    ids = set()
    for root, dirs, files in os.walk(tdir):
        if os.path.basename(root).lower() == 'need_check':
            ids.update(os.path.splitext(f)[0] for f in files)
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tetrimmer-dir', required=True)
    ap.add_argument('--input', required=True, help='library TEtrimmer was given')
    ap.add_argument('--out', required=True)
    ap.add_argument('--report', required=True)
    ap.add_argument('--adopt-cov', type=float, default=0.8)
    ap.add_argument('--adopt-margin', type=float, default=0.2)
    ap.add_argument('--min-identity', type=float, default=80.0)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--evidence', default=None, help="TE_EVIDENCE table (columns 'name', 'tier')")
    ap.add_argument('--reinject-tiers', default='1,2,T',
                    help='evidence tiers whose inputs return when TEtrimmer left no consensus')
    ap.add_argument('--lineage', default=None,
                    help='lineage of --input (id origin parent source_id)')
    ap.add_argument('--out-user', default=None,
                    help="with --lineage: write the user's records here, exactly as given")
    ap.add_argument('--out-lineage', default=None,
                    help='lineage of --out (and --out-user): each consensus with its input as parent')
    args = ap.parse_args()

    tdir = args.tetrimmer_dir
    cons_path = next((p for p in (os.path.join(tdir, 'TEtrimmer_consensus.fasta'),
                                  os.path.join(tdir, 'TEtrimmer_consensus_merged.fasta'))
                      if os.path.exists(p) and os.path.getsize(p) > 0), None)
    if cons_path is None:
        sys.stderr.write('tetrimmer_postprocess: no TEtrimmer consensus\n')
        recs = []
    else:
        recs = read_consensus(cons_path)

    # TEtrimmer renames inputs (Sequence_name_mapping.txt) and reports input -> output
    # per consensus in summary.txt; together they give each consensus its original input.
    inputs, _ = read_fasta(args.input)
    modified_to_orig = {}
    mp = os.path.join(tdir, 'Sequence_name_mapping.txt')
    if os.path.exists(mp):
        with open(mp) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 2:
                    modified_to_orig[f[1].strip()] = f[0].split()[0]
    out_to_in = {}
    sp = os.path.join(tdir, 'summary.txt')
    if os.path.exists(sp):
        with open(sp) as fh:
            for row in csv.DictReader(fh):
                if row.get('output_name') and row.get('input_name'):
                    out_to_in[row['output_name']] = row['input_name']

    def origin_input(name):
        """The input (header token) a TEtrimmer consensus was built from, or None."""
        return modified_to_orig.get(out_to_in.get(name, name))

    lineage = read_lineage(args.lineage)
    user = {i for i, r in lineage.items() if r['origin'] == 'user'} if args.out_user else set()

    legacy = legacy_need_check_ids(tdir)
    kept, report = [], []
    for name, label, tier, seq in recs:
        if 'need_check' in tier.lower().replace(' ', '_') or name in legacy:
            report.append((name, '', tier, 'drop_need_check', label, '', ''))
            continue
        orig = origin_input(name)
        if orig and base_of(orig) in user:
            report.append((name, orig, tier, 'drop_user_rebuilt', label, '',
                           "rebuilt from a user record, which is kept as given"))
            continue
        kept.append([name, label, tier, seq])

    # label re-check (classified consensi with a known input only)
    oq = {f'o{i}': r[3] for i, r in enumerate(kept)
          if group_of(r[1]) and modified_to_orig.get(out_to_in.get(r[0], r[0])) in inputs}
    inames = list(inputs)
    isub = {f'i{j}': inputs[n] for j, n in enumerate(inames)}
    spans_o, spans_i = defaultdict(list), defaultdict(list)
    if oq:
        for q, s, pid, ln, qsp, ssp in blast_pairs(oq, isub, args.threads, args.min_identity):
            spans_o[(q, s)].append(qsp)
            spans_i[(q, s)].append(ssp)

    counts = defaultdict(int)
    for i, r in enumerate(kept):
        name, label, tier, seq = r
        key = f'o{i}'
        orig = modified_to_orig.get(out_to_in.get(name, name))
        if key not in oq:
            report.append((name, orig or '', tier, 'keep', label, label, 'not checked'))
            continue
        ikey = f'i{inames.index(orig)}'
        lo, li = len(seq), len(inputs[orig])
        self_cov = covered_bp(lo, spans_o.get((key, ikey), [])) / lo
        in_cov = covered_bp(li, spans_i.get((key, ikey), [])) / li
        detail = f'self_cov {self_cov:.2f}, input_cov {in_cov:.2f}'
        best, best_cov = None, 0.0
        for (q, s), sps in spans_o.items():
            if q != key or s == ikey:
                continue
            other = inames[int(s[1:])]
            g = group_of(label_of(other))
            if not g or g == group_of(label):
                continue
            c = covered_bp(lo, sps) / lo
            if c > best_cov:
                best, best_cov = other, c
        if best and best_cov >= args.adopt_cov and best_cov >= self_cov + args.adopt_margin:
            r[1] = label_of(best)
            counts['relabel_from_input'] += 1
            report.append((name, orig, tier, 'relabel_from_input', label, r[1],
                           detail + f'; {best_cov:.2f} covered by {best}'))
        else:
            if best:
                detail += f'; best other-class input {best} ({best_cov:.2f})'
            report.append((name, orig, tier, 'keep', label, label, detail))

    # Inputs with evidence but no surviving consensus come back unchanged. Survivors
    # are traced to their inputs through TEtrimmer's renaming table; without it they
    # cannot be, and re-injection would duplicate curated families.
    if args.evidence and recs and not modified_to_orig:
        sys.stderr.write('tetrimmer_postprocess: WARNING no Sequence_name_mapping.txt; '
                         'dropped inputs are not re-injected\n')
    elif args.evidence:
        tiers, actions = {}, {}
        with open(args.evidence) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                tiers[row['name']] = row['tier']
                actions[row['name']] = row.get('action', '')
        wanted = {t.strip() for t in args.reinject_tiers.split(',') if t.strip()}
        survivors = {origin_input(r[0]) or r[0] for r in kept}
        survivors |= {base_of(s) for s in survivors}
        for orig in inputs:
            base = base_of(orig)
            if orig in survivors or base in survivors or base in user:
                continue
            tier = tiers.get(base)
            unassessed = actions.get(base) == 'keep_unassessed'
            if tier in wanted or unassessed:
                kept.append([base, label_of(orig), '', inputs[orig]])
                counts['reinject'] += 1
                why = f'evidence tier {tier}'
                if tier not in wanted:
                    why = f'kept unassessed by the gate (domain search failed; tier {tier})'
                    counts['reinject_unassessed'] += 1
                report.append((base, orig, '', 'reinject', label_of(orig), label_of(orig),
                               f'no surviving TEtrimmer consensus; {why}'))
        if counts['reinject_unassessed']:
            sys.stderr.write(f"tetrimmer_postprocess: WARNING {counts['reinject_unassessed']} inputs "
                             'without domain assessment re-injected (the domain search failed)\n')

    with open(args.out, 'w') as out:
        for name, label, tier, seq in kept:
            out.write(f'>{name}#{label}\n{wrap(seq)}\n')
    out_lineage = {}
    for name, label, tier, seq in kept:
        orig = origin_input(name)
        src = base_of(orig) if orig else name
        row = lineage.get(src, {'origin': 'denovo', 'source_id': src})
        out_lineage[name] = {'origin': row['origin'], 'source_id': row['source_id'],
                             'parent': src if src != name else '-'}
    if args.out_user:
        with open(args.out_user, 'w') as out:
            for orig, seq in inputs.items():
                b = base_of(orig)
                if b in user:
                    out.write(f'>{orig}\n{wrap(seq)}\n')
                    out_lineage[b] = {'origin': 'user', 'parent': '-',
                                      'source_id': lineage[b]['source_id']}
                    report.append((b, orig, '', 'keep_user', label_of(orig), label_of(orig),
                                   'user record, kept as given'))
                    counts['keep_user'] += 1
    if args.out_lineage:
        write_lineage(args.out_lineage, out_lineage)
    with open(args.report, 'w') as fh:
        fh.write('consensus\tinput\ttier\taction\tlabel_in\tlabel_out\tdetail\n')
        for row in report:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    n_drop = sum(1 for r in report if r[3] == 'drop_need_check')
    sys.stderr.write(
        f'tetrimmer_postprocess: {len(recs)} consensi from '
        f'{os.path.basename(cons_path) if cons_path else "no TEtrimmer consensus file"}; '
        f'dropped {n_drop} Need_check; relabelled {counts["relabel_from_input"]}; '
        f're-injected {counts["reinject"]} ({counts["reinject_unassessed"]} unassessed); '
        f'kept {len(kept)}; user records {counts["keep_user"]}\n')


if __name__ == '__main__':
    main()
