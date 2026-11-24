#!/usr/bin/env python3
"""Evidence table and evidence gate for a deduplicated TE consensus library.

Every consensus gets one evidence tier, from the strongest evidence it has:

  0  non-TE label (Satellite, rRNA, tRNA, snRNA, Simple_repeat, Low_complexity, ...)
  1  TE protein domain: a TEsorter / REXdb classification, or a TE-specific Pfam
     domain on any six-frame ORF fragment (--pfam-hits; see te_pfam.py). REXdb's
     clade-level profiles miss divergent TEs of unexplored lineages that Pfam's broad
     family models still detect.
  2  another discovery tool's library contains the same element: hits from that
     tool's sequences cover >= --support-cov of the consensus at >= --support-id %
  3  copy evidence: >= --min-full-copies full-length genomic copies (>= 75% identity
     over >= 90% of the consensus; TEtrimmer's own definition), or, when
     --min-half-copies > 0, that many copies covering >= 50% of it. Old families decay
     by mutation and deletion and keep partial copies only (RIP makes this extreme in
     fungi); in unexplored genomes they also lack a known domain and a second tool.
  4  none of the above
Sequences of a --trusted library (the user's curated families) get tier "T" and
are always kept.

Gate (--mode):
  lite    drop tiers 0 and 4 (default)
  strict  drop tiers 0, 3 and 4
  off     keep everything (the table is still written)
Tier-0 consensi go to --out-nonte (a repeat library, but not a TE one); tier-4 (and,
in strict mode, tier-3) consensi go to --out-removed.

With --relabel-unknown, an Unknown consensus with a TEsorter domain class takes that
class (RepeatMasker vocabulary); "mixture" calls, i.e. domains of several classes,
are not used. With --relabel-conflicts, a classified consensus whose TEsorter class
(LTR / LINE-SINE / TIR / Helitron) differs from its label also takes the domain's
class. Meant for curated consensi (after TEtrimmer), where one consensus is one
element; on uncurated, possibly chimeric consensi a domain describes only a part.

The source tool of each consensus comes from --provenance (name<TAB>tool, written when
the tool libraries are concatenated) or, failing that, from the naming scheme of each
tool. Pieces cut out of a chimera (<name>_partN) inherit the tool of <name>.
"""
import argparse
import re
import sys
from collections import defaultdict

NON_TE = re.compile(r'^(satellite|simple_repeat|simple|low_complexity|rrna|trna|snrna|scrna|'
                    r'srprna|cent|centr|rdna|artefact|artifact|segmental|buffer)', re.I)
UNKNOWN = re.compile(r'^(unknown|unclassified|unspecified|repeat|repeat_region|na|none|'
                     r'other)?$', re.I)
PART = re.compile(r'_part\d+$')

# TEsorter Order/Superfamily -> RepeatMasker vocabulary
TESORTER_TO_RM = {
    'LTR/Copia': 'LTR/Copia', 'LTR/Gypsy': 'LTR/Gypsy', 'LTR/Bel-Pao': 'LTR/Pao',
    'LTR/Retrovirus': 'LTR/ERV', 'LTR/unknown': 'LTR', 'LTR/mixture': 'LTR',
    'pararetrovirus/unknown': 'LTR/Caulimovirus', 'DIRS/unknown': 'LTR/DIRS',
    'LINE/unknown': 'LINE', 'SINE/unknown': 'SINE', 'Penelope/unknown': 'LINE/Penelope',
    'TIR/EnSpm_CACTA': 'DNA/CMC-EnSpm', 'TIR/hAT': 'DNA/hAT', 'TIR/MuDR_Mutator': 'DNA/MULE-MuDR',
    'TIR/PIF_Harbinger': 'DNA/PIF-Harbinger', 'TIR/Tc1_Mariner': 'DNA/TcMar',
    'TIR/PiggyBac': 'DNA/PiggyBac', 'TIR/Merlin': 'DNA/Merlin', 'TIR/Transib': 'DNA/CMC-Transib',
    'TIR/P': 'DNA/P', 'TIR/Sola': 'DNA/Sola', 'TIR/Kolobok': 'DNA/Kolobok',
    'TIR/Academ': 'DNA/Academ', 'TIR/Ginger': 'DNA/Ginger', 'TIR/unknown': 'DNA',
    'Helitron/unknown': 'RC/Helitron', 'Maverick/unknown': 'DNA/Maverick',
    'Crypton/unknown': 'DNA/Crypton',
}


def tesorter_label(order, sf):
    """RepeatMasker label for a TEsorter call, or None for multi-class (mixture) calls."""
    if order == 'mixture' or (sf == 'mixture' and order != 'LTR'):
        return None
    key = f'{order}/{sf or "unknown"}'
    if key in TESORTER_TO_RM:
        return TESORTER_TO_RM[key]
    if order == 'LINE':
        return f'LINE/{sf}'
    if order == 'TIR':
        return f'DNA/{sf}'
    return TESORTER_TO_RM.get(f'{order}/unknown')


def class_group(label):
    """Coarse class of a RepeatMasker-style label: ltr / nonltr / tir / helitron / None."""
    l = label.lower()
    top = l.split('/')[0]
    if 'helitron' in l or top == 'rc':
        return 'helitron'
    if top in ('ltr', 'dirs', 'erv') or re.search(r'gypsy|copia|pao|retrovirus|caulimo', l):
        return 'ltr'
    if top in ('line', 'sine', 'ple', 'penelope', 'retroposon'):
        return 'nonltr'
    if top in ('dna', 'tir', 'mite'):
        return 'tir'
    return None


def tool_from_name(name):
    n = name.replace('_', '-')
    if re.match(r'TE-\d+', n):
        return 'edta'
    if re.match(r'(LTR|TIR|Helitron|Homology-Non-LTR|Denovo-Non-LTR)-\d+', n):
        return 'hite'
    if n.startswith('insertion-HLE'):
        return 'heliano'
    if n.startswith('rnd-'):
        return 'earlgrey'
    return 'unknown'


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def split_header(h):
    base, _, label = h.partition('#')
    return base, label


def covered(length, spans):
    if not spans:
        return 0
    spans = sorted(spans)
    total, (s, e) = 0, spans[0]
    for a, b in spans[1:]:
        if a > e + 1:
            total += e - s + 1
            s, e = a, b
        else:
            e = max(e, b)
    return min(length, total + e - s + 1)


def wrap(seq, width=80):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--library', required=True, help='deduplicated library (name#label)')
    ap.add_argument('--raw', required=True, help='raw combined library the support hits refer to')
    ap.add_argument('--provenance', default=None, help='name<TAB>tool for --raw sequences')
    ap.add_argument('--trusted', default=None, help='curated library whose sequences are never gated')
    ap.add_argument('--tesorter-cls', required=True, help='TEsorter *.cls.tsv (may be empty)')
    ap.add_argument('--pfam-hits', default=None,
                    help='hmmsearch --domtblout of the TE Pfam profiles vs ORF fragments '
                         'named <consensus>|<frame>|<k>')
    ap.add_argument('--pfam-evalue', type=float, default=1e-5,
                    help='independent E-value for a Pfam domain to count (default 1e-5)')
    ap.add_argument('--support-hits', required=True,
                    help='blastn outfmt "6 qseqid sseqid pident length qstart qend qlen", library vs raw')
    ap.add_argument('--genome-hits', required=True,
                    help='blastn outfmt "6 qseqid sseqid pident length qstart qend qlen", library vs '
                         'genome (HSPs covering >= 50%% of the query at >= 70%% identity)')
    ap.add_argument('--mode', choices=('lite', 'strict', 'off'), default='lite')
    ap.add_argument('--relabel-unknown', action='store_true')
    ap.add_argument('--relabel-conflicts', action='store_true',
                    help='a TEsorter class that contradicts the label replaces it (curated consensi)')
    ap.add_argument('--support-id', type=float, default=80.0)
    ap.add_argument('--support-cov', type=float, default=0.5)
    ap.add_argument('--min-full-copies', type=int, default=2)
    ap.add_argument('--min-half-copies', type=int, default=0,
                    help='half-length copies that count as copy evidence (0 = off)')
    ap.add_argument('--out-library', required=True)
    ap.add_argument('--out-nonte', required=True)
    ap.add_argument('--out-removed', required=True)
    ap.add_argument('--evidence', required=True, help='output TSV')
    args = ap.parse_args()

    lib = read_fasta(args.library)

    # source tool of every raw sequence
    tool_of = {}
    if args.provenance:
        with open(args.provenance) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 2 and f[0]:
                    tool_of[split_header(f[0])[0]] = f[1]
    for h, _ in read_fasta(args.raw):
        b = split_header(h)[0]
        tool_of.setdefault(b, tool_from_name(b))

    trusted = set()
    if args.trusted:
        for h, _ in read_fasta(args.trusted):
            b = split_header(h)[0]
            trusted.add(b)
            tool_of[b] = 'user'

    def tool(base):
        return tool_of.get(base) or tool_of.get(PART.sub('', base)) or tool_from_name(base)

    # TE protein domains (TEsorter keys sequences by their full FASTA ID)
    domain = {}
    try:
        with open(args.tesorter_cls) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) >= 3:
                    domain[split_header(f[0])[0]] = (f[1], f[2])
    except FileNotFoundError:
        pass

    # TE-specific Pfam domains on ORF fragments
    pfam = defaultdict(set)
    if args.pfam_hits:
        try:
            with open(args.pfam_hits) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    f = line.split()
                    if len(f) >= 13 and float(f[12]) <= args.pfam_evalue:
                        pfam[f[0].split('|')[0]].add(f[3])
        except FileNotFoundError:
            pass

    # cross-tool support: coverage of each consensus by sequences of OTHER tools
    spans = defaultdict(list)
    qlen = {}
    with open(args.support_hits) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 7:
                continue
            q, s = split_header(f[0])[0], split_header(f[1])[0]
            if float(f[2]) < args.support_id:
                continue
            ts = tool(s)
            if ts == tool(q) or ts == 'unknown':
                continue
            spans[(q, ts)].append((int(f[4]), int(f[5])))
            qlen[q] = int(f[6])
    support = defaultdict(set)
    for (q, ts), sp in spans.items():
        if covered(qlen[q], sp) >= args.support_cov * qlen[q]:
            support[q].add(ts)

    # genomic copies: full-length (>= 75% identity over >= 90%) and half-length (>= 50%;
    # blastn already keeps only HSPs at >= 70% identity covering >= 50% of the query)
    n_full, n_half = defaultdict(int), defaultdict(int)
    with open(args.genome_hits) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 7:
                continue
            q = split_header(f[0])[0]
            cov = (abs(int(f[5]) - int(f[4])) + 1) / int(f[6])
            if cov >= 0.5:
                n_half[q] += 1
            if cov >= 0.9 and float(f[2]) >= 75:
                n_full[q] += 1

    drop_tiers = {'lite': {4}, 'strict': {3, 4}, 'off': set()}[args.mode]
    counts = defaultdict(int)
    with open(args.out_library, 'w') as keep, open(args.out_nonte, 'w') as nonte, \
         open(args.out_removed, 'w') as removed, open(args.evidence, 'w') as ev:
        ev.write('name\ttool\tlabel_in\tlabel_out\tlength\tdomain\tsupport_tools\t'
                 'full_copies\ttier\taction\thalf_copies\tpfam\n')
        for h, seq in lib:
            base, label = split_header(h)
            d = domain.get(base)
            dom_label = tesorter_label(*d) if d else None
            pf = sorted(pfam.get(base, ()))
            sup = sorted(support.get(base, ()))
            nf = n_full.get(base, 0)
            nh = n_half.get(base, 0)
            if base in trusted or PART.sub('', base) in trusted:
                tier = 'T'
            elif NON_TE.match(label):
                tier = 0
            elif d or pf:
                tier = 1
            elif sup:
                tier = 2
            elif nf >= args.min_full_copies or (args.min_half_copies > 0
                                                 and nh >= args.min_half_copies):
                tier = 3
            else:
                tier = 4
            new_label = label
            if args.relabel_unknown and dom_label and UNKNOWN.match(label.split('/')[0]):
                new_label = dom_label
            elif (args.relabel_conflicts and dom_label and class_group(dom_label)
                  and class_group(label) and class_group(dom_label) != class_group(label)):
                new_label = dom_label
            if tier == 0:
                action, fh_out = ('to_nonte', nonte) if args.mode != 'off' else ('keep', keep)
            elif tier in drop_tiers:
                action, fh_out = 'remove', removed
            else:
                action = 'keep_relabel' if new_label != label else 'keep'
                fh_out = keep
            counts[action] += 1
            counts[f'tier{tier}'] += 1
            fh_out.write(f'>{base}#{new_label}\n{wrap(seq)}\n')
            ev.write(f'{base}\t{tool(base)}\t{label}\t{new_label}\t{len(seq)}\t'
                     f'{"/".join(d) if d else ""}\t{",".join(sup)}\t{nf}\t{tier}\t{action}\t{nh}\t'
                     f'{",".join(pf)}\n')

    sys.stderr.write(
        'te_evidence: %d consensi; tiers %s; mode %s -> kept %d (relabelled %d), '
        'non-TE %d, removed %d\n' % (
            len(lib), ' '.join(f'{t}:{counts[f"tier{t}"]}' for t in ('T', 0, 1, 2, 3, 4)), args.mode,
            counts['keep'] + counts['keep_relabel'], counts['keep_relabel'],
            counts['to_nonte'], counts['remove']))


if __name__ == '__main__':
    main()
