#!/usr/bin/env python3
"""RIP (repeat-induced point mutation) scan of a masked fungal genome.

RIP mutates C to T in duplicated DNA, mostly at CpA, so RIP-affected sequence is rich
in TpA and poor in CpA/TpG. Three complementary views are reported:

  1. Mutations per TE copy (RIPCAL-style; Hane & Oliver 2008). Every RepeatMasker
     alignment of a genomic copy to its consensus is walked base by base. A C->T at a
     consensus CpA (or G->A at TpG, the same mutation on the other strand) is a RIP-type
     transition; every other C->T / G->A is not. The enrichment is the RIP-type rate per
     CpA/TpG site over the rate per other C/G site (1 = no context preference). A
     consensus built from copies that were all mutated before they diverged already
     carries the TpA, so old, heavily RIP-mutated families can show few RIP-type
     changes against it; the indices below do not depend on the consensus.
  2. RIP indices per TE copy (Margolin et al. 1998; Lewis et al. 2009):
       product   = TpA / ApT
       substrate = (CpA + TpG) / (ApC + GpT)
       composite = product - substrate
     computed on the genomic sequence of the copy (copies shorter than --min-copy-len
     are not indexed). A copy is RIP-affected with the thresholds The RIPper uses for
     windows (van Wyk et al. 2019): product >= 1.15, substrate <= 0.75 and composite
     > 0. A positive composite index alone also passes many copies of AT-rich genomes
     without RIP (16% in Caenorhabditis elegans, against 4% with these thresholds).
  3. Large RIP-affected genomic regions, LRARs (as in The RIPper): windows of --window
     bp every --step bp are RIP-positive with the same thresholds; overlapping positive
     windows are merged and a region of at least --min-lrar bp is an LRAR.

Families are called RIP-affected when at least half of their indexed copies (and at
least --min-family-copies of them) are RIP-affected; the mutation enrichment is
reported beside it as independent evidence. Simple, low-complexity,
satellite and RNA repeats are left out.

Outputs (<prefix>.*): rip_copies.tsv, rip_families.tsv, rip_lrar.bed,
rip_cri.bedgraph (composite index per window) and rip_summary.tsv.
"""
import argparse
import gzip
import re
import statistics
import sys
from collections import defaultdict

NON_TE = {'Simple_repeat', 'Low_complexity', 'Satellite', 'rRNA', 'tRNA', 'snRNA',
          'scRNA', 'srpRNA', 'Structural_RNA', 'ARTEFACT'}
HEADER = re.compile(r'^\d+\s+\d+(\.\d+)?\s+\d+(\.\d+)?\s+\d+(\.\d+)?\s')
ALN_LINE = re.compile(r'^(?:C |  )(\S+)\s+(\d+)\s+([A-Za-z-]+)\s+(\d+)\s*$')
DINUC = ('TA', 'AT', 'CA', 'TG', 'AC', 'GT')


def read_fasta(path):
    opener = gzip.open if path.endswith('.gz') else open
    seqs, name, parts = {}, None, []
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(parts).upper()
    return seqs


def dinuc_counts(seq):
    # the six dinucleotides have two different letters, so str.count (non-overlapping)
    # equals the overlapping count
    return [seq.count(d) for d in DINUC]


def rip_indices(c):
    ta, at, ca, tg, ac, gt = c
    product = ta / at if at else None
    substrate = (ca + tg) / (ac + gt) if (ac + gt) else None
    composite = product - substrate if product is not None and substrate is not None else None
    return product, substrate, composite


def fmt(x, nd=4):
    return 'NA' if x is None else ('%.*f' % (nd, x))


def is_te(cls):
    return bool(cls) and cls.split('/')[0] not in NON_TE


def count_mutations(gen, con):
    """RIP-type and other C->T / G->A transitions of one alignment (consensus context)."""
    cols = [i for i, b in enumerate(con) if b != '-']
    cpa_sites = other_sites = rip = other = 0
    for k, i in enumerate(cols):
        c, g = con[i], gen[i]
        if c not in 'CG' or g not in 'ACGT':
            continue
        nxt = con[cols[k + 1]] if k + 1 < len(cols) else ''
        prv = con[cols[k - 1]] if k > 0 else ''
        rip_site = (c == 'C' and nxt == 'A') or (c == 'G' and prv == 'T')
        transition = (c == 'C' and g == 'T') or (c == 'G' and g == 'A')
        if rip_site:
            cpa_sites += 1
            rip += transition
        else:
            other_sites += 1
            other += transition
    return cpa_sites, rip, other_sites, other


def parse_align(path):
    """Yield (query, qstart, qend, family, class, rm_id, div, genome_aln, consensus_aln)."""
    hdr, gen, con, n = None, [], [], 0

    def flush():
        if hdr is not None and gen:
            return hdr + (''.join(gen).upper(), ''.join(con).upper())
        return None

    with open(path) as fh:
        for line in fh:
            if HEADER.match(line):
                rec = flush()
                if rec:
                    yield rec
                t = line.split()
                rep = next((x for x in t[8:] if '#' in x), t[9] if t[8] == 'C' else t[8])
                fam, _, cls = rep.partition('#')
                hdr = (t[4], int(t[5]), int(t[6]), fam, cls, t[-1], float(t[1]))
                gen, con, n = [], [], 0
                continue
            m = ALN_LINE.match(line.rstrip('\n'))
            if m and hdr is not None:
                (gen if n % 2 == 0 else con).append(m.group(3))
                n += 1
    rec = flush()
    if rec:
        yield rec


def ratio(a, b):
    return a / b if b else None


def rip_positive(p, sub, cri, a):
    return (p is not None and sub is not None and cri is not None and
            p >= a.min_product and sub <= a.max_substrate and cri > 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--genome', required=True)
    ap.add_argument('--align', required=True, help='RepeatMasker .align (-a)')
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--window', type=int, default=1000)
    ap.add_argument('--step', type=int, default=500)
    ap.add_argument('--min-product', type=float, default=1.15)
    ap.add_argument('--max-substrate', type=float, default=0.75)
    ap.add_argument('--min-lrar', type=int, default=4000)
    ap.add_argument('--min-copy-len', type=int, default=200)
    ap.add_argument('--min-family-copies', type=int, default=3)
    a = ap.parse_args()

    genome = read_fasta(a.genome)
    total_bp = sum(len(s) for s in genome.values())
    acgt_bp = sum(len(s) - s.count('N') for s in genome.values())
    gc_bp = sum(s.count('G') + s.count('C') for s in genome.values())

    # ── 3. Genome windows and LRARs ──────────────────────────────────────────
    lrars, rip_windows_bp, windows_eval = [], 0, 0
    with open(a.prefix + '.rip_cri.bedgraph', 'w') as bg:
        for ctg, seq in genome.items():
            runs, cur = [], None          # runs of overlapping RIP-positive windows
            for s in range(0, len(seq) - a.window + 1, a.step):
                w = seq[s:s + a.window]
                if len(w) - w.count('N') < a.window // 2:
                    continue
                windows_eval += 1
                p, sub, cri = rip_indices(dinuc_counts(w))
                if cri is not None:
                    bg.write('%s\t%d\t%d\t%.4f\n' % (ctg, s, s + a.window, cri))
                if rip_positive(p, sub, cri, a):
                    if cur and s <= cur[1]:
                        cur[1] = s + a.window
                        cur[2].append(cri)
                    else:
                        cur = [s, s + a.window, [cri]]
                        runs.append(cur)
            for st, en, cris in runs:
                rip_windows_bp += en - st
                if en - st >= a.min_lrar:
                    reg = seq[st:en]
                    acgt = len(reg) - reg.count('N')
                    gc = (reg.count('G') + reg.count('C')) / acgt if acgt else 0
                    lrars.append((ctg, st, en, statistics.mean(cris), gc))
    with open(a.prefix + '.rip_lrar.bed', 'w') as fh:
        for i, (ctg, st, en, cri, gc) in enumerate(lrars, 1):
            fh.write('%s\t%d\t%d\tLRAR_%d\t%.4f\t%.4f\n' % (ctg, st, en, i, cri, gc))

    # ── 1-2. Copies: mutations from the alignments, indices from the genome ──
    copies = {}
    missing = set()
    for q, qs, qe, fam, cls, rid, div, gen, con in parse_align(a.align):
        if not is_te(cls):
            continue
        key = (q, rid)
        cp = copies.get(key)
        if cp is None:
            cp = copies[key] = {'contig': q, 'start': qs, 'end': qe, 'family': fam,
                                'class': cls, 'frags': 0, 'bp': 0, 'div': [],
                                'dn': [0] * 6, 'mut': [0, 0, 0, 0]}
        cp['start'], cp['end'] = min(cp['start'], qs), max(cp['end'], qe)
        cp['frags'] += 1
        cp['bp'] += qe - qs + 1
        cp['div'].append(div)
        if q in genome:
            cp['dn'] = [x + y for x, y in zip(cp['dn'], dinuc_counts(genome[q][qs - 1:qe]))]
        else:
            missing.add(q)
        cp['mut'] = [x + y for x, y in zip(cp['mut'], count_mutations(gen, con))]
    if missing:
        sys.stderr.write('WARN rip_scan: %d aligned sequence(s) not in the genome FASTA '
                         '(no indices for their copies), e.g. %s\n' % (len(missing), sorted(missing)[0]))

    fams = defaultdict(lambda: {'class': '', 'copies': 0, 'bp': 0, 'cri': [], 'pos': 0,
                                'mut': [0, 0, 0, 0]})
    te_bp = rip_te_bp = indexed = pos_copies = 0
    tot_mut = [0, 0, 0, 0]
    with open(a.prefix + '.rip_copies.tsv', 'w') as fh:
        fh.write('contig\tstart\tend\trm_id\tfamily\tclass\tfragments\taligned_bp\tdiv\t'
                 'tpa_apt\tcpa_tpg_apc_gpt\tcomposite\trip_index\t'
                 'cpa_sites\tcpa_to_tpa\tother_cg_sites\tother_c_to_t\tenrichment\n')
        for (q, rid), cp in sorted(copies.items(), key=lambda kv: (kv[1]['contig'], kv[1]['start'])):
            p, sub, cri = rip_indices(cp['dn']) if cp['bp'] >= a.min_copy_len else (None, None, None)
            pos = rip_positive(p, sub, cri, a)
            call = 'NA' if cri is None else ('yes' if pos else 'no')
            cs, rp, osites, ot = cp['mut']
            enr = ratio(ratio(rp, cs) or 0, ratio(ot, osites)) if cs and osites and ot else None
            fh.write('%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%.2f\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n' % (
                q, cp['start'], cp['end'], rid, cp['family'], cp['class'], cp['frags'], cp['bp'],
                statistics.mean(cp['div']), fmt(p), fmt(sub), fmt(cri), call, cs, rp, osites, ot, fmt(enr, 2)))
            f = fams[cp['family']]
            f['class'] = cp['class']
            f['copies'] += 1
            f['bp'] += cp['bp']
            f['mut'] = [x + y for x, y in zip(f['mut'], cp['mut'])]
            tot_mut = [x + y for x, y in zip(tot_mut, cp['mut'])]
            te_bp += cp['bp']
            if cri is not None:
                f['cri'].append(cri)
                indexed += 1
                if pos:
                    f['pos'] += 1
                    pos_copies += 1
                    rip_te_bp += cp['bp']

    n_rip_fams = 0
    with open(a.prefix + '.rip_families.tsv', 'w') as fh:
        fh.write('family\tclass\tcopies\taligned_bp\tindexed_copies\trip_index_copies\t'
                 'rip_index_fraction\tmedian_composite\tcpa_to_tpa\tother_c_to_t\t'
                 'rip_transition_share\tenrichment\trip_affected\n')
        for fam, f in sorted(fams.items(), key=lambda kv: -kv[1]['bp']):
            n_idx = len(f['cri'])
            frac = ratio(f['pos'], n_idx)
            cs, rp, osites, ot = f['mut']
            share = ratio(rp, rp + ot)
            enr = ratio(ratio(rp, cs) or 0, ratio(ot, osites)) if cs and osites and ot else None
            if n_idx < a.min_family_copies:
                call = 'NA'
            else:
                call = 'yes' if frac >= 0.5 else 'no'
                n_rip_fams += call == 'yes'
            fh.write('%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n' % (
                fam, f['class'], f['copies'], f['bp'], n_idx, f['pos'], fmt(frac, 3),
                fmt(statistics.median(f['cri']) if f['cri'] else None), rp, ot,
                fmt(share, 3), fmt(enr, 2), call))

    cs, rp, osites, ot = tot_mut
    lrar_bp = sum(en - st for _, st, en, _, _ in lrars)
    lrar_gc_num = sum((en - st) * gc for _, st, en, _, gc in lrars)
    rows = [
        ('genome_bp', total_bp),
        ('genome_gc', fmt(gc_bp / acgt_bp if acgt_bp else None)),
        ('windows_evaluated', windows_eval),
        ('rip_window_bp', rip_windows_bp),
        ('rip_window_genome_pct', fmt(100 * rip_windows_bp / total_bp if total_bp else None, 2)),
        ('lrar_count', len(lrars)),
        ('lrar_bp', lrar_bp),
        ('lrar_genome_pct', fmt(100 * lrar_bp / total_bp if total_bp else None, 2)),
        ('lrar_gc', fmt(lrar_gc_num / lrar_bp if lrar_bp else None)),
        ('te_copies', len(copies)),
        ('te_copies_indexed', indexed),
        ('te_copies_rip_index', pos_copies),
        ('te_copies_rip_index_pct', fmt(100 * pos_copies / indexed if indexed else None, 2)),
        ('te_aligned_bp', te_bp),
        ('te_bp_rip_index_pct', fmt(100 * rip_te_bp / te_bp if te_bp else None, 2)),
        ('families', len(fams)),
        ('families_rip_affected', n_rip_fams),
        ('cpa_to_tpa', rp),
        ('other_c_to_t', ot),
        ('rip_transition_share', fmt(ratio(rp, rp + ot), 3)),
        ('rip_transition_enrichment', fmt(ratio(ratio(rp, cs) or 0, ratio(ot, osites)) if cs and osites and ot else None, 2)),
    ]
    with open(a.prefix + '.rip_summary.tsv', 'w') as fh:
        fh.write('metric\tvalue\n')
        for k, v in rows:
            fh.write('%s\t%s\n' % (k, v))
    sys.stderr.write('INFO rip_scan [%s]: %d LRARs (%s%% of genome); %s of %d indexed TE copies '
                     'RIP-positive; %d families RIP-affected; RIP-type transition enrichment %s\n' % (
                         a.prefix, len(lrars), rows[7][1], pos_copies, indexed, n_rip_fams, rows[-1][1]))


if __name__ == '__main__':
    main()
