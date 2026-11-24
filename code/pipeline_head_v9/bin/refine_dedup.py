#!/usr/bin/env python3
"""Class-aware refinement of a CD-HIT-EST deduplicated TE library.

CD-HIT keeps the longest sequence of each cluster as its representative, with that
sequence's own label. When a cluster mixes TE classes, the minority class silently
disappears. De novo tools also emit chimeric consensi, e.g. an LTR-int that spans a
LINE inserted next to it, and those LINE segments are then masked as LTR. This script
repairs both, in two steps.

Cluster step (reads the cd-hit-est .clstr):
  relabel   an Unknown representative takes the class of its classified members
            (length-weighted majority; the full label of that class's longest member).
  rescue    a classified member whose class differs from the representative's is kept
            as its own family (the longest member per foreign class).

Alignment step (all-against-all blastn of the classified library; skip with --no-split):
  For every pair of classified sequences from different classes that align over at
  least --min-segment bp at >= --min-identity %, the sequence in which the shared
  segment is the SMALLER fraction of its length is the host (the one carrying a foreign
  insertion); the other is the donor family.
  split     a host whose foreign segments cover at most --max-foreign-frac of it is a
            chimera: the segments are cut out and the pieces of at least
            --min-fragment bp are kept as <name>_partN with the host's label.
  conflict  when the shared segment is the larger part of BOTH sequences, the two carry
            different labels for essentially the same element. Which label is right is
            not decidable here, so the pair is only reported; with --unlabel-conflicts
            both become Unknown and the re-classification ensemble decides.

Classes are compared at the level of major groups (LTR, LINE, SINE, PLE, DNA, RC, ...).
EDTA's DNAauto/DNAnona, MITE and TIR count as DNA and Helitron counts as RC, so
tool-specific spellings of the same class are not conflicts.

Writes the refined library and a TSV report of every action.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

UNKNOWN_EXACT = {'', 'na', 'none'}
UNKNOWN_PREFIX = ('unknown', 'unspecified', 'unclassified')
GROUP = {
    'ltr': 'LTR', 'dirs': 'LTR', 'erv': 'LTR',
    'line': 'LINE', 'sine': 'SINE', 'ple': 'PLE', 'penelope': 'PLE',
    'dna': 'DNA', 'dnaauto': 'DNA', 'dnanona': 'DNA', 'mite': 'DNA', 'tir': 'DNA',
    'rc': 'RC', 'helitron': 'RC',
}
MEMBER_RE = re.compile(
    r'(\d+)(?:nt|aa), >(.+?)\.\.\. '
    r'(?:(\*)|at (?:(\d+):(\d+):(\d+):(\d+)/)?(?:([+-])/)?([\d.]+)%)')


def label_of(name):
    return name.split('#', 1)[1] if '#' in name else ''


def base_of(name):
    return name.split('#', 1)[0]


def is_unknown(label):
    lab = label.strip().lower()
    top = lab.split('/', 1)[0]
    return (lab in UNKNOWN_EXACT or top in UNKNOWN_EXACT
            or lab.startswith(UNKNOWN_PREFIX) or top.startswith(UNKNOWN_PREFIX))


def group_of(label):
    if is_unknown(label):
        return None
    # Helitrons carry different top-level tokens per tool (EDTA "DNA/Helitron",
    # HELIANO and RepeatMasker "RC/Helitron"): one class, whatever the prefix.
    if 'helitron' in label.lower():
        return 'RC'
    top = label.split('/', 1)[0].strip()
    return GROUP.get(top.lower(), top)


def read_fasta(path):
    seqs, order, name, buf = {}, [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                order.append(name)
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs, order


def read_clusters(path):
    clusters, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith('>Cluster'):
                cur = []
                clusters.append(cur)
                continue
            m = MEMBER_RE.search(line)
            if not m:
                continue
            ln, name, star, qs, qe, rs, re_, strand, ident = m.groups()
            cur.append({
                'name': name.split()[0], 'len': int(ln), 'rep': star == '*',
                # member coordinates on the representative (1-based, inclusive)
                'rep_span': (min(int(rs), int(re_)), max(int(rs), int(re_))) if rs else None,
                'ident': ident,
            })
    return clusters


def wrap(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def merge_intervals(length, cuts):
    """Merged 1-based inclusive intervals covered by `cuts` (clipped to the length)."""
    merged = []
    for s, e in sorted((max(1, s), min(length, e)) for s, e in cuts):
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def covered_bp(length, spans):
    return sum(e - s + 1 for s, e in merge_intervals(length, spans))


def subtract(length, cuts):
    """1-based inclusive pieces of [1, length] left after removing `cuts`."""
    pieces, pos = [], 1
    for s, e in merge_intervals(length, cuts):
        if s > pos:
            pieces.append((pos, s - 1))
        pos = max(pos, e + 1)
    if pos <= length:
        pieces.append((pos, length))
    return pieces


def blast_pairs(query, subject, threads, min_identity, task='dc-megablast', evalue=1e-10):
    """blastn `query` against `subject` (both {id: seq}); yields
    (qid, sid, pident, length, (qs, qe), (ss, se)) with 1-based inclusive spans."""
    with tempfile.TemporaryDirectory() as tmp:
        qf, sf = os.path.join(tmp, 'q.fa'), os.path.join(tmp, 's.fa')
        for path, recs in ((qf, query), (sf, subject)):
            with open(path, 'w') as fh:
                for k, v in recs.items():
                    fh.write(f'>{k}\n{wrap(v)}\n')
        subprocess.run(['makeblastdb', '-in', sf, '-dbtype', 'nucl', '-out', sf],
                       check=True, stdout=subprocess.DEVNULL)
        res = subprocess.run(
            ['blastn', '-task', task, '-query', qf, '-db', sf, '-evalue', str(evalue),
             '-perc_identity', str(min_identity), '-num_threads', str(threads),
             '-max_target_seqs', '1000', '-dust', 'no',
             '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send'],
            check=True, capture_output=True, text=True)
    for line in res.stdout.splitlines():
        q, s, pid, ln, qs, qe, ss, se = line.split('\t')
        qs, qe, ss, se = int(qs), int(qe), int(ss), int(se)
        yield q, s, float(pid), int(ln), (min(qs, qe), max(qs, qe)), (min(ss, se), max(ss, se))


def count_full_copies(seqs, genome, threads):
    """{name: genomic hits covering >= 90% of the sequence at >= 75% identity}.
    The BLAST database goes in the working directory: a genome can be large."""
    with tempfile.TemporaryDirectory(dir='.') as tmp:
        qf, db = os.path.join(tmp, 'q.fa'), os.path.join(tmp, 'genome')
        ids = {}
        with open(qf, 'w') as fh:
            for i, (n, s) in enumerate(seqs.items()):
                ids[f'm{i}'] = n
                fh.write(f'>m{i}\n{wrap(s)}\n')
        subprocess.run(['makeblastdb', '-in', genome, '-dbtype', 'nucl', '-out', db],
                       check=True, stdout=subprocess.DEVNULL)
        res = subprocess.run(
            ['blastn', '-task', 'megablast', '-query', qf, '-db', db, '-evalue', '1e-20',
             '-perc_identity', '75', '-qcov_hsp_perc', '90', '-max_target_seqs', '100',
             '-max_hsps', '10', '-num_threads', str(threads), '-outfmt', '6 qseqid'],
            check=True, capture_output=True, text=True)
    copies = defaultdict(int)
    for line in res.stdout.splitlines():
        if line.strip():
            copies[ids[line.strip()]] += 1
    return copies


def representative_step(combined, dedup, order, clusters, genome, threads, report, counts):
    """Make the member with the most full-length genomic copies each cluster's
    representative (within the representative's class). CD-HIT keeps the longest
    member, which can be an over-extended or chimeric consensus that absorbed an
    exact one: a correct consensus matches real insertions end to end, extra sequence
    does not. Ties keep CD-HIT's choice. Edits `dedup`, `order` and the clusters."""
    names = {m['name'] for c in clusters if len(c) > 1 for m in c if m['name'] in combined}
    if not names:
        return
    copies = count_full_copies({n: combined[n] for n in names}, genome, threads)
    pos = {n: i for i, n in enumerate(order)}
    for ci, members in enumerate(clusters):
        reps = [m for m in members if m['rep']]
        if len(members) < 2 or not reps or reps[0]['name'] not in pos:
            continue
        rep = reps[0]
        g = group_of(label_of(rep['name']))
        cands = [m for m in members
                 if m['name'] in combined and group_of(label_of(m['name'])) == g]
        best = max(cands, key=lambda m: (copies.get(m['name'], 0), m['rep'], m['len']))
        if best is rep or copies.get(best['name'], 0) <= copies.get(rep['name'], 0):
            continue
        i = pos.pop(rep['name'])
        order[i], pos[best['name']] = best['name'], i
        dedup[best['name']] = combined[best['name']]
        dedup.pop(rep['name'], None)
        rep['rep'], best['rep'] = False, True
        rep['rep_span'] = best['rep_span'] = None
        report.append((f'cluster{ci}', 'representative', rep['name'], best['name'],
                       f"{copies.get(best['name'], 0)} vs {copies.get(rep['name'], 0)} "
                       f"full-length copies ({best['len']} vs {rep['len']} bp)"))
        counts['representative'] += 1


def cluster_step(combined, dedup, order, clusters, report, counts):
    """relabel / rescue; returns the library as an ordered list of (name, seq)."""
    out_reps, rescued = {}, []
    for ci, members in enumerate(clusters):
        reps = [m for m in members if m['rep']]
        if not reps:
            continue
        rname = reps[0]['name']
        rgroup = group_of(label_of(rname))
        classified = [m for m in members if not m['rep'] and group_of(label_of(m['name']))]

        if rgroup is None and classified:
            weight = defaultdict(int)
            for m in classified:
                weight[group_of(label_of(m['name']))] += m['len']
            rgroup = max(weight, key=weight.get)
            donor = max((m for m in classified if group_of(label_of(m['name'])) == rgroup),
                        key=lambda m: m['len'])
            new = f"{base_of(rname)}#{label_of(donor['name'])}"
            out_reps[rname] = new
            report.append((f'cluster{ci}', 'relabel', rname, new,
                           f"from {donor['name']} ({weight[rgroup]} bp of {rgroup} members)"))
            counts['relabel'] += 1

        foreign = defaultdict(list)
        for m in classified:
            g = group_of(label_of(m['name']))
            if rgroup is not None and g != rgroup:
                foreign[g].append(m)
        for g, ms in foreign.items():
            keep = max(ms, key=lambda m: m['len'])
            if keep['name'] not in combined:
                sys.stderr.write(f"refine_dedup: {keep['name']} not in combined library; skipped\n")
                continue
            rescued.append((keep['name'], combined[keep['name']]))
            report.append((f'cluster{ci}', 'rescue', keep['name'], keep['name'],
                           f"{g} member of {rgroup} rep {rname} "
                           f"({keep['ident']}%, span {keep['rep_span']}); "
                           f"{len(ms) - 1} more {g} member(s) not kept"))
            counts['rescue'] += 1

    lib = [(out_reps.get(n, n), dedup[n]) for n in order]
    seen = {n for n, _ in lib}
    lib += [(n, s) for n, s in rescued if n not in seen]
    return lib


def chimera_step(lib, args, report, counts):
    """split chimeric hosts / unlabel cross-class conflicts; returns the new library."""
    ids = {f's{i}': (n, s) for i, (n, s) in enumerate(lib)}
    cls = {k: (n, s) for k, (n, s) in ids.items() if group_of(label_of(n))}
    if len(cls) < 2:
        return lib
    # per ordered pair: spans on the query and on the subject
    qspans, sspans = defaultdict(list), defaultdict(list)
    seqs = {k: s for k, (n, s) in cls.items()}
    for q, s, pid, ln, qsp, ssp in blast_pairs(seqs, seqs, args.threads, args.min_identity):
        if q == s or ln < args.min_segment:
            continue
        if group_of(label_of(cls[q][0])) == group_of(label_of(cls[s][0])):
            continue
        key = tuple(sorted((q, s)))
        # store spans keyed by the unordered pair, per member
        qspans[(key, q)].append(qsp)
        sspans[(key, s)].append(ssp)
    spans = defaultdict(list)
    for d in (qspans, sspans):
        for k, v in d.items():
            spans[k].extend(v)
    pairs = {k for k, _ in spans}

    cuts, conflict = defaultdict(list), set()
    for a, b in pairs:
        la, lb = len(cls[a][1]), len(cls[b][1])
        fa = covered_bp(la, spans[((a, b), a)]) / la
        fb = covered_bp(lb, spans[((a, b), b)]) / lb
        host, donor = (a, b) if fa <= fb else (b, a)
        fh = min(fa, fb)
        if fh <= args.max_foreign_frac:
            cuts[host].extend(spans[((a, b), host)])
            report.append(('align', 'segment', cls[host][0], cls[donor][0],
                           f'{fh:.0%} of host matches the {group_of(label_of(cls[donor][0]))} family'))
        else:
            if args.unlabel_conflicts:
                conflict.update((a, b))
            report.append(('align', 'conflict_pair', cls[a][0], cls[b][0],
                           f'shared segment covers {fa:.0%} / {fb:.0%}; labels disagree'))
            counts['conflict_pair'] += 1

    new = []
    for k, (n, s) in ids.items():
        if k in conflict:
            un = f'{base_of(n)}#Unknown'
            new.append((un, s))
            report.append(('align', 'conflict', n, un, 'label set to Unknown for the ensemble'))
            counts['conflict'] += 1
        elif k in cuts:
            pieces = [(a, b) for a, b in subtract(len(s), cuts[k]) if b - a + 1 >= args.min_fragment]
            parts = [(f'{base_of(n)}_part{i}#{label_of(n)}', s[a - 1:b])
                     for i, (a, b) in enumerate(pieces, 1)]
            new.extend(parts)
            cut = covered_bp(len(s), cuts[k])
            report.append(('align', 'split', n, ','.join(h for h, _ in parts) or '(dropped)',
                           f'chimera: {cut}/{len(s)} bp of other-class sequence cut; pieces {pieces}'))
            counts['split'] += 1
        else:
            new.append((n, s))
    return new


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--combined', required=True, help='all sequences given to CD-HIT')
    ap.add_argument('--dedup', required=True, help='CD-HIT representatives (FASTA)')
    ap.add_argument('--clstr', required=True, help='CD-HIT .clstr (run with -p 1)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--report', required=True)
    ap.add_argument('--no-split', action='store_true', help='skip the alignment step')
    ap.add_argument('--max-foreign-frac', type=float, default=0.5)
    ap.add_argument('--min-segment', type=int, default=500,
                    help='min aligned length (bp) of a cross-class segment')
    ap.add_argument('--min-identity', type=float, default=80.0)
    ap.add_argument('--min-fragment', type=int, default=80)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--unlabel-conflicts', action='store_true',
                    help='set both members of a cross-class conflict pair to Unknown')
    ap.add_argument('--genome', default=None,
                    help='choose each cluster\'s representative by full-length genomic copies')
    args = ap.parse_args()

    combined, _ = read_fasta(args.combined)
    dedup, order = read_fasta(args.dedup)
    clusters = read_clusters(args.clstr)
    report, counts = [], defaultdict(int)

    if args.genome:
        representative_step(combined, dedup, order, clusters, args.genome, args.threads,
                            report, counts)
    lib = cluster_step(combined, dedup, order, clusters, report, counts)
    if not args.no_split:
        lib = chimera_step(lib, args, report, counts)

    with open(args.out, 'w') as out:
        for header, seq in lib:
            out.write(f'>{header}\n{wrap(seq)}\n')
    with open(args.report, 'w') as fh:
        fh.write('step\taction\tsequence\tresult\tdetail\n')
        for row in report:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    sys.stderr.write(
        f"refine_dedup: {len(order)} representatives -> {len(lib)} sequences; "
        f"representative by copies {counts['representative']}, "
        f"relabel {counts['relabel']}, rescue {counts['rescue']}, "
        f"split {counts['split']}, conflict pairs {counts['conflict_pair']} "
        f"(unlabelled {counts['conflict']})\n")


if __name__ == '__main__':
    main()
