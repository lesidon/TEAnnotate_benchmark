#!/usr/bin/env python3
"""Write an MMseqs2 clustering of a TE library in CD-HIT-EST's output format.

CD_HIT_EST clusters large libraries with MMseqs2 instead of CD-HIT (see that module).
Everything downstream reads CD-HIT's files, so this script writes them:

  <out>          representatives, in input order (as cd-hit-est writes them)
  <out>.clstr    one block per cluster, as cd-hit-est -p 1 -d 0 writes it:
                   0   2535nt, >rep#LTR/Gypsy... *
                   1   1088nt, >member#LTR/Gypsy... at 1:1088:120:1207/+/91.20%
                 member start:end on the member, then on the representative, the
                 strand and the identity. DEDUP_REFINE reads the representative span.

Inputs: the clustered FASTA, MMseqs2's cluster table (createtsv: rep<TAB>member) and an
alignment table of members against representatives (convertalis with
query,target,fident,qstart,qend,tstart,tend; query = member, target = representative).
A member without an alignment row (rare: the search did not report the pair) is
written without coordinates, which DEDUP_REFINE treats as "span unknown".
"""
import argparse
import sys


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='the library that was clustered')
    ap.add_argument('--clusters', required=True, help='mmseqs createtsv: rep<TAB>member')
    ap.add_argument('--alignments', required=True,
                    help='mmseqs convertalis: member, rep, fident, qstart, qend, tstart, tend')
    ap.add_argument('--out', required=True, help='representative FASTA; <out>.clstr is written too')
    args = ap.parse_args()

    seqs, order = read_fasta(args.fasta)
    rank = {n: i for i, n in enumerate(order)}

    members = {}                                  # rep -> [member, ...]
    rep_of = {}
    with open(args.clusters) as fh:
        for line in fh:
            rep, mem = line.rstrip('\n').split('\t')[:2]
            if rep not in seqs or mem not in seqs:
                sys.exit(f'mmseqs_to_clstr: {rep if rep not in seqs else mem} is not in {args.fasta}')
            if mem in rep_of:
                sys.exit(f'mmseqs_to_clstr: {mem} is in two clusters')
            rep_of[mem] = rep
            members.setdefault(rep, [])
            if mem != rep:
                members[rep].append(mem)
    missing = [n for n in order if n not in rep_of]
    if missing:
        sys.exit(f'mmseqs_to_clstr: {len(missing)} sequences are in no cluster, e.g. {missing[0]}')

    # best-scoring row (first reported) per member on its own representative
    aln = {}
    with open(args.alignments) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            mem, rep = f[0], f[1]
            if rep_of.get(mem) != rep or mem == rep or mem in aln:
                continue
            ident = float(f[2]) * 100
            qs, qe, ts, te = (int(x) for x in f[3:7])
            # MMseqs2 reports a minus-strand hit with one span reversed
            strand = '-' if (qs > qe) != (ts > te) else '+'
            aln[mem] = (min(qs, qe), max(qs, qe), min(ts, te), max(ts, te), strand, ident)

    reps = sorted(members, key=lambda r: rank[r])
    with open(args.out, 'w') as fa, open(args.out + '.clstr', 'w') as cl:
        for c, rep in enumerate(reps):
            seq = seqs[rep]
            fa.write(f'>{rep}\n')
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i + 60] + '\n')
            cl.write(f'>Cluster {c}\n')
            cl.write(f'0\t{len(seq)}nt, >{rep}... *\n')
            # cd-hit lists members by length, longest first
            for k, mem in enumerate(sorted(members[rep], key=lambda m: (-len(seqs[m]), rank[m])), 1):
                a = aln.get(mem)
                where = (f'at {a[0]}:{a[1]}:{a[2]}:{a[3]}/{a[4]}/{a[5]:.2f}%' if a else 'at 0.00%')
                cl.write(f'{k}\t{len(seqs[mem])}nt, >{mem}... {where}\n')
    n_mem = sum(len(v) for v in members.values())
    print(f'mmseqs_to_clstr: {len(order)} sequences -> {len(reps)} clusters; '
          f'{n_mem - sum(1 for m in rep_of if m != rep_of[m] and m in aln)} of {n_mem} members without an alignment row',
          file=sys.stderr)


if __name__ == '__main__':
    main()
