#!/usr/bin/env python3
"""Pfam TE-domain evidence helpers for TE_EVIDENCE.

  te_pfam.py subset <Pfam-A.hmm> <out.hmm>     extract the TE-protein profiles listed below
  te_pfam.py orfs   <library.fa> <out.faa>     six-frame ORF fragments (>= --min-aa aa)

REXdb (TEsorter) profiles are clade-level models built from well-studied plant and
animal TEs; Pfam's TE families are built from very diverse sequences and detect more
distant relatives, which matters for unexplored lineages. Only families specific to
transposable elements are used: reverse transcriptase, RNase H of retroelements,
integrase, retroelement proteases and Gag, DDE transposases and their DNA-binding
domains, Helitron Rep/Helicase, tyrosine recombinase (DIRS/Crypton/Ngaro) and the
protein-primed DNA polymerase of Mavericks. Families shared with cellular genes
(cellular RNase H, zinc knuckles, chromodomains, PIF1 helicase, GIY-YIG, THAP) are
left out, so a gene family mistaken for a TE does not pass on them.
ORFs are stop-to-stop fragments, not Met-started, so degraded copies still count.
"""
import sys

TE_PFAM = [
    # reverse transcriptase and retroelement RNase H
    'RVT_1', 'RVT_2', 'RVT_3', 'RVT_N', 'RVT_thumb', 'RT_RNaseH', 'RT_RNaseH_2', 'zf-RVT',
    # integrase
    'rve', 'rve_2', 'rve_3', 'Integrase_H2C2', 'zf-H2C2', 'Integrase_Zn',
    # retroelement protease and Gag
    'RVP', 'RVP_2', 'Asp_protease_2', 'gag-asp_proteas', 'Retrotrans_gag', 'Retrotran_gag_2',
    'Retrotran_gag_3', 'UBN2', 'UBN2_2', 'UBN2_3', 'gag_pre-integrs',
    # non-LTR (LINE) endonuclease / ORF1
    'Exo_endo_phos_2', 'Transposase_22',
    # tyrosine recombinase (DIRS, Ngaro, Crypton)
    'Phage_integrase',
    # DDE transposases and their DNA-binding domains
    'DDE_1', 'DDE_3', 'DDE_Tnp_1', 'DDE_Tnp_1_2', 'DDE_Tnp_1_3', 'DDE_Tnp_1_4', 'DDE_Tnp_1_5',
    'DDE_Tnp_1_6', 'DDE_Tnp_1_7', 'DDE_Tnp_2', 'DDE_Tnp_4', 'DDE_Tnp_IS1595', 'DDE_Tnp_ISL3',
    'Transposase_1', 'Transposase_11', 'Transposase_20', 'Transposase_21', 'Transposase_23',
    'Transposase_24', 'Transposase_28', 'Transposase_mut', 'Transpos_assoc', 'Dimer_Tnp_hAT',
    'DUF4413', 'Tnp_P_element', 'Tnp_P_element_C', 'MULE', 'Plant_tran', 'HTH_Tnp_Tc5',
    'HTH_Tnp_Tc3_2', 'HTH_Tnp_Tc3_1', 'HTH_Tnp_1', 'HTH_Tnp_4', 'Tc5_transposase',
    # Helitron Rep / helicase
    'Helitron_like_N', 'Herpes_Helicase',
    # Maverick / Polinton protein-primed DNA polymerase
    'DNA_pol_B_2',
]


def subset(pfam, out):
    want, found = set(TE_PFAM), set()
    keep, rec = False, []
    with open(pfam) as fh, open(out, 'w') as o:
        for line in fh:
            if line.startswith('HMMER3'):
                rec, keep = [line], False
                continue
            rec.append(line)
            if line.startswith('NAME  '):
                name = line.split()[1]
                keep = name in want
                if keep:
                    found.add(name)
            elif line.startswith('//'):
                if keep:
                    o.writelines(rec)
                rec, keep = [], False
    missing = sorted(want - found)
    sys.stderr.write('te_pfam: %d TE profiles extracted%s\n' % (
        len(found), '; not in this Pfam release: ' + ','.join(missing) if missing else ''))


CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def translate(s):
    return ''.join(CODON.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def orfs(lib, out, min_aa=80):
    n = 0
    name, buf = None, []

    def emit(o, name, seq):
        nonlocal n
        seq = seq.upper()
        rc = seq.translate(COMP)[::-1]
        for strand, s in (('+', seq), ('-', rc)):
            for frame in range(3):
                for k, frag in enumerate(translate(s[frame:]).split('*')):
                    if len(frag) >= min_aa:
                        o.write(f'>{name}|{strand}{frame}|{k}\n{frag}\n')
                        n += 1

    with open(lib) as fh, open(out, 'w') as o:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    emit(o, name, ''.join(buf))
                name, buf = line[1:].split()[0].split('#')[0], []
            else:
                buf.append(line.strip())
        if name is not None:
            emit(o, name, ''.join(buf))
    sys.stderr.write(f'te_pfam: {n} ORF fragments >= {min_aa} aa\n')


if __name__ == '__main__':
    if len(sys.argv) < 4 or sys.argv[1] not in ('subset', 'orfs'):
        sys.exit(__doc__)
    if sys.argv[1] == 'subset':
        subset(sys.argv[2], sys.argv[3])
    else:
        orfs(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 80)
