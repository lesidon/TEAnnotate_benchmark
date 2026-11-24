#!/usr/bin/env python3
"""Rewrite FASTA headers for CREATE input.

CREATE's create_test.py reads each sequence's TRUE label from the text BEFORE '#'
(a benchmark convention) and crashes (find_node_path IndexError) when that isn't a
taxonomy node — which is exactly our '#Unknown'. Rewrite every header to
">TE#<orig_id>": the valid root node "TE" goes before '#', and the real family id
after it, where it survives into CREATE's relabeled_seqs.fasta so
normalize_create.py can recover it.

    create_prep_headers.py IN.fa OUT.fa
"""
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit('usage: create_prep_headers.py IN.fa OUT.fa')
    inp, out = sys.argv[1], sys.argv[2]
    with open(inp) as fh, open(out, 'w') as o:
        for line in fh:
            if line.startswith('>'):
                base = line[1:].split()[0].split('#')[0]
                o.write('>TE#%s\n' % base)
            else:
                o.write(line)


if __name__ == '__main__':
    main()
