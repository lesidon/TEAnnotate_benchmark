#!/usr/bin/env python3
"""make_combos.py <ablation_dir> <species>... - raw libraries for every combination of the
discovery tools run for a species (libs/<sp>/<tool>.fa), as the pipeline's CAT_TE_LIBRARIES
would build them: tools in name order, a repeated ID renamed <id>_<tool> (then
<id>_<tool>2, ...), and a provenance table (name<TAB>tool) for the evidence gate.
HELIANO alone is skipped: it finds Helitrons only, so it is never a route on its own.
Writes combos/<sp>/<combo>/{raw.fa,provenance.tsv}; <combo> = tools joined by '+'.
"""
import itertools
import os
import sys

A = sys.argv[1]
TOOLS = ['earlgrey', 'edta', 'heliano', 'hite']   # CAT_TE_LIBRARIES order (sorted names)

for sp in sys.argv[2:]:
    # every tool must have run (a library file, possibly empty); a tool that found
    # nothing is left out of the combinations, which would only repeat those without it
    missing = [t for t in TOOLS if not os.path.exists(f'{A}/libs/{sp}/{t}.fa')]
    if missing:
        sys.exit(f'{sp}: tool libraries missing: {",".join(missing)}')
    have = [t for t in TOOLS if os.path.getsize(f'{A}/libs/{sp}/{t}.fa') > 0]
    empty = [t for t in TOOLS if t not in have]
    if empty:
        print(f'{sp}: {",".join(empty)} found nothing; left out of the combinations')
    n = 0
    for k in range(1, len(have) + 1):
        for combo in itertools.combinations(have, k):
            if combo == ('heliano',):
                continue
            name = '+'.join(combo)
            d = f'{A}/combos/{sp}/{name}'
            os.makedirs(d, exist_ok=True)
            seen = set()
            with open(f'{d}/raw.fa', 'w') as fa, open(f'{d}/provenance.tsv', 'w') as pv:
                for t in combo:
                    for line in open(f'{A}/libs/{sp}/{t}.fa'):
                        if line.startswith('>'):
                            tok, rest = line[1:].rstrip('\n'), ''
                            if ' ' in tok:
                                tok, rest = tok.split(' ', 1)
                                rest = ' ' + rest
                            base, sep, lab = tok.partition('#')
                            new, i = base, 1
                            while new in seen:
                                new = f'{base}_{t}' + (str(i) if i > 1 else '')
                                i += 1
                            seen.add(new)
                            pv.write(f'{new}\t{t}\n')
                            fa.write(f'>{new}{sep}{lab}{rest}\n')
                        else:
                            fa.write(line)
            n += 1
    print(f'{sp}: {n} combinations from {",".join(have)}')
