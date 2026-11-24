#!/usr/bin/env python3
"""patch_module.py <pipeline checkout> - benchmark-only patch of the checkout's TETRIMMER module:
--tetrimmer_precomputed <dir> (a TEtrimmer output dir merged by tt_chunks.py) is copied into
the task and finished with TEtrimmer --continue_analysis. Without the param the task runs as
before (only its script text, and so its cache hash, changes). Keeps main.nf.pre_ttchunks."""
import os, shutil, sys

H = sys.argv[1]
p = os.path.join(H, 'modules/local/tetrimmer/main.nf')
src = open(p).read()
if 'tetrimmer_precomputed' in src:
    sys.exit(f'{p}: already patched')
old = """    \\${TETRIMMER_CMD} \\\\
        --input_file  ${library} \\\\
        --genome_file ${genome} \\\\
        --output_dir  ${meta.id}_tetrimmer \\\\
        --pfam_dir    ${pfam_dir} \\\\
        --num_threads ${task.cpus}
"""
new = """    # BENCHMARK CHECKOUT ONLY (not in the pipeline repository): --tetrimmer_precomputed <dir>
    # is a TEtrimmer output directory whose per-sequence work was done in chunks outside
    # Nextflow (tetrimmer_chunks/tt_chunks.py). TEtrimmer --continue_analysis on a copy of it
    # finds every sequence done and writes only its library-wide final outputs.
    TT_CA=
    if [[ -n "${params.tetrimmer_precomputed ?: ''}" ]]; then
        cp -r "${params.tetrimmer_precomputed ?: ''}" ${meta.id}_tetrimmer
        TT_CA=--continue_analysis
    fi
    \\${TETRIMMER_CMD} \\\\
        --input_file  ${library} \\\\
        --genome_file ${genome} \\\\
        --output_dir  ${meta.id}_tetrimmer \\\\
        --pfam_dir    ${pfam_dir} \\\\
        --num_threads ${task.cpus} \\${TT_CA}
"""
assert src.count(old) == 1, 'TEtrimmer call not found exactly once'
shutil.copy2(p, p + '.pre_ttchunks')
open(p, 'w').write(src.replace(old, new))
with open(os.path.join(H, 'BENCHMARK_PATCH.txt'), 'a') as f:
    f.write('modules/local/tetrimmer/main.nf: --tetrimmer_precomputed (tetrimmer_chunks/patch_module.py, '
            '2026-11-04): finish a chunked TEtrimmer run with --continue_analysis. Benchmark only; '
            'original in main.nf.pre_ttchunks.\n')
print(f'patched {p}')
