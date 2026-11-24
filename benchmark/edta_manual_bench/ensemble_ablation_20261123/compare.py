"""Paired ensemble ablation using the published vote and its exact input library."""
import argparse
import csv
import hashlib
import importlib.util
from pathlib import Path
import shutil
import sys


B = Path('/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark')
SP = 'caenorhabditis_briggsae'
ROOT = Path(__file__).resolve().parent
ORIGINAL = B / 'edta_manual_bench' / SP
VOTE = B / 'work/heldout_v10' / SP / 'work/e8/ff48941de188fabb886541df9de624'
GENOME = B / 'genomes/animal' / SP / f'{SP}.fna'
ARMS = ('ensemble_off', 'ensemble_on')


def fasta(path):
    records = {}
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name, _, label = line[1:].split()[0].partition('#')
            if name in records:
                raise ValueError(f'Duplicate FASTA record: {name}')
            records[name] = (label, [])
        elif line.strip():
            records[name][1].append(line.strip())
    return {name: (label, ''.join(seq)) for name, (label, seq) in records.items()}


def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_tsv(path, header, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def prepare():
    off = VOTE / f'{SP}.domain_relabelled.fa'
    on = VOTE / f'{SP}.classified_library.fa'
    voting = VOTE / f'{SP}.voting_table.tsv'
    if (VOTE / '.exitcode').read_text().strip() != '0':
        raise RuntimeError('Saved ensemble vote did not finish successfully')
    if digest(on) != digest(ORIGINAL / 'teannotate_ref10/library.fa'):
        raise RuntimeError('Saved vote does not match the published masking library')
    a, b = fasta(off), fasta(on)
    if list(a) != list(b) or any(a[k][1] != b[k][1] for k in a):
        raise RuntimeError('Paired libraries must have identical sequence IDs, order and bases')
    with voting.open(newline='') as handle:
        calls = {r['repeat_id']: r['Ensemble'] for r in csv.DictReader(handle, delimiter='\t')}
    changes = [(k, a[k][0], b[k][0]) for k in a if a[k][0] != b[k][0]]
    for name, old, new in changes:
        if calls.get(name) != new or new.split('/')[0] not in ('DNA', 'RC'):
            raise RuntimeError(f'Unexpected label change: {name}: {old} -> {new}')
    if not changes:
        raise RuntimeError('No label changes found in the cached result')
    write_tsv(ROOT / 'relabelled_consensi.tsv', ('repeat_id', 'before', 'after'), changes)
    write_tsv(ROOT / 'library_counts.tsv', ('consensi', 'sent_to_vote', 'relabelled'), [(len(a), len(calls), len(changes))])
    copies = {
        off: ROOT / SP / ARMS[0] / 'library.fa',
        on: ROOT / SP / ARMS[1] / 'library.fa',
        ORIGINAL / 'teannotate_ref10' / f'{SP}.fna.out': ROOT / SP / ARMS[1] / f'{SP}.fna.out',
        ORIGINAL / 'reference' / f'{SP}.fna.out': ROOT / SP / 'reference' / f'{SP}.fna.out',
    }
    for src, dst in copies.items():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        if digest(src) != digest(dst):
            raise RuntimeError(f'Snapshot checksum mismatch: {dst}')
    (ROOT / 'results').mkdir(exist_ok=True)
    write_tsv(ROOT / 'inputs.tsv', ('source', 'snapshot', 'sha256'),
              [(str(src.resolve()), str(dst.relative_to(ROOT)), digest(dst)) for src, dst in copies.items()]
              + [(str(GENOME.resolve()), 'shared genome', digest(GENOME))])
    print(f'{len(a)} identical consensi, {len(calls)} voted, {len(changes)} relabelled', flush=True)


def scorer():
    path = B / 'edta_manual_bench/libtest_harmonised.py'
    spec = importlib.util.spec_from_file_location('original_harmonised', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.E = str(ROOT)
    return module


def coverage(out, size, library, changed, canon):
    """Union coverage in half-open genomic intervals, using RepeatMasker labels."""
    import numpy as np
    def blank():
        return {chrom: np.zeros(n, dtype=bool) for chrom, n in size.items()}
    all_bp, classified, unknown, affected = blank(), blank(), blank(), blank()
    observed = set()
    with out.open() as handle:
        for line in handle:
            f = line.split()
            if not f or not f[0].isdigit():
                continue
            chrom, start, end, name, label = f[4], int(f[5])-1, int(f[6]), f[9], f[10]
            if name not in library or not 0 <= start < end <= size[chrom]:
                raise ValueError(f'Unexpected RepeatMasker record: {line.strip()}')
            observed.add(name)
            all_bp[chrom][start:end] = True
            if canon(label) == 'Unknown':
                unknown[chrom][start:end] = True
            elif canon(label) != 'Unspecified':
                classified[chrom][start:end] = True
            if name in changed:
                affected[chrom][start:end] = True
    if not observed:
        raise RuntimeError(f'No RepeatMasker hits in {out}')
    return all_bp, classified, unknown, affected


def score():
    module = scorer()
    _, size = module.sizes(SP)
    with (ROOT / 'relabelled_consensi.tsv').open(newline='') as handle:
        changed = {r['repeat_id'] for r in csv.DictReader(handle, delimiter='\t')}
    counts, totals = {}, {}
    for arm in ARMS:
        path = ROOT / SP / arm
        library = fasta(path / 'library.fa')
        all_bp, classified, unknown, affected = coverage(path / f'{SP}.fna.out', size, library, changed, module.canon)
        totals[arm] = all_bp
        counts[arm] = [sum(int(mask[c].sum()) for c in size) for mask in (all_bp, classified, affected)]
        counts[arm].append(sum(int((unknown[c] & ~classified[c]).sum()) for c in size))
        if arm == ARMS[1]:
            on_classified, on_affected = classified, affected
        else:
            off_classified, off_affected = classified, affected
    mask_header = ('metric', *ARMS, 'delta_on_minus_off')
    write_tsv(ROOT / 'masked_bases.tsv', mask_header,
              [(k, a, b, b-a) for k, a, b in zip(('any_masked_bp', 'classified_te_bp', 'relabelled_consensus_bp', 'unknown_without_classified_bp'),
                                                  counts[ARMS[0]], counts[ARMS[1]])])
    off, on = (totals[arm] for arm in ARMS)
    write_tsv(ROOT / 'coverage_overlap.tsv', ('metric', 'bp'), [
        ('both_masked', sum(int((off[c] & on[c]).sum()) for c in size)),
        ('off_only_masked', sum(int((off[c] & ~on[c]).sum()) for c in size)),
        ('on_only_masked', sum(int((on[c] & ~off[c]).sum()) for c in size)),
        ('on_only_classified', sum(int((on_classified[c] & ~off_classified[c]).sum()) for c in size)),
        ('off_only_classified', sum(int((off_classified[c] & ~on_classified[c]).sum()) for c in size)),
        ('relabelled_consensus_bp_in_both_arms', sum(int((off_affected[c] & on_affected[c]).sum()) for c in size)),
    ])
    sys.argv = [module.__file__, SP, *ARMS]
    module.main()
    scored = ROOT / 'results' / 'Caenorhabditis_briggsae_bench_harmonised.csv'
    with scored.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    types = ('ltr', 'nonltr', 'line', 'sine', 'tir', 'mite', 'helitron', 'any', 'any_dontcare')
    paired = []
    for key in types:
        selected = []
        for arm in ARMS:
            matches = [r for r in rows if r['Method'] == arm and r['TE type'] == key]
            if len(matches) != 1:
                raise RuntimeError(f'Missing or duplicate score for {arm}/{key}')
            selected.append(matches[0])
        a, b = selected
        ref = int(a['TP']) + int(a['FN'])
        if ref != int(b['TP']) + int(b['FN']):
            raise RuntimeError(f'Reference coverage differs between arms for {key}')
        paired.append((key, ref, a['F1'], b['F1'], float(b['F1']) - float(a['F1'])))
    write_tsv(ROOT / 'paired_f1.tsv', ('category', 'reference_bp', *ARMS, 'delta_on_minus_off'), paired)
    (ROOT / 'COMPLETE').touch()


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('stage', choices=('prepare', 'score'))
    {'prepare': prepare, 'score': score}[cli.parse_args().stage]()
