"""Run the established maize scoring functions in isolated per-method roots."""
import argparse
import csv
import fcntl
import importlib.util
import os
from pathlib import Path
import shutil

B = Path('/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark').resolve()
ROOT = Path(__file__).resolve().parent
METHODS = ['earlgrey', 'hite', 'teannotate_wo_v13', 'edta2', 'heliano', 'teannotate_v13']
SUFFIXES = ['bench_edta_manual', 'bench_harmonised', 'family_level']


def original():
    spec = importlib.util.spec_from_file_location('maize_original', B / 'maize_mask/maize_scoring_20261117.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def scope(mod, root, method):
    selected = [pair for pair in mod.METHODS if pair[0] in ('reference', method)]
    if {pair[0] for pair in selected} != {'reference', method}:
        raise ValueError('Unknown method: ' + method)
    mod.METHODS = selected
    mod.TESTS = [method]
    mod.LABEL = {method: 'edta' if method == 'edta2' else method}
    mod.E = root / method / 'benchmark'
    mod.S = mod.E / 'zea_mays'
    mod.R = root / method / 'state'
    return mod


def redirect_harmoniser(source):
    old = 'E = f"{B}/edta_manual_bench"'
    if source.count(old) != 1:
        raise RuntimeError('Harmoniser destination declaration changed; refusing to patch')
    return source.replace(old, 'E = os.environ["MAIZE_SCORE_ROOT"]')


def prepare():
    source = B / 'edta_manual_bench'
    family = (source / 'family_level.sh').read_text()
    old = 'E=$B/edta_manual_bench'
    if family.count(old) != 1:
        raise RuntimeError('Family destination declaration changed')
    family = family.replace(old, 'E=${MAIZE_SCORE_ROOT:?isolated scoring root required}')
    harmoniser = redirect_harmoniser((source / 'libtest_harmonised.py').read_text())
    for method in METHODS:
        mod = scope(original(), ROOT, method)
        (mod.E / 'results').mkdir(parents=True, exist_ok=True)
        mod.R.mkdir(parents=True, exist_ok=True)
        (mod.E / 'family_level.sh').write_text(family)
        (mod.E / 'libtest_harmonised.py').write_text(harmoniser)
        shutil.copyfile(source / 'family_level.py', mod.E / 'family_level.py')
    print('Prepared six isolated scoring roots; only destination paths adapted', flush=True)


def validate_rows(rows, methods, field, values):
    for method in methods:
        for value in values:
            matches = [r for r in rows if r['Method'] == method and r[field] == value]
            if len(matches) != 1:
                raise RuntimeError('Missing or duplicate result: %s/%s' % (method, value))


def run(stage, method):
    mod = scope(original(), ROOT, method)
    # Check both reference and this method; another method never blocks this job.
    mod.validate_masks()
    os.environ['MAIZE_SCORE_ROOT'] = str(mod.E)
    with (mod.R / (stage + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        {'base': mod.base, 'family': mod.family}[stage]()


def merge():
    mod = original()
    if set(mod.TESTS) != set(METHODS):
        raise RuntimeError('Unexpected method inventory')
    collected = {}
    for suffix in SUFFIXES:
        rows, fields = [], None
        for method in METHODS:
            scoped = scope(original(), ROOT, method)
            stage = 'family' if suffix == 'family_level' else 'base'
            if not (scoped.R / (stage + '_complete.json')).is_file():
                raise RuntimeError('Incomplete stage: %s/%s' % (method, stage))
            with (scoped.E / 'results' / ('Zea_mays_' + suffix + '.csv')).open() as f:
                reader = csv.DictReader(f)
                if fields is not None and fields != reader.fieldnames:
                    raise RuntimeError('Inconsistent CSV columns')
                fields = reader.fieldnames
                part = list(reader)
            label = method if suffix == 'family_level' else scoped.LABEL[method]
            if {r['Method'] for r in part} != {label}:
                raise RuntimeError('Cross-method result contamination')
            rows.extend(part)
        labels = METHODS if suffix == 'family_level' else [mod.LABEL[m] for m in METHODS]
        key = 'class' if suffix == 'family_level' else 'TE type'
        needed = ['all'] if suffix == 'family_level' else ['any']
        if suffix == 'bench_harmonised':
            needed += ['any_dontcare']
        validate_rows(rows, labels, key, needed)
        collected[suffix] = (fields, rows)
    # Publish only after all six method results exist. Existing final verification
    # runs after this job and checks all masks, family outputs and score ranges.
    mod.attach_masks()
    for method in METHODS:
        isolated = scope(original(), ROOT, method)
        if mod.digest(isolated.S / method / 'library.fa') != mod.digest(mod.S / method / 'library.fa'):
            raise RuntimeError('Scored library differs from masking snapshot: ' + method)
        src = isolated.S / method / 'famlevel'
        dst = mod.S / method / 'famlevel'
        if dst.exists() or dst.is_symlink():
            if dst.resolve() != src.resolve():
                raise RuntimeError('Existing different family result: ' + method)
        else:
            dst.symlink_to(src, target_is_directory=True)
    (mod.E / 'results').mkdir(exist_ok=True)
    for suffix, (fields, rows) in collected.items():
        target = mod.E / 'results' / ('Zea_mays_' + suffix + '.csv')
        tmp = target.with_name(target.name + '.parallel.tmp')
        with tmp.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, target)
    mod.validate_base()
    mod.validate_family()
    print('Published complete results for all six methods', flush=True)


def index():
    import pyfastx
    genome = pyfastx.Fasta(str(original().G), full_index=True)
    print('Genome bases:', genome.size, 'composition:', genome.composition, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'index', 'base', 'family', 'merge'])
    parser.add_argument('--method', choices=METHODS)
    parser.add_argument('--task', type=int, choices=range(6))
    args = parser.parse_args()
    if args.stage in ('base', 'family'):
        method = args.method or (METHODS[args.task] if args.task is not None else None)
        if method is None:
            parser.error('method or task is required')
        run(args.stage, method)
    else:
        {'prepare': prepare, 'index': index, 'merge': merge}[args.stage]()
