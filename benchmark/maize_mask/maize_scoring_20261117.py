"""Finish maize scoring with the existing benchmark algorithms and strict input gates."""
import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

B = Path('/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark').resolve()
M = B / 'maize_mask'
E = B / 'edta_manual_bench'
S = E / 'zea_mays'
G = B / 'genomes/plant/zea_mays/zea_mays.fna'
PY = '/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python'
METHODS = [line.rstrip().split('\t') for line in (M / 'methods.tsv').read_text().splitlines() if line.strip() and not line.startswith('#')]
TESTS = [m for m, _ in METHODS if m != 'reference']
LABEL = {m: ('edta' if m == 'edta2' else m) for m in TESTS}
R = M / 'scoring_20261117'

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def table(path):
    with open(path) as f:
        return list(csv.DictReader(f))

def atomic_json(path, data):
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    os.replace(tmp, path)

def register():
    R.mkdir(exist_ok=True)
    with open(R / 'register.lock', 'a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        records = {}
        for method, original in METHODS:
            snap = M / 'libs' / (method + '.fa')
            src = snap if snap.exists() else Path(original).resolve()
            require(src.is_file() and src.stat().st_size > 0, 'Library not ready: ' + method)
            sha = digest(src)
            dst = S / method / 'library.fa'
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                require(digest(dst) == sha, 'Registered library differs from masking input: ' + method)
            else:
                tmp = dst.with_name('library.fa.tmp.' + str(os.getpid()))
                shutil.copyfile(src, tmp)
                os.replace(tmp, dst)
            records[method] = {'source': str(src), 'original': original, 'sha256': sha}
        atomic_json(R / 'libraries.json', records)
    return records

def validate_masks():
    chunks = [line.split('\t')[0] for line in (M / 'chunks.tsv').read_text().splitlines()[1:] if line.strip()]
    for method, _ in METHODS:
        absent = [c for c in chunks if not (M / 'queue/done' / (method + '__' + c)).is_file()]
        require(not absent, '%s masking incomplete: %d/%d chunks missing' % (method, len(absent), len(chunks)))
        src = M / 'merged' / method / 'zea_mays.fna.out'
        require(src.is_file() and src.stat().st_size > 0, 'Merged masking output not ready: ' + method)
        require((M / 'libs' / (method + '.fa')).is_file(), 'Masking library snapshot missing: ' + method)
    print('All seven merged masks and every chunk completion marker are present', flush=True)

def attach_masks():
    chunks = [line.split('\t')[0] for line in (M / 'chunks.tsv').read_text().splitlines()[1:] if line.strip()]
    for method, _ in METHODS:
        require(all((M / 'queue/done' / (method + '__' + c)).is_file() for c in chunks), 'Cannot score incomplete chunks: ' + method)
    # The last worker can still be merging after the feeder observes every done marker.
    deadline = time.monotonic() + 7200
    while not all((M / 'merged' / method / 'zea_mays.fna.out').is_file() for method, _ in METHODS):
        require(time.monotonic() < deadline, 'Timed out waiting for final atomic mask merges')
        print('Waiting for final atomic mask merges', flush=True)
        time.sleep(30)
    validate_masks()
    register()
    for method, _ in METHODS:
        src = M / 'merged' / method / 'zea_mays.fna.out'
        dst = S / method / 'zea_mays.fna.out'
        if dst.exists():
            require(digest(src) == digest(dst), 'Registered mask differs: ' + method)
        else:
            try:
                os.link(src, dst)
            except OSError:
                shutil.copyfile(src, dst)

def validate_base():
    for suffix, needed in [('bench_edta_manual', ['any']), ('bench_harmonised', ['any', 'any_dontcare'])]:
        path = E / 'results' / ('Zea_mays_' + suffix + '.csv')
        rows = table(path)
        for label in LABEL.values():
            for kind in needed:
                found = [r for r in rows if r['Method'] == label and r['TE type'] == kind]
                require(len(found) == 1, 'Missing/duplicate score: %s %s %s' % (suffix, label, kind))
                require(0 <= float(found[0]['F1']) <= 1, 'Invalid F1: ' + label)

def validate_family():
    path = E / 'results/Zea_mays_family_level.csv'
    rows = table(path)
    for method in TESTS:
        found = [r for r in rows if r['Method'] == method and r['class'] == 'all']
        require(len(found) == 1 and int(found[0]['n_reference']) > 0, 'Family result incomplete: ' + method)
        require((S / method / 'famlevel/library.fa.out').is_file(), 'Family masking output missing: ' + method)

def family():
    register()
    subprocess.run(['bash', '-euo', 'pipefail', str(E / 'family_level.sh'), 'zea_mays'], check=True)
    validate_family()
    atomic_json(R / 'family_complete.json', {'completed_at': datetime.datetime.now().isoformat()})

def base():
    attach_masks()
    rows = []
    for method in TESTS:
        d = S / method
        out = d / 'zea_mays.fna.out'
        with open(d / 'libtest.log', 'w') as log:
            subprocess.run([PY, str(B / 'tools/EDTA/lib-test.py'), '--genome', str(G), '--reference', str(S / 'reference/zea_mays.fna.out'), '--test', str(out), '--extended_report'], cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
        report = d / 'zea_mays.fna.out.repeatmasker.lib.report'
        for line in report.read_text().splitlines():
            fields = line.split('\t')
            match = re.search(r'\.(\w+)\.lib\.report$', fields[0])
            if match and len(fields) >= 11:
                rows.append([LABEL[method], match.group(1)] + fields[1:11])
        print('Default lib-test complete: ' + method, flush=True)
    dest = E / 'results/Zea_mays_bench_edta_manual.csv'
    tmp = dest.with_suffix('.csv.tmp')
    with tmp.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Method', 'TE type', 'sens', 'spec', 'accu', 'prec', 'FDR', 'F1', 'TP', 'TN', 'FP', 'FN'])
        writer.writerows(rows)
    os.replace(tmp, dest)
    subprocess.run([PY, str(E / 'libtest_harmonised.py'), 'zea_mays'] + TESTS, check=True)
    validate_base()
    atomic_json(R / 'base_complete.json', {'completed_at': datetime.datetime.now().isoformat()})

def finalize():
    validate_masks()
    validate_base()
    validate_family()
    scores = table(E / 'results/Zea_mays_bench_harmonised.csv')
    families = table(E / 'results/Zea_mays_family_level.csv')
    fields = ['method', 'harmonised_F1', 'sensitivity', 'precision', 'perfect', 'good', 'present_only', 'not_found', 'n_library', 'n_reference']
    with (R / 'maize_summary.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for method in TESTS:
            sc = next(r for r in scores if r['Method'] == LABEL[method] and r['TE type'] == 'any_dontcare')
            fa = next(r for r in families if r['Method'] == method and r['class'] == 'all')
            writer.writerow(dict(method=method, harmonised_F1=sc['F1'], sensitivity=sc['sens'], precision=sc['prec'], **{k: fa[k] for k in fields[4:]}))
    files = [E / 'results' / ('Zea_mays_' + suffix + '.csv') for suffix in ('bench_edta_manual', 'bench_harmonised', 'family_level')]
    atomic_json(R / 'COMPLETE.json', {'completed_at': datetime.datetime.now().isoformat(), 'results_sha256': {str(p): digest(p) for p in files}})
    print((R / 'maize_summary.tsv').read_text(), flush=True)

if __name__ == '__main__':
    mode = argparse.ArgumentParser()
    mode.add_argument('stage', choices=['verify-masks', 'register', 'family', 'base', 'finalize'])
    args = mode.parse_args()
    {'verify-masks': validate_masks, 'register': register, 'family': family, 'base': base, 'finalize': finalize}[args.stage]()
