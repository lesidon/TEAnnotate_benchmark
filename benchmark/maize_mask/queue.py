#!/usr/bin/env python3
"""Work queue of the chunked maize masking (benchmark only).

An item is one <method, chunk> pair: RepeatMasker of one genome chunk with one library. State is
kept in files under queue/ (claims/, done/, attempts/), and claims are made inside a mkdir mutex,
so workers on different nodes never take the same item. A claim whose job is no longer in squeue
is stale and the item is retried, at most MAX_ATTEMPTS times. Methods are taken in the order of
methods.tsv, chunks largest first. A library is snapshotted into libs/ when its first item is
claimed, and only once its file is at least READY_AGE seconds old (still being written otherwise).

  queue.py claim <job> <cpus> <seconds_left>  -> "method chunk library chunk_fasta est_s", or nothing
  queue.py done <method> <chunk> <seconds> <cpus>
  queue.py fail <method> <chunk>
  queue.py claimable     items a new worker could start now (library ready, unclaimed, attempts left)
  queue.py remaining     items not done and not exhausted (including methods whose library is not ready)
  queue.py merge_ready   methods whose chunks are all done and that are not merged yet
  queue.py status
"""
import getpass
import hashlib
import os
import shutil
import subprocess
import sys
import time

M = os.path.dirname(os.path.abspath(__file__))
Q = os.path.join(M, 'queue')
MAX_ATTEMPTS = 3
DEFAULT_RATE = 115.0     # CPU-h per Mb of library per Gb of genome (Danio Earl Grey masking, 64 CPUs)
SAFETY, MARGIN = 1.3, 900
READY_AGE = 600
for d in ('claims', 'done', 'attempts'):
    os.makedirs(os.path.join(Q, d), exist_ok=True)
os.makedirs(os.path.join(M, 'libs'), exist_ok=True)


def p(*a):
    return os.path.join(Q, *a)


class Mutex:
    def __enter__(self):
        self.path = p('.mutex')
        while True:
            try:
                os.mkdir(self.path)
                return self
            except FileExistsError:
                try:
                    if time.time() - os.stat(self.path).st_mtime > 300:
                        os.rmdir(self.path)
                        continue
                except FileNotFoundError:
                    continue
                time.sleep(0.5 + (os.getpid() % 10) / 10)

    def __exit__(self, *a):
        os.rmdir(self.path)


def methods():
    out = []
    for line in open(os.path.join(M, 'methods.tsv')):
        if line.strip() and not line.startswith('#'):
            m, lib = line.rstrip('\n').split('\t')[:2]
            out.append((m, lib))
    return out


def chunks():
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(M, 'chunks.tsv'))][1:]
    return sorted(((r[0], int(r[1])) for r in rows), key=lambda x: (-x[1], x[0]))


def live_jobs():
    try:
        r = subprocess.run(['squeue', '-h', '-u', getpass.getuser(), '-o', '%i'], capture_output=True,
                           text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return set(r.stdout.split()) if r.returncode == 0 else None   # None: treat every claim as live


def ready(src):
    return os.path.isfile(src) and os.path.getsize(src) > 0 and time.time() - os.path.getmtime(src) > READY_AGE


def snapshot(m, src):
    dst = os.path.join(M, 'libs', m + '.fa')
    if os.path.exists(dst):
        return dst
    if not ready(src):
        return None
    tmp = '%s.tmp%d' % (dst, os.getpid())
    shutil.copyfile(src, tmp)
    with open(tmp, 'rb') as fh:
        md5 = hashlib.md5(fh.read()).hexdigest()
    os.replace(tmp, dst)
    with open(os.path.join(M, 'libs', m + '.src'), 'w') as fh:
        fh.write('%s\t%s\t%s\n' % (src, md5, time.ctime(os.path.getmtime(src))))
    return dst


def attempts(k):
    f = p('attempts', k)
    return sum(1 for _ in open(f)) if os.path.exists(f) else 0


def done_rows(m=None):
    out = []
    for k in os.listdir(p('done')):
        mm, c = k.split('__')
        if m is None or mm == m:
            sec, cpus, bp = open(p('done', k)).read().split()[:3]
            out.append((mm, c, float(sec), int(cpus), int(bp)))
    return out


def rate(m, lib):
    """CPU-seconds per genome bp for this method: observed if any chunk is done, else the default."""
    d = done_rows(m)
    if d:
        return sum(s * c for _, _, s, c, _ in d) / sum(bp for *_, bp in d)
    return DEFAULT_RATE * 3600 * (os.path.getsize(lib) / 1e6) / 1e9


def is_live_claim(k, live):
    f = p('claims', k)
    if not os.path.exists(f):
        return False
    owner = open(f).read().split()[0]
    return live is None or owner in live


def claim(job, cpus, left):
    with Mutex():
        live = live_jobs()
        for m, src in methods():
            lib = None
            for c, bp in chunks():
                k = '%s__%s' % (m, c)
                if os.path.exists(p('done', k)) or is_live_claim(k, live):
                    continue
                if attempts(k) >= MAX_ATTEMPTS:
                    continue
                lib = lib or snapshot(m, src)
                if lib is None:
                    break                          # library not ready yet: next method
                est = rate(m, lib) * bp / cpus
                if est * SAFETY + MARGIN > left:
                    continue
                with open(p('claims', k), 'w') as fh:
                    fh.write('%s %d\n' % (job, time.time()))
                with open(p('attempts', k), 'a') as fh:
                    fh.write('%s %d\n' % (job, time.time()))
                print(m, c, lib, os.path.join(M, 'chunks', c + '.fa'), int(est))
                return


def finish(m, c, sec=None, cpus=None):
    k = '%s__%s' % (m, c)
    with Mutex():
        if sec is not None:
            bp = dict(chunks())[c]
            with open(p('done', k), 'w') as fh:
                fh.write('%s %s %d\n' % (sec, cpus, bp))
        if os.path.exists(p('claims', k)):
            os.remove(p('claims', k))


def count(need_ready):
    live = live_jobs()
    n = 0
    for m, src in methods():
        if need_ready and not (os.path.exists(os.path.join(M, 'libs', m + '.fa')) or ready(src)):
            continue
        for c, _ in chunks():
            k = '%s__%s' % (m, c)
            if os.path.exists(p('done', k)) or attempts(k) >= MAX_ATTEMPTS and not is_live_claim(k, live):
                continue
            if need_ready and is_live_claim(k, live):
                continue
            n += 1
    return n


def merge_ready():
    cs = [c for c, _ in chunks()]
    for m, _ in methods():
        if os.path.exists(os.path.join(M, 'merged', m, 'zea_mays.fna.out')):
            continue
        if all(os.path.exists(p('done', '%s__%s' % (m, c))) for c in cs):
            print(m)


def status():
    live = live_jobs()
    cs = chunks()
    gb = sum(bp for _, bp in cs) / 1e9
    print('%-22s %6s %7s %7s %9s %10s %s' % ('method', 'done', 'running', 'failed', 'CPU-h', 'CPU-h left', 'merged'))
    for m, src in methods():
        d = done_rows(m)
        run = sum(1 for c, _ in cs if is_live_claim('%s__%s' % (m, c), live))
        bad = sum(1 for c, _ in cs if not os.path.exists(p('done', '%s__%s' % (m, c)))
                  and attempts('%s__%s' % (m, c)) >= MAX_ATTEMPTS and not is_live_claim('%s__%s' % (m, c), live))
        used = sum(s * c for _, _, s, c, _ in d) / 3600
        lib = os.path.join(M, 'libs', m + '.fa')
        lib = lib if os.path.exists(lib) else (src if os.path.exists(src) else None)
        left = '-'
        if lib:
            left = '%.0f' % (rate(m, lib) * (gb * 1e9 - sum(bp for *_, bp in d)) / 3600)
        merged = 'yes' if os.path.exists(os.path.join(M, 'merged', m, 'zea_mays.fna.out')) else ''
        print('%-22s %3d/%-3d %7d %7d %9.0f %10s %s%s' % (m, len(d), len(cs), run, bad, used, left, merged,
                                                     '' if lib else '  (library not ready)'))


if __name__ == '__main__':
    cmd, a = sys.argv[1], sys.argv[2:]
    if cmd == 'claim':
        claim(a[0], int(a[1]), float(a[2]))
    elif cmd == 'done':
        finish(a[0], a[1], a[2], a[3])
    elif cmd == 'fail':
        finish(a[0], a[1])
    elif cmd == 'claimable':
        print(count(True))
    elif cmd == 'remaining':
        print(count(False))
    elif cmd == 'merge_ready':
        merge_ready()
    elif cmd == 'status':
        status()
    else:
        sys.exit('unknown command ' + cmd)
