"""Count existing queue items that fit a worker allocation, without claiming any."""
import importlib.util
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('maize_queue', str(Path(__file__).with_name('queue.py')))
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)
cpus, seconds = int(sys.argv[1]), int(sys.argv[2])
live, count = q.live_jobs(), 0
for method, source in q.methods():
    snap = Path(q.M) / 'libs' / (method + '.fa')
    if snap.is_file():
        lib = str(snap)
    elif q.ready(source):
        lib = source
    else:
        continue
    for chunk, bp in q.chunks():
        key = method + '__' + chunk
        if Path(q.p('done', key)).exists() or q.is_live_claim(key, live) or q.attempts(key) >= q.MAX_ATTEMPTS:
            continue
        if q.rate(method, lib) * bp / cpus * q.SAFETY + q.MARGIN <= seconds:
            count += 1
print(count)
