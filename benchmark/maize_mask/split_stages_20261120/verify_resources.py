"""Check the resolved recovery config before submitting any scientific stages."""
from pathlib import Path
import sys


def verify(path):
    cfg = {}
    for line in Path(path).read_text().splitlines():
        if ' = ' in line:
            key, value = line.split(' = ', 1)
            cfg[key] = value.strip("'")
    assert cfg['process.executor'] == 'slurm', 'Stages must use separate Slurm jobs'
    assert cfg['process.queue'] == 'htc'
    stages = {
        'CURATION:TETRIMMER_POSTPROCESS': (4, None, '2d', 'mem'),
        'CURATION:CD_HIT_EST_CURATED': (16, '32 GB', '1d', 'htc'),
        'CURATION:DEDUP_REFINE_CURATED': (16, '64 GB', '3d', 'mem'),
        'CURATION:DOMAIN_RELABEL': (4, '16 GB', '1d', 'htc'),
        'CURATION:EXTRACT_UNKNOWN': (1, '4 GB', '4h', 'htc'),
        'CURATION:TECLASS2': (8, '64 GB', '1d', 'htc'),
        'CURATION:DEEPTE': (4, '16 GB', '1d', 'htc'),
        'CURATION:CREATE': (8, '32 GB', '1d', 'htc'),
        'CURATION:NEURALTE': (8, '32 GB', '1d', 'htc'),
        'CURATION:TERRIER': (4, '16 GB', '1d', 'htc'),
        'CURATION:ENSEMBLE_VOTE': (1, '8 GB', '4h', 'htc'),
        'FINAL_LIBRARY': (1, '4 GB', '4h', 'htc'),
    }
    for stage, (cpus, memory, duration, queue) in stages.items():
        prefix = "process.'withName:" + stage + "'."
        assert int(cfg[prefix + 'cpus']) == cpus, stage
        if memory:
            assert cfg[prefix + 'memory'] == memory, stage
        assert cfg[prefix + 'time'] == duration, stage
        assert cfg.get(prefix + 'queue', cfg['process.queue']) == queue, stage
        if queue == 'htc':
            assert cpus <= 16 and int(memory.split()[0]) <= 64
            assert duration in ('1d', '4h')
    pp = "process.'withName:CURATION:TETRIMMER_POSTPROCESS'."
    assert cfg[pp + 'memory'] == '{ 32.GB * task.attempt }'
    assert cfg[pp + 'maxRetries'] == '1'
    assert '137' in cfg[pp + 'errorStrategy']
    print('Verified twelve stage allocations, HTC limits, and postprocess32-to64GB OOM retry')


if __name__ == '__main__':
    verify(sys.argv[1])
