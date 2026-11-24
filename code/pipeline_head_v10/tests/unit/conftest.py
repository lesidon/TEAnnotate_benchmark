"""Shared helpers for the unit tests of the pipeline's bin/ scripts.

Run from the repository root:  python -m pytest tests/unit
"""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / 'bin'


def load_script(name):
    """Import bin/<name>.py as a module (the scripts are not a package)."""
    spec = importlib.util.spec_from_file_location(name, BIN / f'{name}.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def run_script():
    """Run bin/<name>.py with arguments; fail the test on a non-zero exit."""
    def run(name, *args):
        r = subprocess.run([sys.executable, str(BIN / f'{name}.py'), *map(str, args)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        return r
    return run


def headers(path):
    """FASTA header lines of a file, without '>'."""
    return [l[1:].strip() for l in open(path) if l.startswith('>')]
