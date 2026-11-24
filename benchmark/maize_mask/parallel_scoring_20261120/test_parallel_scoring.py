import importlib.util
from pathlib import Path
import tempfile
import types
import unittest
import csv
import subprocess
import sys
from unittest import mock

SPEC = importlib.util.spec_from_file_location('parallel', Path(__file__).with_name('parallel_scoring.py'))
p = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(p)


class Tests(unittest.TestCase):
    def test_scope_selects_only_reference_and_requested_method(self):
        with tempfile.TemporaryDirectory() as tmp:
            mod = types.SimpleNamespace(METHODS=[('reference', 'ref'), ('hite', 'h'), ('heliano', 'x')])
            p.scope(mod, Path(tmp), 'hite')
            self.assertEqual(mod.METHODS, [('reference', 'ref'), ('hite', 'h')])
            self.assertEqual(mod.TESTS, ['hite'])
            self.assertEqual(mod.S, Path(tmp) / 'hite/benchmark/zea_mays')

    def test_scope_refuses_unknown_method(self):
        with self.assertRaises(ValueError):
            p.scope(types.SimpleNamespace(METHODS=[('reference', 'r')]), Path('.'), 'missing')

    def test_merge_refuses_missing_or_duplicate_method_rows(self):
        rows = [{'Method': 'hite', 'TE type': 'any'}]
        with self.assertRaises(RuntimeError):
            p.validate_rows(rows, ['hite', 'heliano'], 'TE type', ['any'])
        with self.assertRaises(RuntimeError):
            p.validate_rows(rows * 2, ['hite'], 'TE type', ['any'])

    def test_isolated_harmoniser_changes_only_destination_assignment(self):
        source = 'E = f"{B}/edta_manual_bench"\ndef canon(x): return x\n'
        patched = p.redirect_harmoniser(source)
        self.assertIn('MAIZE_SCORE_ROOT', patched)
        self.assertEqual(patched.splitlines()[1:], source.splitlines()[1:])

    def test_refuses_unrecognised_harmoniser_source(self):
        with self.assertRaises(RuntimeError):
            p.redirect_harmoniser('unexpected source')

    @unittest.skipUnless((p.B / 'tools/EDTA/lib-test.py').exists(), 'cluster integration fixture')
    def test_real_default_and_harmonised_scores_on_tiny_identical_masks(self):
        with tempfile.TemporaryDirectory(dir=p.ROOT) as tmp:
            root = Path(tmp)
            mod = p.scope(p.original(), root, 'hite')
            mod.E.joinpath('results').mkdir(parents=True)
            mod.R.mkdir(parents=True)
            mod.M = root / 'masks'
            mod.G = root / 'genome.fa'
            mod.G.write_text('>chr1\n' + 'ACGT' * 250 + '\n')
            (mod.M / 'queue/done').mkdir(parents=True)
            (mod.M / 'libs').mkdir()
            (mod.M / 'chunks.tsv').write_text('chunk\nchunk0\n')
            for method in ('reference', 'hite'):
                (mod.M / 'queue/done' / (method + '__chunk0')).touch()
                (mod.M / 'libs' / (method + '.fa')).write_text('>repeat#LTR/Gypsy\n' + 'ACGT' * 25 + '\n')
                d = mod.M / 'merged' / method
                d.mkdir(parents=True)
                (d / 'zea_mays.fna.out').write_text('header\nheader\n\n 1000 0.0 0.0 0.0 chr1 1 100 (900) + repeat LTR/Gypsy 1 100 (0) 1\n')
            hspec = importlib.util.spec_from_file_location('fixture_harmoniser', p.B / 'edta_manual_bench/libtest_harmonised.py')
            h = importlib.util.module_from_spec(hspec)
            hspec.loader.exec_module(h)
            h.E = str(mod.E)
            h.sizes = lambda sp: (str(mod.G), {'chr1': 1000})
            real_run = subprocess.run

            def invoke(command, **kwargs):
                if str(command[1]).endswith('libtest_harmonised.py'):
                    with mock.patch.object(sys, 'argv', ['harmoniser', 'zea_mays', 'hite']):
                        h.main()
                    return subprocess.CompletedProcess(command, 0)
                return real_run(command, **kwargs)

            with mock.patch.object(subprocess, 'run', side_effect=invoke):
                mod.base()
            for suffix in ('bench_edta_manual', 'bench_harmonised'):
                with (mod.E / 'results' / ('Zea_mays_' + suffix + '.csv')).open() as f:
                    rows = list(csv.DictReader(f))
                any_row = next(r for r in rows if r['TE type'] == 'any')
                self.assertEqual(float(any_row['F1']), 1.0)
                self.assertEqual({r['Method'] for r in rows}, {'hite'})


if __name__ == '__main__':
    unittest.main()
