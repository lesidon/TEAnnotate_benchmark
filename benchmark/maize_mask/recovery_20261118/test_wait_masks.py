import unittest
from pathlib import Path
from types import SimpleNamespace
from wait_masks import wait_for_masks


class WaitTests(unittest.TestCase):
    def fake(self, failures=1, registration_error=False):
        state = dict(checks=0, registrations=0, written=None, now=0)
        def check():
            state['checks'] += 1
            if state['checks'] <= failures:
                raise RuntimeError('teannotate_v13 masking incomplete')
        def register():
            state['registrations'] += 1
            if registration_error:
                raise RuntimeError('Registered library differs from masking input')
            return {'reference': {'sha256': 'verified'}}
        def write(path, data):
            state['written'] = data
        scoring = SimpleNamespace(validate_masks=check, register=register, atomic_json=write, R=Path('/fixture'))
        def sleep(seconds):
            state['now'] += seconds
        return scoring, state, sleep

    def test_waits_for_every_mask_before_registering_or_signalling(self):
        scoring, state, sleep = self.fake()
        wait_for_masks(scoring, timeout=100, interval=10, clock=lambda: state['now'], sleeper=sleep)
        self.assertEqual(state['checks'], 2)
        self.assertEqual(state['registrations'], 1)
        self.assertEqual(state['written']['status'], 'ALL_MASKS_VERIFIED')

    def test_library_hash_mismatch_is_fatal_not_retried(self):
        scoring, state, sleep = self.fake(failures=0, registration_error=True)
        with self.assertRaisesRegex(RuntimeError, 'differs'):
            wait_for_masks(scoring, timeout=100, interval=10, clock=lambda: state['now'], sleeper=sleep)
        self.assertEqual(state['checks'], 1)
        self.assertIsNone(state['written'])

    def test_incomplete_masks_never_emit_ready_marker(self):
        scoring, state, sleep = self.fake(failures=100)
        with self.assertRaises(TimeoutError):
            wait_for_masks(scoring, timeout=20, interval=10, clock=lambda: state['now'], sleeper=sleep)
        self.assertEqual(state['registrations'], 0)
        self.assertIsNone(state['written'])


if __name__ == '__main__':
    unittest.main()
