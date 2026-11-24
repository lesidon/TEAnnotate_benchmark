"""Signal readiness only after all seven masks and immutable libraries validate."""
import datetime
import importlib.util
from pathlib import Path
import time


def wait_for_masks(scoring, timeout=5 * 86400 + 23 * 3600, interval=300, clock=time.monotonic, sleeper=time.sleep):
    deadline = clock() + timeout
    previous = None
    while True:
        try:
            scoring.validate_masks()
        except RuntimeError as error:
            remaining = deadline - clock()
            if remaining <= 0:
                raise TimeoutError('Maize masks did not become ready: ' + str(error)) from error
            message = str(error)
            if message != previous:
                print(datetime.datetime.now(datetime.timezone.utc).isoformat(), 'Waiting:', message, flush=True)
                previous = message
            sleeper(min(interval, remaining))
        else:
            # A library mismatch is a real failure, not something waiting can fix.
            libraries = scoring.register()
            scoring.atomic_json(scoring.R / 'MASKS_READY_20261118.json', {
                'status': 'ALL_MASKS_VERIFIED',
                'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'libraries': libraries,
            })
            print('All maize masking outputs and library identities verified', flush=True)
            return


if __name__ == '__main__':
    root = Path(__file__).resolve().parent.parent
    spec = importlib.util.spec_from_file_location('maize_scoring', root / 'maize_scoring_20261117.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    expected = {'reference', 'earlgrey', 'hite', 'teannotate_wo_v13', 'teannotate_v13', 'edta2', 'heliano'}
    if {method for method, _ in module.METHODS} != expected:
        raise RuntimeError('Unexpected maize method inventory; refusing partial scoring')
    wait_for_masks(module)
