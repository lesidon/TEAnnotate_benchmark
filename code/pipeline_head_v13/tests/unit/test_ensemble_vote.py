"""Two-step ensemble vote (bin/ensemble_vote.py): Order first, then a superfamily
inside the winning Order."""
from conftest import load_script

ev = load_script('ensemble_vote')


def cfg(kingdom):
    return ev.load_two_step(kingdom)


def test_superfamily_comes_from_the_winning_order():
    # Fungus: LTR wins the Order step (create 1 + neuralte 1.5 > terrier 2). Terrier's
    # hAT call is outside LTR and must not decide the superfamily.
    calls = {'create': 'LTR/Gypsy', 'neuralte': 'LTR/Copia', 'terrier': 'DNA/hAT'}
    assert ev.decide_two_step(calls, cfg('fungus'), ['create', 'neuralte', 'terrier']) == 'LTR/Copia'


def test_order_tie_is_no_call():
    calls = {'neuralte': 'LINE/L1', 'terrier': 'DNA/hAT'}   # plant: 2 vs 2
    assert ev.decide_two_step(calls, cfg('plant'), ['teclass2', 'create', 'neuralte', 'terrier']) is None


def test_order_gate_counts_abstaining_models():
    # Only CREATE (weight 1 of 4.5) votes; 1 < 0.3 * 4.5, so the consensus stays Unknown.
    calls = {'create': 'LTR/Gypsy', 'neuralte': 'Unknown', 'terrier': 'Unknown'}
    assert ev.decide_two_step(calls, cfg('fungus'), ['create', 'neuralte', 'terrier']) is None


def test_crashed_model_does_not_count_against_the_gate():
    # NeuralTE and Terrier crashed (no vote file): CREATE alone carries the Order step.
    # It has no fungal superfamily weight, so the call stops at the Order.
    calls = {'create': 'LTR/Gypsy'}
    assert ev.decide_two_step(calls, cfg('fungus'), ['create']) == 'LTR'


def test_l2_and_cr1_pool_their_superfamily_votes():
    calls = {'neuralte': 'LINE/L2', 'terrier': 'LINE/CR1'}
    assert ev.decide_two_step(calls, cfg('animal'), ['teclass2', 'neuralte', 'terrier']) == 'LINE/CR1'


def test_dirs_and_penelope_are_orders_of_their_own():
    assert ev.vote_order('LTR/DIRS') == 'DIRS'
    assert ev.vote_order('LINE/Penelope') == 'PLE'
    assert ev.vote_order('LTR/Gypsy') == 'LTR'


def test_kingdom_aware_subset_drops_other_models():
    c = ev.load_two_step('plant', include={'neuralte', 'terrier'})
    assert set(c['order']['weights']) == {'neuralte', 'terrier'}
    assert set(c['superfamily']['weights']) == {'neuralte', 'terrier'}


def test_unknown_kingdom_uses_animal_settings():
    assert ev.load_two_step('protist')['order']['weights'] == cfg('animal')['order']['weights']
