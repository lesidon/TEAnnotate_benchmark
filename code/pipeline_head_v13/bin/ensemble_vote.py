#!/usr/bin/env python3
"""Ensemble for TE re-classification (hard / weighted voting).

Each classifier (TEclass2, DeepTE, CREATE, NeuralTE, Terrier) contributes a
[base_id <TAB> RM_label] TSV of votes (one vote per classifier per consensus; a
classifier that abstains simply has no row). For each consensus whose current
library label is Unknown, the winning class is chosen by one of:

  two_step  TWO-STEP VOTE (default). Step 1 votes on the Order with per-kingdom
            model weights and an abstention gate. Step 2 picks a superfamily INSIDE
            the winning Order with its own per-model weights: only models whose call
            lies in that Order take part, so a superfamily of another Order can never
            win; if no superfamily clears its gate the call stays at Order level.
            Orders follow the classification rather than the RepeatMasker prefix
            (DIRS/Ngaro/VIPER vote for DIRS, Penelope for PLE). The per-organism_type
            settings (TWO_STEP_DEFAULTS) were tuned on leakage-filtered curated
            consensi (PanTEon v1.5.1 + MycoMobilome v1.1 held out from RepBase/Dfam,
            cluster-split train set, support-weighted F1) and can be overridden with
            --two-step-config. The abstention gate and min_agree of the options below
            do not apply to it.

  hard      MAJORITY RULE. Vote at the Order level; the Order with the most votes
            wins iff it clears the abstention gate AND there is no tie. Refine to
            Order/Superfamily when the winning-order votes agree on one.

  weighted  WEIGHTED VOTE (default). Same logic, but each model contributes its
            weight instead of 1. Weights prioritise the strongest models — by
            default NeuralTE and Terrier, the top-2 at the Order level in our
            confusion-matrix analysis. A superfamily is adopted if it carries
            >= half the winning Order's weight. Weights can be set explicitly
            (--weights) or derived from an order-level metrics table
            (--weights-from).

  all       Emit ALL ensemble calls as separate voting-table columns
            (Ensemble_hard, Ensemble_weighted, Ensemble_two_step) for side-by-side
            comparison; the relabelled FASTA uses the two-step call.

Only calls whose Order is in --relabel-orders (default DNA,RC) relabel a consensus;
any other call leaves it Unknown. The Unknown fraction that reaches the ensemble has
already been through domain-based labelling, so it is mostly consensi without a
detectable coding domain: expected for non-autonomous DNA transposons (MITEs and other
TIR elements) and Helitrons, but for an "LTR" or "LINE" call usually a fragment or a
non-TE repeat that the classifiers, which have no "not a TE" answer, push into their
majority class. On the development genomes, Unknown -> LTR relabels fell almost
entirely outside reference LTRs (LTR F1 Zymoseptoria -0.04, Septoria -0.03) while
DNA / Helitron relabels mostly hit their class (Septoria Helitron F1 0 -> 0.88).
--relabel-orders all lifts the restriction.

For hard / weighted, the abstention gate is a FRACTION (--min-frac) of what the
voting subset could maximally cast, not an absolute number — see resolve_gates() for
why an absolute bar cannot work once the subset size varies by kingdom.

Confident (non-Unknown) labels are never overridden, and every input sequence is
kept. Soft voting is not used because DeepTE reports no per-class probabilities.

Output headers are RepeatModeler/RepeatMasker style ">name#Class/Subfamily" (the
winning RM label after '#'), so `RepeatMasker -lib` classifies the masked repeats
directly. Sequences that stay Unknown keep their original ">name#Unknown ..." header.
"""
import argparse
import copy
import csv
import json
import math
import os
import sys
from collections import Counter, defaultdict

UNKNOWN = {'', 'unknown', 'unspecified', 'unclassified', 'na', 'none'}

# Pretty model names for the voting-table columns, in a stable display order.
MODEL_DISPLAY = {
    'teclass2': 'TEclass2', 'deepte': 'DeepTE', 'create': 'CREATE',
    'neuralte': 'NeuralTE', 'terrier': 'Terrier',
}
MODEL_ORDER = ['teclass2', 'deepte', 'create', 'neuralte', 'terrier']

# Default per-model weights for the weighted ensemble. NeuralTE and Terrier are
# the strongest at the Order level (confusion-matrix analysis), so they lead;
# CREATE/TEclass2/DeepTE follow. Override with --weights or --weights-from.
DEFAULT_WEIGHTS = {
    'neuralte': 1.5, 'terrier': 1.5, 'create': 1.3,
    'teclass2': 1.2, 'deepte': 1.1,
}

# Two-step vote settings per organism_type: Order-step weights (a model absent here
# does not vote), gate as a fraction of the weight the models that ran could cast, and
# the minimum number of models on the winning Order; superfamily-step weights and gate
# (fraction of the superfamily weight of the models that voted for the winning Order;
# 0 = plurality). Tuned on the cluster-split train set of leakage-filtered PanTEon v1.5.1
# (+ MycoMobilome v1.1 for fungi) consensi for support-weighted F1. On the held-out test
# set (Order / superfamily weighted F1) against the best single model:
#   plant  0.963 / 0.941  (NeuralTE 0.943 / 0.937)
#   animal 0.910 / 0.890  (Terrier  0.873 / 0.877)
#   fungus 0.918 / 0.913  (NeuralTE 0.887 / 0.904)
# DeepTE was never selected. These model sets are also the default
# ensemble_models_<kingdom>, i.e. the classifiers that run.
TWO_STEP_DEFAULTS = {
    'plant': {
        'order':       {'weights': {'teclass2': 1.5, 'create': 1.0, 'neuralte': 2.0, 'terrier': 2.0},
                        'min_frac': 0.3, 'min_agree': 1},
        'superfamily': {'weights': {'teclass2': 1.0, 'create': 1.0, 'neuralte': 3.0, 'terrier': 1.0},
                        'min_frac': 0.5},
    },
    'animal': {
        'order':       {'weights': {'teclass2': 1.5, 'neuralte': 1.0, 'terrier': 2.0},
                        'min_frac': 0.3, 'min_agree': 1},
        'superfamily': {'weights': {'neuralte': 1.0, 'terrier': 2.0},
                        'min_frac': 0.5},
    },
    'fungus': {
        'order':       {'weights': {'create': 1.0, 'neuralte': 1.5, 'terrier': 2.0},
                        'min_frac': 0.3, 'min_agree': 1},
        'superfamily': {'weights': {'neuralte': 2.0, 'terrier': 1.0},
                        'min_frac': 0.0},
    },
}

# The Order a vote counts for: lineages RepeatMasker files under another prefix
# ('LTR/DIRS', 'LINE/Penelope') are Orders of their own, so a DIRS call never
# supports an LTR/Gypsy call.
SF_ORDER = {'dirs': 'DIRS', 'ngaro': 'DIRS', 'viper': 'DIRS', 'penelope': 'PLE'}

# Superfamilies the curated libraries behind the tuning treat as one (L2 within CR1,
# Jockey within I): their superfamily votes add up, and the output keeps the
# most-weighted label of the group.
SF_GROUP = {'l2': 'cr1', 'jockey': 'i'}

# Which ensemble columns each --method writes (column name, decision method).
ENSEMBLE_COLUMNS = {
    'two_step': [('Ensemble', 'two_step')],
    'hard':     [('Ensemble', 'hard')],
    'weighted': [('Ensemble', 'weighted')],
    'all':      [('Ensemble_hard', 'hard'), ('Ensemble_weighted', 'weighted'),
                 ('Ensemble_two_step', 'two_step')],
}
# Method used to relabel the FASTA (the library's Unknown fraction).
PRIMARY_METHOD = {'two_step': 'two_step', 'hard': 'hard', 'weighted': 'weighted',
                  'all': 'two_step'}


def model_from_path(path):
    """Derive the model name from a '<id>.<model>.preds.tsv' vote filename."""
    name = os.path.basename(path)
    for suffix in ('.preds.tsv', '.tsv'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    return name.rsplit('.', 1)[-1].lower()


def load_per_model(paths):
    """Return ({base_id: {model: RM_label}}, [models in encounter order])."""
    per_model = {}
    models = []
    for p in paths:
        model = model_from_path(p)
        if model not in models:
            models.append(model)
        with open(p) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 2 or not parts[0] or not parts[1]:
                    continue
                per_model.setdefault(parts[0], {})[model] = parts[1]
    return per_model, models


def order_models(models):
    """Canonical models first (TEclass2..Terrier), then any extras seen."""
    ordered = [m for m in MODEL_ORDER if m in models]
    ordered += [m for m in models if m not in ordered]
    return ordered


# --------------------------------------------------------------------------- #
# Weight configuration
# --------------------------------------------------------------------------- #
def parse_weights(spec):
    """Parse a 'name=w,name=w' string into {model_lower: float}."""
    out = {}
    if not spec:
        return out
    for tok in spec.replace(';', ',').split(','):
        tok = tok.strip()
        if '=' not in tok:
            continue
        k, v = tok.split('=', 1)
        try:
            out[k.strip().lower()] = float(v)
        except ValueError:
            sys.stderr.write("ensemble_vote: ignoring bad weight '%s'\n" % tok)
    return out


def _find_col(cols, candidates):
    lut = {str(c).lower(): c for c in (cols or [])}
    for cand in candidates:
        if cand.lower() in lut:
            return lut[cand.lower()]
    return None


def load_weights_from(path, metric):
    """Derive {model: weight} from an order-level metrics TSV.

    Reads evaluate_classification.py's '<tag>.all_models.summary_metrics.tsv'
    (or any TSV with a model column + the chosen metric column, e.g. macro_f1 /
    macro_precision / macro_recall) and normalises the metric to mean 1.0, so a
    better model gets proportionally more say while --min-weight stays comparable
    to a count of unit votes.
    """
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        cols = reader.fieldnames or []
        mcol = _find_col(cols, ['model', 'classifier', 'name'])
        vcol = _find_col(cols, [metric, 'macro_f1', 'f1',
                                'macro_precision', 'macro_recall'])
        if not mcol or not vcol:
            sys.stderr.write(
                "ensemble_vote: --weights-from needs model + metric columns; "
                "found %s\n" % cols)
            return {}
        for row in reader:
            try:
                val = float(row[vcol])
            except (ValueError, TypeError, KeyError):
                continue
            rows.append((str(row[mcol]).strip().lower(), val))
    vals = [v for _, v in rows if v > 0]
    if not vals:
        return {}
    mean = sum(vals) / len(vals)
    weights = {m: (v / mean if mean > 0 else 1.0) for m, v in rows if v > 0}
    sys.stderr.write("ensemble_vote: weights from %s (%s): %s\n" % (
        path, metric, ', '.join('%s=%.2f' % (m, w)
                                for m, w in sorted(weights.items()))))
    return weights


def resolve_weights(args):
    """Built-in defaults <- --weights-from <- --weights (later wins)."""
    weights = dict(DEFAULT_WEIGHTS)
    if args.weights_from:
        weights.update(load_weights_from(args.weights_from, args.weight_metric))
    if args.weights:
        weights.update(parse_weights(args.weights))
    return weights


def load_two_step(kingdom, path=None, include=None):
    """Two-step settings for one organism_type.

    TWO_STEP_DEFAULTS, overlaid with an optional JSON of the same shape (any subset
    of kingdoms / steps / keys; a 'weights' map replaces the default one). Unknown
    organism types use the animal settings, like the rest of the ensemble. With a
    kingdom-aware `include` set, models outside it do not vote in either step.
    """
    cfg = copy.deepcopy(TWO_STEP_DEFAULTS)
    if path:
        with open(path) as fh:
            user = json.load(fh)
        for k, steps in user.items():
            for step, vals in (steps or {}).items():
                cfg.setdefault(k, {}).setdefault(step, {}).update(vals)
    k = (kingdom or 'animal').lower()
    out = cfg.get(k) or cfg['animal']
    for step in ('order', 'superfamily'):
        w = {m.lower(): float(v) for m, v in out[step].get('weights', {}).items()}
        if include is not None:
            w = {m: v for m, v in w.items() if m in include}
        out[step]['weights'] = w
    return out


# --------------------------------------------------------------------------- #
# Abstention gate
# --------------------------------------------------------------------------- #
def eligible_models(models, include):
    """Models allowed to vote AND that actually ran.

    A model that crashed (errorStrategy 'ignore' -> no vote file) never reaches
    `models`, so it must not count against the gate: it had no chance to speak.
    A model that ran but abstained on a given repeat IS in `models` and does
    count, because declining to call a confident label is itself evidence of
    uncertainty about that repeat.
    """
    if include is None:
        return list(models)
    return [m for m in models if m.lower() in include]


def resolve_gates(args, models, weights, include):
    """Turn --min-frac into concrete (required_votes, required_weight) gates.

    The gate scales with what the voting subset could maximally cast rather than
    being an absolute number, because an absolute bar means different things to
    different subsets. Kingdom-aware voting (the default) hands a different model
    set to each organism_type, so a bar tuned for five models is nonsense for one:
    --min-weight 3.0 is a sane "most of the ensemble agrees" bar across all five
    (max 6.6), but the fungus subset is Terrier alone (max weight 1.5) and could
    never reach 3.0 — every fungal consensus would silently stay Unknown. That is
    why this used to disable the gate entirely for kingdom-aware runs, which left
    production with no abstention gate at all: an unopposed single model, or a
    three-way disagreement settled by 0.1 of weight, both produced a confident
    call.

    Expressing the bar as a fraction of the subset's own maximum fixes both ends
    and is near-identical to the old behaviour when all five models vote:
      all five  -> ceil(0.5*5) = 3 votes (the documented 3-of-5) and 0.5*6.6 = 3.3
      plant     -> ceil(0.5*4) = 2 votes and 0.5*5.1 = 2.55
      fungus    -> ceil(0.5*1) = 1 vote  and 0.5*1.5 = 0.75  (Terrier can call)

    --min-votes / --min-weight remain as optional ABSOLUTE floors (0 = off) on
    top of the fraction. A floor that the subset cannot reach is ignored with a
    warning rather than silently muting the whole stage.
    """
    elig = eligible_models(models, include)
    if not elig:
        return max(1, args.min_votes), args.min_weight

    max_w = sum(weights.get(m.lower(), 1.0) for m in elig)
    req_votes = max(1, int(math.ceil(args.min_frac * len(elig))))
    req_w = args.min_frac * max_w

    if args.min_votes > 0:
        if args.min_votes > len(elig):
            sys.stderr.write(
                "ensemble_vote: WARNING --min-votes %d exceeds the %d model(s) "
                "that can vote (%s); no consensus could ever be called, so the "
                "floor is ignored.\n"
                % (args.min_votes, len(elig), ",".join(sorted(elig))))
        else:
            req_votes = max(req_votes, args.min_votes)

    if args.min_weight > 0:
        if args.min_weight > max_w:
            sys.stderr.write(
                "ensemble_vote: WARNING --min-weight %.2f exceeds the maximum "
                "weight %.2f the voting subset (%s) can cast; no consensus "
                "could ever be called, so the floor is ignored.\n"
                % (args.min_weight, max_w, ",".join(sorted(elig))))
        else:
            req_w = max(req_w, args.min_weight)

    sys.stderr.write(
        "ensemble_vote: gate = %.2f of subset max (%s; %d model(s), max weight "
        "%.2f) -> needs >= %d vote(s) / >= %.2f weight\n"
        % (args.min_frac, ",".join(sorted(elig)), len(elig), max_w,
           req_votes, req_w))
    return req_votes, req_w


# --------------------------------------------------------------------------- #
# Decision rules
# --------------------------------------------------------------------------- #
def decide_hard(labels, min_votes):
    """Majority-rule winner (RM label) from a list of votes, or None."""
    if len(labels) < min_votes:
        return None
    orders = Counter(l.split('/', 1)[0] for l in labels)
    ranked = orders.most_common()
    top_order, top_n = ranked[0]
    if top_n < min_votes:
        return None
    if len(ranked) > 1 and ranked[1][1] == top_n:   # tie at the top -> no majority
        return None
    # Refine to a superfamily if the winning-order votes agree on one.
    sfs = Counter(l for l in labels if l.split('/', 1)[0] == top_order and '/' in l)
    if sfs:
        sf, sf_n = sfs.most_common(1)[0]
        if sf_n >= min_votes:
            return sf
    return top_order


def decide_weighted(model_labels, weights, min_weight):
    """Weighted-vote winner from {model: RM_label}, or None.

    Each model adds its weight to the Order it predicts. The top Order wins iff
    its summed weight >= min_weight and it is not tied. A superfamily is adopted
    when it carries >= half of the winning Order's weight.
    """
    order_w = defaultdict(float)
    sf_w = defaultdict(float)
    sf_order = {}
    for model, label in model_labels.items():
        if not label or is_unknown(label):
            continue
        w = weights.get(model.lower(), 1.0)
        order = label.split('/', 1)[0]
        order_w[order] += w
        if '/' in label:
            sf_w[label] += w
            sf_order[label] = order
    if not order_w:
        return None
    top_order, top_w = max(order_w.items(), key=lambda kv: kv[1])
    if top_w < min_weight:
        return None
    ranked_w = sorted(order_w.values(), reverse=True)
    if len(ranked_w) > 1 and abs(ranked_w[1] - top_w) < 1e-9:  # tie -> no call
        return None
    cand = {sf: w for sf, w in sf_w.items() if sf_order[sf] == top_order}
    if cand:
        top_sf, sf_wv = max(cand.items(), key=lambda kv: kv[1])
        if sf_wv >= 0.5 * top_w:
            return top_sf
    return top_order


def call_method(name, model_labels, min_votes, weights, min_weight, include=None):
    """Dispatch a single consensus's per-model calls to a decision rule.

    `include` (a set of lower-cased model names) restricts the ensemble decision
    to that kingdom-specific subset; models outside it keep their column in the
    table but do not vote.
    """
    if include is not None:
        model_labels = {m: l for m, l in model_labels.items() if m.lower() in include}
    if name == 'hard':
        return decide_hard(list(model_labels.values()), min_votes)
    if name == 'weighted':
        return decide_weighted(model_labels, weights, min_weight)
    raise ValueError("unknown ensemble method: %s" % name)


def call_with_agreement(name, model_labels, min_votes, weights, min_weight, include,
                        min_agree):
    """call_method, then require >= min_agree models of the subset on the winning Order.

    Unlike --min-votes, this floor is never relaxed for a small subset: a kingdom
    whose subset is a single model (fungus = Terrier) then relabels nothing. In the
    benchmark, single-model calls on Unknown consensi were mostly wrong: in Septoria,
    the seven Unknown -> LTR calls masked 0.50 Mb outside the reference against
    0.02 Mb inside it. Unknown costs a class only its recall; a wrong class costs
    precision in one class and recall in another.
    """
    winner = call_method(name, model_labels, min_votes, weights, min_weight, include)
    if winner and min_agree > 1:
        order = winner.split('/', 1)[0]
        agree = sum(1 for m, l in model_labels.items()
                    if (include is None or m.lower() in include)
                    and l and not is_unknown(l) and l.split('/', 1)[0] == order)
        if agree < min_agree:
            return None
    return winner


def vote_order(label):
    """The Order a call votes for ('LTR/DIRS' -> DIRS, 'LINE/Penelope' -> PLE)."""
    order, _, sf = label.partition('/')
    return SF_ORDER.get(sf.lower(), order) if sf else order


def sf_key(label):
    """Superfamily key of a call; SF_GROUP members share one key."""
    sf = label.partition('/')[2].lower()
    return SF_GROUP.get(sf, sf)


def decide_two_step(model_labels, cfg, ran):
    """Two-step winner (RM label) from {model: RM_label}, or None.

    `ran` lists the models that produced a vote file. A model that ran but
    abstained on this repeat counts against the Order gate; a crashed model does
    not (as in resolve_gates).
    """
    ow = cfg['order']['weights']
    sw = cfg['superfamily']['weights']
    max_w = sum(ow.get(m, 0.0) for m in ran)
    if max_w <= 0:
        return None
    calls = {m.lower(): l for m, l in model_labels.items() if l and not is_unknown(l)}

    # Step 1: Order.
    order_w = defaultdict(float)
    order_n = Counter()
    for m, label in calls.items():
        w = ow.get(m, 0.0)
        if w > 0:
            order = vote_order(label)
            order_w[order] += w
            order_n[order] += 1
    if not order_w:
        return None
    ranked = sorted(order_w.values(), reverse=True)
    top, top_w = max(order_w.items(), key=lambda kv: kv[1])
    if top_w < cfg['order']['min_frac'] * max_w - 1e-9:
        return None
    if len(ranked) > 1 and abs(ranked[1] - top_w) < 1e-9:      # tie -> no call
        return None
    if order_n[top] < cfg['order'].get('min_agree', 1):
        return None

    # Step 2: superfamily inside the winning Order.
    voted = 0.0
    key_w = defaultdict(float)
    label_w = defaultdict(lambda: defaultdict(float))
    for m, label in calls.items():
        w = sw.get(m, 0.0)
        if w <= 0 or vote_order(label) != top:
            continue
        voted += w
        if '/' in label:
            key = sf_key(label)
            key_w[key] += w
            label_w[key][label] += w
    if not key_w:
        return top
    ranked = sorted(key_w.values(), reverse=True)
    key, kw = max(key_w.items(), key=lambda kv: kv[1])
    if len(ranked) > 1 and abs(ranked[1] - kw) < 1e-9:
        return top
    if kw < cfg['superfamily'].get('min_frac', 0.5) * voted - 1e-9:
        return top
    labels = label_w[key]
    return max(sorted(labels), key=lambda l: labels[l])


# --------------------------------------------------------------------------- #
# Table + FASTA writers
# --------------------------------------------------------------------------- #
def write_table(path, per_model, models, method, callers):
    """Wide voting matrix: rows = repeat ids, cols = models + ensemble winner(s).

    A blank ('NA') cell means that model abstained on that repeat. Every model
    keeps its column; `callers` maps each decision method to its call function.
    """
    ordered = order_models(models)
    ens_cols = ENSEMBLE_COLUMNS[method]
    header = (['repeat_id'] + [MODEL_DISPLAY.get(m, m) for m in ordered]
              + [name for name, _ in ens_cols])
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for base in sorted(per_model):
            ml = per_model[base]
            row = [base] + [ml.get(m, 'NA') for m in ordered]
            for _name, meth in ens_cols:
                winner = callers[meth](ml)
                row.append(winner if winner else 'Unknown')
            fh.write('\t'.join(row) + '\n')


def iter_fasta(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, seq
                header, seq = line, []
            else:
                seq.append(line)
    if header is not None:
        yield header, seq


def is_unknown(label):
    lab = label.strip().lower()
    return lab in UNKNOWN or lab.startswith('unknown')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--library', required=True, help='curated consensus FASTA')
    ap.add_argument('--out', required=True, help='output relabelled FASTA')
    ap.add_argument('--pred', nargs='*', default=[],
                    help='one or more classifier vote TSVs')
    ap.add_argument('--method', default='two_step',
                    choices=['two_step', 'hard', 'weighted', 'all'],
                    help='ensemble strategy (default two_step)')
    ap.add_argument('--kingdom', default='animal',
                    help='organism_type whose two-step settings apply '
                         '(plant / animal / fungus; others use animal)')
    ap.add_argument('--two-step-config', default=None,
                    help='JSON overriding TWO_STEP_DEFAULTS, e.g. {"fungus": {"order": '
                         '{"weights": {"neuralte": 1, "terrier": 2}}}}')
    ap.add_argument('--relabel-orders', default='DNA,RC',
                    help="Orders the ensemble may assign (comma list, matched against "
                         "the call's Order, e.g. 'DNA,RC'); other calls leave the "
                         "consensus Unknown. 'all' = no restriction (default DNA,RC)")
    ap.add_argument('--min-frac', type=float, default=0.5,
                    help='abstention gate: the winner must carry at least this '
                         'fraction of what the voting subset could maximally '
                         'cast (default 0.5; 0 = no gate, plain argmax)')
    ap.add_argument('--min-votes', type=int, default=0,
                    help='hard voting: optional ABSOLUTE floor on agreeing '
                         'votes, applied on top of --min-frac (default 0 = off)')
    ap.add_argument('--min-weight', type=float, default=0.0,
                    help='weighted voting: optional ABSOLUTE floor on summed '
                         'weight, applied on top of --min-frac (default 0 = off)')
    ap.add_argument('--weights', default=None,
                    help="explicit weights, e.g. 'neuralte=1.5,terrier=1.5,"
                         "create=1.3,teclass2=1.2,deepte=1.1' (overrides defaults)")
    ap.add_argument('--weights-from', default=None,
                    help='derive weights from an order-level metrics TSV '
                         '(evaluate_classification all_models.summary_metrics.tsv)')
    ap.add_argument('--weight-metric', default='macro_f1',
                    help='metric column used with --weights-from (default macro_f1)')
    ap.add_argument('--include-models', default=None,
                    help='comma/space list of models whose votes form the ensemble '
                         '(kingdom-aware subset, e.g. "deepte,terrier"). Other '
                         'models still appear in the table but do not vote. Empty '
                         '= all models.')
    ap.add_argument('--table', default=None,
                    help='write a per-model voting matrix TSV to this path')
    ap.add_argument('--min-agree', type=int, default=0,
                    help='a call needs at least this many models of the voting subset '
                         'on the winning Order; never relaxed for small subsets, so a '
                         'single-model subset relabels nothing (default 0 = off)')
    args = ap.parse_args()

    per_model, models = load_per_model(args.pred)
    weights = resolve_weights(args)

    # Kingdom-aware subset: only these models decide the ensemble call. The gate
    # is resolved against whatever this subset can cast (resolve_gates), so a
    # small subset gets a proportionally smaller bar rather than no bar at all.
    include = None
    if args.include_models:
        raw = args.include_models.replace(';', ',').replace(' ', ',')
        include = {x.strip().lower() for x in raw.split(',') if x.strip()}
        missing = include - {m.lower() for m in models}
        sys.stderr.write("ensemble_vote: kingdom-aware subset = %s%s\n" % (
            ",".join(sorted(include)),
            "; not present: %s" % ",".join(sorted(missing)) if missing else ""))

    methods = {meth for _, meth in ENSEMBLE_COLUMNS[args.method]}
    callers = {}
    min_votes, min_weight = max(1, args.min_votes), args.min_weight
    if methods & {'hard', 'weighted'}:
        if 'weighted' in methods:
            sys.stderr.write("ensemble_vote: weights = %s\n" % ', '.join(
                '%s=%.2f' % (m, weights.get(m, 1.0)) for m in order_models(models)))
        min_votes, min_weight = resolve_gates(args, models, weights, include)
        for meth in methods & {'hard', 'weighted'}:
            callers[meth] = (lambda ml, meth=meth: call_with_agreement(
                meth, ml, min_votes, weights, min_weight, include, args.min_agree))
    if 'two_step' in methods:
        ts = load_two_step(args.kingdom, args.two_step_config, include)
        ran = [m.lower() for m in models]
        fmt = lambda w: ','.join('%s=%g' % (m, v) for m, v in w.items()) or 'none'
        sys.stderr.write(
            "ensemble_vote: two-step (%s): order %s, min_frac=%.2f, min_agree=%d; "
            "superfamily %s, min_frac=%.2f\n" % (
                args.kingdom, fmt(ts['order']['weights']), ts['order']['min_frac'],
                ts['order'].get('min_agree', 1), fmt(ts['superfamily']['weights']),
                ts['superfamily'].get('min_frac', 0.5)))
        missing = [m for m in ts['order']['weights'] if m not in ran]
        if missing:
            sys.stderr.write("ensemble_vote: two-step models without votes: %s\n"
                             % ','.join(missing))
        callers['two_step'] = lambda ml: decide_two_step(ml, ts, ran)

    # Order allowlist: a call outside it leaves the consensus Unknown (table and FASTA).
    allowed = None
    if args.relabel_orders and args.relabel_orders.strip().lower() != 'all':
        allowed = {o.strip().upper() for o in args.relabel_orders.replace(';', ',').split(',')
                   if o.strip()}
        blocked = Counter()

        def restrict(call):
            def wrapped(ml):
                w = call(ml)
                if w and vote_order(w).upper() not in allowed:
                    blocked[vote_order(w)] += 1
                    return None
                return w
            return wrapped
        callers = {m: restrict(c) for m, c in callers.items()}
        sys.stderr.write("ensemble_vote: relabel only to Orders %s\n" % ','.join(sorted(allowed)))

    # Voting-results table: rows = repeats, cols = models, cells = each call.
    if args.table:
        write_table(args.table, per_model, models, args.method, callers)
    if allowed is not None:
        blocked.clear()          # count the library's consensi only

    primary = callers[PRIMARY_METHOD[args.method]]
    relabelled = kept = 0
    with open(args.out, 'w') as out:
        for header, seq in iter_fasta(args.library):
            token = header[1:].split()[0] if len(header) > 1 else ''
            desc = header[1:].split(None, 1)
            desc = desc[1] if len(desc) > 1 else ''
            base = token.split('#', 1)[0]
            existing = token.split('#', 1)[1] if '#' in token else ''

            winner = primary(per_model.get(base, {})) if is_unknown(existing) else None

            if winner:
                new_header = '>{}#{}'.format(base, winner)
                if desc:
                    new_header += ' ' + desc
                out.write(new_header + '\n')
                relabelled += 1
            else:
                out.write(header + '\n')
                kept += 1
            for s in seq:
                out.write(s + '\n')

    if allowed is not None:
        sys.stderr.write('ensemble_vote: Unknown consensi kept Unknown by the Order '
                         'allowlist: %s\n'
                         % (', '.join('%s=%d' % kv for kv in sorted(blocked.items())) or 'none'))
    if PRIMARY_METHOD[args.method] == 'two_step':
        sys.stderr.write('ensemble_vote: method=%s (two-step, %s), %d reclassified, '
                         '%d kept\n' % (args.method, args.kingdom, relabelled, kept))
    else:
        sys.stderr.write(
            'ensemble_vote: method=%s, %d reclassified, %d kept '
            '(gate: min_frac=%.2f -> min_votes=%d, min_weight=%.2f; min_agree=%d)\n' % (
                args.method, relabelled, kept, args.min_frac, min_votes, min_weight,
                args.min_agree))


if __name__ == '__main__':
    main()
