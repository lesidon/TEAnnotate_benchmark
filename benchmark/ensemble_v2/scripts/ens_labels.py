"""ens_labels.py - one label resolution for truth and predictions in ensemble_v2.

te_label_utils resolves the order from the first order token, so RepeatMasker-style
'LTR/DIRS' becomes order LTR while PanTEon's 'CLASSI/DIRS/DIRS' becomes DIRS, and
'LINE/Penelope' becomes LINE while 'PLE' becomes PLE. A correct call would then be
scored wrong depending on the source's convention. Here the superfamily decides the
order for those lineages, identically for every source and every classifier.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from te_label_utils import load_label_overrides, parse_raw_label  # noqa: E402,F401

SF_ORDER = {"DIRS": "DIRS", "Ngaro": "DIRS", "VIPER": "DIRS", "Penelope": "PLE"}
NON_TE = {"Unknown", "Other"}

# held-out genomes of the whole-pipeline benchmark (one-word entry = whole genus)
HELDOUT_SPECIES = [
    "Caenorhabditis elegans", "Caenorhabditis briggsae", "Blumeria",
    "Chlamydomonas reinhardtii", "Ectocarpus", "Gallus gallus", "Taeniopygia guttata",
    "Mus musculus",
]


def order_sf(raw):
    """(order, superfamily) with the DIRS / Penelope order fixed by the superfamily."""
    o, sf = parse_raw_label(raw)
    if sf in SF_ORDER:
        o = SF_ORDER[sf]
    return o, sf


def to_level(raw, level):
    o, sf = order_sf(raw)
    return sf if level == "superfamily" else o


def is_heldout(species):
    s = " ".join((species or "").strip().lower().replace("_", " ").split())
    for b in HELDOUT_SPECIES:
        b = b.lower()
        if s == b or s.startswith(b + " "):
            return True
    return False

