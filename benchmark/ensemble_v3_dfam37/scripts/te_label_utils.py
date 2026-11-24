#!/usr/bin/env python3
"""
Shared utilities for the teannotate classification-benchmark suite.

The four benchmark scripts (dfam_partitioner, evaluate_classification,
calculate_curves, plot_benchmarks) all need to reconcile the *different*
class vocabularies emitted by heterogeneous tools:

  * Dfam / RepeatMasker  ->  "LTR/Gypsy", "DNA/TcMar-Tc1", "SINE/MIR"
  * DeepTE               ->  "ClassI_LTR_Copia", "ClassII_DNA_Mutator", "nMITE"
  * TEclass2 / CREATE    ->  Wicker 3-letter codes  ("RLG", "DTT", "DHH")
  * RepeatModeler style   ->  "id#LTR/Gypsy"

This module normalises any of those into a common controlled vocabulary at two
granularities:

  * level="order"        ->  LTR, LINE, SINE, DNA, RC, Maverick, DIRS, PLE,
                             Unknown, Other
  * level="superfamily"  ->  Gypsy, Copia, L1, hAT, Tc1-Mariner, Helitron, ...

Everything here is dependency-light (std-lib only) so it imports cheaply on a
compute node before the heavy scientific stack is touched.
"""
from __future__ import annotations

import gzip
import io
import json
import logging
import os
import re
import sys
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

# --------------------------------------------------------------------------- #
# Controlled vocabulary
# --------------------------------------------------------------------------- #
CANONICAL_ORDERS: List[str] = [
    "LTR", "LINE", "SINE", "DNA", "RC", "Maverick", "DIRS", "PLE",
    "Unknown", "Other",
]

# TE class = Wicker top level: Class I (retrotransposons, RNA intermediate) vs
# Class II (DNA transposons). Used for --level class.
CANONICAL_CLASSES: List[str] = ["ClassI", "ClassII"]

ORDER_TO_CLASS: Dict[str, str] = {
    "LTR": "ClassI", "LINE": "ClassI", "SINE": "ClassI",
    "DIRS": "ClassI", "PLE": "ClassI",
    "DNA": "ClassII", "RC": "ClassII", "Maverick": "ClassII",
    "Unknown": "Unknown", "Other": "Other",
}

# Direct class tokens, for labels that name the class but resolve to no order
# (e.g. a bare "ClassI"/"CLASSI"/"Retrotransposon"). Order-derived class wins;
# these are only a fallback.
_CLASS_TOKEN_MAP: Dict[str, str] = {
    "classi": "ClassI", "class1": "ClassI", "classone": "ClassI",
    "retrotransposon": "ClassI", "retrotransposons": "ClassI",
    "retroelement": "ClassI",
    "classii": "ClassII", "class2": "ClassII", "classtwo": "ClassII",
    "dnatransposon": "ClassII", "dnatransposons": "ClassII",
}

# superfamily token -> (canonical superfamily, canonical order)
SUPERFAMILY_MAP: Dict[str, Tuple[str, str]] = {
    # ---- LTR retrotransposons -------------------------------------------- #
    "gypsy": ("Gypsy", "LTR"), "ty3": ("Gypsy", "LTR"), "ty3_gypsy": ("Gypsy", "LTR"),
    "metavir": ("Gypsy", "LTR"),
    "copia": ("Copia", "LTR"), "ty1": ("Copia", "LTR"), "ty1_copia": ("Copia", "LTR"),
    "pseudovir": ("Copia", "LTR"),
    "belpao": ("Bel-Pao", "LTR"), "bel": ("Bel-Pao", "LTR"), "pao": ("Bel-Pao", "LTR"),
    "erv": ("ERV", "LTR"), "ervk": ("ERV", "LTR"), "ervl": ("ERV", "LTR"),
    "erv1": ("ERV", "LTR"), "retrovirus": ("ERV", "LTR"), "lentivirus": ("ERV", "LTR"),
    "caulimovirus": ("Caulimovirus", "LTR"), "trim": ("TRIM", "LTR"),
    "lard": ("LARD", "LTR"),
    # ---- DIRS ------------------------------------------------------------- #
    "dirs": ("DIRS", "DIRS"), "ngaro": ("Ngaro", "DIRS"), "viper": ("VIPER", "DIRS"),
    # ---- Penelope --------------------------------------------------------- #
    "penelope": ("Penelope", "PLE"), "ple": ("Penelope", "PLE"),
    # ---- LINEs ------------------------------------------------------------ #
    "l1": ("L1", "LINE"), "l1md": ("L1", "LINE"), "tx1": ("L1", "LINE"),
    "l2": ("L2", "LINE"), "cr1": ("CR1", "LINE"), "rex": ("CR1", "LINE"),
    "rte": ("RTE", "LINE"), "bovb": ("RTE", "LINE"),
    "jockey": ("Jockey", "LINE"), "cr1_zenon": ("CR1", "LINE"),
    "i": ("I", "LINE"), "r1": ("R1", "LINE"), "r2": ("R2", "LINE"),
    "loa": ("I", "LINE"), "nlr": ("LINE", "LINE"), "nlr_clade": ("LINE", "LINE"),
    "proto2": ("RTE", "LINE"),
    # ---- SINEs ------------------------------------------------------------ #
    "alu": ("Alu", "SINE"), "b1": ("SINE", "SINE"), "b2": ("SINE", "SINE"),
    "mir": ("MIR", "SINE"), "trna": ("tRNA", "SINE"), "7sl": ("7SL", "SINE"),
    "5s": ("5S", "SINE"), "id": ("SINE", "SINE"), "b4": ("SINE", "SINE"),
    "sine2": ("tRNA", "SINE"), "sine3": ("5S", "SINE"),
    # ---- DNA / TIR transposons ------------------------------------------- #
    "hat": ("hAT", "DNA"), "hat_ac": ("hAT", "DNA"), "hat_charlie": ("hAT", "DNA"),
    "hat_tag1": ("hAT", "DNA"), "hat_tip100": ("hAT", "DNA"), "charlie": ("hAT", "DNA"),
    "tc1": ("Tc1-Mariner", "DNA"), "mariner": ("Tc1-Mariner", "DNA"),
    "tcmar": ("Tc1-Mariner", "DNA"), "tc1_mariner": ("Tc1-Mariner", "DNA"),
    "stowaway": ("Tc1-Mariner", "DNA"), "pogo": ("Tc1-Mariner", "DNA"),
    "mutator": ("Mutator", "DNA"), "mule": ("Mutator", "DNA"), "mudr": ("Mutator", "DNA"),
    "cacta": ("CACTA", "DNA"), "enspm": ("CACTA", "DNA"), "cmc": ("CACTA", "DNA"),
    "harbinger": ("PIF-Harbinger", "DNA"), "pif": ("PIF-Harbinger", "DNA"),
    "tourist": ("PIF-Harbinger", "DNA"), "pif_harbinger": ("PIF-Harbinger", "DNA"),
    "p": ("P", "DNA"), "piggybac": ("PiggyBac", "DNA"), "pgbd": ("PiggyBac", "DNA"),
    "merlin": ("Merlin", "DNA"), "transib": ("Transib", "DNA"),
    "kolobok": ("Kolobok", "DNA"), "academ": ("Academ", "DNA"),
    "ginger": ("Ginger", "DNA"), "sola": ("Sola", "DNA"), "zator": ("Zator", "DNA"),
    "novosib": ("Novosib", "DNA"), "dada": ("Dada", "DNA"),
    "crypton": ("Crypton", "DNA"), "mite": ("MITE", "DNA"),
    # ---- Rolling-circle (Helitron) --------------------------------------- #
    "helitron": ("Helitron", "RC"), "heli": ("Helitron", "RC"),
    # ---- Maverick / Polinton --------------------------------------------- #
    "maverick": ("Maverick", "Maverick"), "polinton": ("Maverick", "Maverick"),
}

# order token -> canonical order
ORDER_MAP: Dict[str, str] = {
    "ltr": "LTR", "line": "LINE", "sine": "SINE",
    "dna": "DNA", "tir": "DNA", "dnatransposon": "DNA",
    "rc": "RC", "helitron": "RC", "rollingcircle": "RC", "rolling_circle": "RC",
    "maverick": "Maverick", "polinton": "Maverick",
    "dirs": "DIRS", "ple": "PLE", "penelope": "PLE",
    "retrotransposon": "LTR", "retroelement": "LINE",
    "nonltr": "LINE", "nltr": "LINE",
    # NB: 'classii' is a CLASS, not an order — Class II spans DNA(TIR), RC(Helitron)
    # and Maverick. It must NOT set the order, or a nested label like PanTEon's
    # "CLASSII/HELITRON/HELITRON" would resolve to DNA and leave the RC row empty.
    # The specific superfamily (helitron->RC, maverick->Maverick, tir->DNA) sets the
    # order instead; class-level resolution uses _CLASS_TOKEN_MAP (see normalize_label).
    "unknown": "Unknown", "unclassified": "Unknown", "unspecified": "Unknown",
    "na": "Unknown", "none": "Unknown", "nan": "Unknown", "": "Unknown",
    # non-TE repeats: mapped to Other so they never inflate TE scores
    "satellite": "Other", "simple_repeat": "Other", "simple": "Other",
    "low_complexity": "Other", "rrna": "Other", "trna_gene": "Other",
    "snrna": "Other", "scrna": "Other", "artefact": "Other", "other": "Other",
    "segmental": "Other", "tandem": "Other",
}

# Wicker (2007) 3-letter codes -> (superfamily, order)
WICKER_MAP: Dict[str, Tuple[str, str]] = {
    "rlc": ("Copia", "LTR"), "rlg": ("Gypsy", "LTR"), "rlb": ("Bel-Pao", "LTR"),
    "rlr": ("ERV", "LTR"), "rle": ("ERV", "LTR"), "rlx": ("LTR", "LTR"),
    "ryd": ("DIRS", "DIRS"), "ryn": ("Ngaro", "DIRS"), "ryv": ("VIPER", "DIRS"),
    "ryx": ("DIRS", "DIRS"),
    "rpp": ("Penelope", "PLE"), "rpx": ("PLE", "PLE"),
    "rir": ("R2", "LINE"), "rit": ("RTE", "LINE"), "rij": ("Jockey", "LINE"),
    "ril": ("L1", "LINE"), "rii": ("I", "LINE"), "rix": ("LINE", "LINE"),
    "rst": ("tRNA", "SINE"), "rsl": ("7SL", "SINE"), "rss": ("5S", "SINE"),
    "rsx": ("SINE", "SINE"),
    "dta": ("hAT", "DNA"), "dtc": ("CACTA", "DNA"), "dte": ("Merlin", "DNA"),
    "dth": ("PIF-Harbinger", "DNA"), "dtm": ("Mutator", "DNA"), "dtp": ("P", "DNA"),
    "dtt": ("Tc1-Mariner", "DNA"), "dtb": ("PiggyBac", "DNA"), "dtr": ("Transib", "DNA"),
    "dtx": ("DNA", "DNA"), "dyc": ("Crypton", "DNA"),
    "dhh": ("Helitron", "RC"), "dmm": ("Maverick", "Maverick"),
    "xxx": ("Unknown", "Unknown"),
}

_UNKNOWN_TOKENS = {
    "", "na", "n/a", "none", "nan", "null", "unknown", "unclassified",
    "unspecified", ".", "-", "?", "abstain", "no_call",
}

_SPLIT_RE = re.compile(r"[\/\|_\-;:,\.\s]+")
_WICKER_RE = re.compile(r"^[rd][a-z][a-z]$")


# --------------------------------------------------------------------------- #
# Label normalisation
# --------------------------------------------------------------------------- #
def parse_raw_label(raw: object) -> Tuple[str, str]:
    """Return (order, superfamily) canonical strings for an arbitrary raw label.

    Robust to FASTA headers ('id#LTR/Gypsy'), DeepTE ('ClassI_LTR_Copia'),
    Wicker codes ('RLG'), plain superfamilies ('Gypsy') and RepeatMasker
    'Type/SubType'. Missing / abstain values collapse to ('Unknown','Unknown').
    """
    if raw is None:
        return ("Unknown", "Unknown")
    text = str(raw).strip()
    if not text:
        return ("Unknown", "Unknown")

    # RepeatModeler-style "family_id#Class/Subfamily": keep the part after '#'.
    if "#" in text:
        text = text.split("#", 1)[1]
    # Drop a leading '>' and take the first whitespace-delimited token only.
    text = text.lstrip(">").strip()
    if not text:
        return ("Unknown", "Unknown")
    first = text.split()[0]
    low = first.lower()

    if low in _UNKNOWN_TOKENS:
        return ("Unknown", "Unknown")

    # Whole-token Wicker code (e.g. "RLG").
    if _WICKER_RE.match(low) and low in WICKER_MAP:
        return (WICKER_MAP[low][1], WICKER_MAP[low][0])

    parts = [p for p in _SPLIT_RE.split(low) if p]

    order: Optional[str] = None
    superfam: Optional[str] = None

    for p in parts:
        if _WICKER_RE.match(p) and p in WICKER_MAP:
            sf, od = WICKER_MAP[p]
            superfam = superfam or sf
            order = order or od
        if superfam is None and p in SUPERFAMILY_MAP:
            sf, od = SUPERFAMILY_MAP[p]
            superfam = sf
            order = order or od
        if order is None and p in ORDER_MAP:
            order = ORDER_MAP[p]

    # Infer a missing order from a resolved superfamily, and vice-versa.
    if order is None and superfam is not None:
        for sf, od in SUPERFAMILY_MAP.values():
            if sf == superfam:
                order = od
                break
    if order is None:
        order = "Other"
    if superfam is None:
        superfam = order if order in ("Unknown", "Other") else order
    return (order, superfam)


def _class_from_tokens(raw: object) -> Optional[str]:
    """Detect an explicit TE class token in a raw label, else None."""
    if raw is None:
        return None
    text = str(raw).strip()
    if "#" in text:
        text = text.split("#", 1)[1]
    text = text.lstrip(">").strip().lower()
    for p in _SPLIT_RE.split(text):
        if p in _CLASS_TOKEN_MAP:
            return _CLASS_TOKEN_MAP[p]
    return None


def normalize_label(raw: object, level: str = "order") -> str:
    """Normalise a raw label to the requested granularity.

    level='class'       -> ClassI, ClassII, Unknown, Other
    level='order'       -> LTR, LINE, SINE, DNA, RC, Maverick, DIRS, PLE, ...
    level='superfamily' -> Gypsy, Copia, L1, hAT, Tc1-Mariner, Helitron, ...
    """
    order, superfam = parse_raw_label(raw)
    if level == "superfamily":
        return superfam
    if level == "class":
        cls = ORDER_TO_CLASS.get(order, "Other")
        # If the order didn't pin a class, fall back to an explicit class token.
        if cls in ("Unknown", "Other"):
            tok = _class_from_tokens(raw)
            if tok:
                return tok
        return cls
    return order


def parse_panteon_header(header: str) -> Tuple[str, str, str]:
    """Split a PanTEon FASTA header into (family_id, raw_label, species).

    PanTEon headers follow  ``>ID#Class/Order/Superfamily @Species binomial``
    (Orozco-Arias et al. 2026), e.g.::

        >PDB00018676#CLASSI/LINE/CR1 @Bucorvus abyssinicus

    Returns:
      family_id  -- the id before '#'                      ("PDB00018676")
      raw_label  -- the first whitespace token, '#'-tagged  ("PDB00018676#CLASSI/LINE/CR1")
                    so normalize_label()/parse_raw_label() resolve order+superfamily.
      species    -- text after '@', stripped               ("Bucorvus abyssinicus")

    Missing '@' -> species is "". Missing '#' -> raw_label is the bare token.
    """
    h = (header or "").lstrip(">").strip()
    if not h:
        return ("", "", "")
    species = ""
    if "@" in h:
        left, species = h.split("@", 1)
        species = species.strip()
    else:
        left = h
    toks = left.split()
    token = toks[0] if toks else left.strip()
    family_id = token.split("#", 1)[0]
    return (family_id, token, species)


def load_label_overrides(path: Optional[str]) -> None:
    """Merge a user JSON of extra synonyms into the global maps in-place.

    JSON schema (all keys optional):
      {"order": {"token": "LTR"},
       "superfamily": {"token": ["Gypsy", "LTR"]},
       "wicker": {"rlz": ["Gypsy", "LTR"]}}
    """
    if not path:
        return
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    for tok, val in (data.get("order") or {}).items():
        ORDER_MAP[tok.lower()] = str(val)
    for tok, val in (data.get("superfamily") or {}).items():
        SUPERFAMILY_MAP[tok.lower()] = (str(val[0]), str(val[1]))
    for tok, val in (data.get("wicker") or {}).items():
        WICKER_MAP[tok.lower()] = (str(val[0]), str(val[1]))


# --------------------------------------------------------------------------- #
# I/O helpers (gzip-aware, streaming)
# --------------------------------------------------------------------------- #
def open_text(path: str, mode: str = "rt") -> io.TextIOBase:
    """Open a plain or .gz text file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def detect_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    """Case-insensitive best-match of a column name from a candidate list."""
    lut = {str(c).lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lut:
            return lut[cand.lower()]
    return None


def iter_tsv_chunks(path: str, chunksize: int = 200_000,
                    usecols=None, sep: str = "\t", dtype=str):
    """Yield DataFrame chunks from a (optionally gzipped) delimited file.

    Kept in this module so every script streams large tables identically. Imports
    pandas lazily so std-lib-only consumers (dfam_partitioner) stay light.
    """
    import pandas as pd  # local import: heavy dependency
    compression = "gzip" if str(path).endswith(".gz") else None
    reader = pd.read_csv(
        path, sep=sep, usecols=usecols, dtype=dtype, chunksize=chunksize,
        compression=compression, keep_default_na=False, na_values=[],
    )
    for chunk in reader:
        yield chunk


def setup_logger(name: str, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s", "%H:%M:%S"))
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


# Default lineage keyword sets (lower-cased substrings tested against taxonomy).
LINEAGE_KEYWORDS: Dict[str, Tuple[str, ...]] = {
    "plant": ("viridiplantae", "streptophyta", "chlorophyta", "embryophyta"),
    "animal": ("metazoa", "eumetazoa", "bilateria"),
    "fungi": ("fungi", "dikarya", "ascomycota", "basidiomycota", "mucoromycota"),
}


def classify_lineage(taxonomy_blob: str,
                     keywords: Dict[str, Tuple[str, ...]] = LINEAGE_KEYWORDS
                     ) -> Optional[str]:
    """Return 'plant'|'animal'|'fungi' for a taxonomy string, else None."""
    low = (taxonomy_blob or "").lower()
    for lineage, keys in keywords.items():
        if any(k in low for k in keys):
            return lineage
    return None
