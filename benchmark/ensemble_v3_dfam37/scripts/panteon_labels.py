#!/usr/bin/env python3
"""
panteon_labels.py — build a ground-truth label table from a PanTEon FASTA.

Run this on the *leakage-filtered* held-out FASTA produced by leakage_filter.py
(``<kingdom>.heldout.fasta``) to get the matching truth table for
evaluate_classification.py. Because the labels come from exactly the sequences
you will classify, prediction input and ground truth stay perfectly aligned —
no id bookkeeping or --keep-list needed.

Parses each PanTEon header ``>ID#Class/Order/Superfamily @Species`` and writes:
  family_id, species, label, label_order, label_superfamily, seq_length
(+ a constant 'kingdom' column when --kingdom is given).

The output plugs straight into:
  evaluate_classification.py --truth <out>.tsv --truth-id-col family_id \
      --truth-label-col label --level superfamily ...

Example
-------
  python panteon_labels.py \
      --fasta leakage_out/panteon_fungi.heldout.fasta \
      --kingdom fungi --out truth/panteon_fungi.heldout.truth.tsv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from typing import Iterator, List, Optional, Tuple

import te_label_utils as tu

LOG = tu.setup_logger("panteon_labels")

BASE_COLUMNS = [
    "family_id", "species", "label", "label_class", "label_order",
    "label_superfamily", "seq_length",
]


def iter_fasta_headers(path: str) -> Iterator[Tuple[str, int]]:
    """Yield (header, seq_length) streaming — no need to hold sequences."""
    header: Optional[str] = None
    length = 0
    with tu.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, length
                header = line[1:].rstrip("\n")
                length = 0
            elif header is not None:
                length += len(line.strip())
    if header is not None:
        yield header, length


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Extract a ground-truth label TSV from a PanTEon FASTA "
                    "(e.g. the leakage-filtered held-out set).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--fasta", required=True,
                   help="PanTEon-format FASTA (.fasta/.gz).")
    p.add_argument("--out", required=True, help="Output truth TSV.")
    p.add_argument("--kingdom", default=None,
                   help="Constant kingdom value to add as a column (optional).")
    p.add_argument("--keep-nonte", action="store_true",
                   help="Keep sequences whose label normalises to Unknown/Other "
                        "(dropped by default, to match dfam_partitioner).")
    p.add_argument("--label-map", help="JSON of extra label synonyms.")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)

    outdir = os.path.dirname(os.path.abspath(args.out))
    tu.ensure_dir(outdir)

    columns = list(BASE_COLUMNS)
    if args.kingdom:
        columns.append("kingdom")

    stats: "Counter[str]" = Counter()
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(columns)
        for header, length in iter_fasta_headers(args.fasta):
            stats["total"] += 1
            fam_id, raw_label, species = tu.parse_panteon_header(header)
            if not fam_id:
                stats["no_id"] += 1
                continue
            te_class = tu.normalize_label(raw_label, "class")
            order = tu.normalize_label(raw_label, "order")
            superfam = tu.normalize_label(raw_label, "superfamily")
            if not args.keep_nonte and order in ("Unknown", "Other"):
                stats["nonte_dropped"] += 1
                continue
            row: List[object] = [fam_id, species, raw_label, te_class, order,
                                 superfam, length]
            if args.kingdom:
                row.append(args.kingdom)
            w.writerow(row)
            stats["written"] += 1

    LOG.info("Wrote %d labels (total=%d, no_id=%d, nonte_dropped=%d) -> %s",
             stats["written"], stats["total"], stats["no_id"],
             stats["nonte_dropped"], args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
