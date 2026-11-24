#!/usr/bin/env python3
"""
dfam_partitioner.py — partition a local Dfam release into Plant / Animal / Fungi
ground-truth sets for TE-classification benchmarking.

Two input modes:

  1. EMBL  (Dfam.embl / Dfam_curatedonly.embl[.gz])  -- default, streaming.
     Extracts per family: accession, name, RepeatMasker Type/SubType (from the
     `CC   Type:` / `CC   SubType:` annotations, falling back to `KW`), the
     taxonomic lineage (from `OC` / `OS` / `CC Clades`), and the consensus
     sequence (from the `SQ` block).

  2. TSV   (--tsv summary.tsv [--fasta consensi.fa])  -- for pre-exported tables.
     A generic delimited summary; sequences (optional) are pulled by id from a
     companion FASTA using a streaming index.

For each requested lineage it writes:
  <prefix>_<lineage>.fasta   consensus sequences (id = family accession)
  <prefix>_<lineage>.tsv     ground-truth table (one row per family)
and a global <prefix>_partition_summary.tsv (family counts per lineage x order).

The record is streamed one family at a time, so peak memory is O(one family),
independent of the (multi-GB) Dfam release size.

Example
-------
  python dfam_partitioner.py --embl Dfam_curatedonly.embl.gz \
      --outdir dfam_partitions --prefix dfam --min-length 100
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

import te_label_utils as tu

LOG = tu.setup_logger("dfam_partitioner")

TSV_COLUMNS = [
    "family_id", "family_name", "dfam_accession", "rm_type", "rm_subtype",
    "label", "label_class", "label_order", "label_superfamily", "lineage",
    "seq_length",
]


# --------------------------------------------------------------------------- #
# EMBL streaming parser
# --------------------------------------------------------------------------- #
def iter_embl_records(path: str) -> Iterator[Dict[str, object]]:
    """Yield one dict per EMBL family record. Streams line-by-line.

    Keys: id, name, accession, type, subtype, keywords, taxonomy, sequence.
    """
    with tu.open_text(path) as fh:
        rec: Dict[str, object] = {}
        oc_parts: List[str] = []
        cc_parts: List[str] = []
        seq_chunks: List[str] = []
        in_seq = False

        def flush() -> Optional[Dict[str, object]]:
            if not rec and not seq_chunks:
                return None
            rec["taxonomy"] = " ".join(oc_parts + cc_parts).strip()
            rec["sequence"] = "".join(seq_chunks)
            return rec

        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("//"):
                out = flush()
                if out is not None:
                    yield out
                rec, oc_parts, cc_parts, seq_chunks, in_seq = {}, [], [], [], False
                continue

            if in_seq:
                # Sequence lines: strip whitespace and trailing coordinate digits.
                seq_chunks.append("".join(c for c in line if c.isalpha()))
                continue

            tag = line[:2]
            body = line[5:].strip() if len(line) > 5 else ""

            if tag == "ID":
                rec["id"] = body.split(";", 1)[0].strip() or body.strip()
            elif tag == "AC":
                acc = body.split(";", 1)[0].strip()
                if acc:
                    rec["accession"] = acc
            elif tag == "NM":
                rec["name"] = body.strip()
            elif tag == "DE":
                rec.setdefault("description", body.strip())
            elif tag == "KW":
                rec.setdefault("keywords", "")
                rec["keywords"] = (str(rec["keywords"]) + " " + body).strip()
            elif tag == "OS":
                oc_parts.append(body)
            elif tag == "OC":
                oc_parts.append(body)
            elif tag == "CC":
                cc_parts.append(body)
                low = body.lower()
                if low.startswith("type:"):
                    rec["type"] = body.split(":", 1)[1].strip()
                elif low.startswith("subtype:"):
                    rec["subtype"] = body.split(":", 1)[1].strip()
            elif tag == "SQ":
                in_seq = True

        out = flush()
        if out is not None:
            yield out


def resolve_class(rec: Dict[str, object]) -> Tuple[str, str, str]:
    """Return (rm_type, rm_subtype, raw_label) for an EMBL record."""
    rm_type = str(rec.get("type", "") or "").strip()
    rm_subtype = str(rec.get("subtype", "") or "").strip()
    if not rm_type:
        # Fall back to the first '/'-bearing keyword, e.g. "SINE/MIR."
        for kw in str(rec.get("keywords", "")).replace(".", " ").split():
            if "/" in kw:
                parts = kw.split("/", 1)
                rm_type, rm_subtype = parts[0], parts[1]
                break
        else:
            kw0 = str(rec.get("keywords", "")).replace(".", " ").split()
            rm_type = kw0[0] if kw0 else "Unknown"
    raw = f"{rm_type}/{rm_subtype}" if rm_subtype else rm_type
    return rm_type, rm_subtype, raw


# --------------------------------------------------------------------------- #
# TSV input mode
# --------------------------------------------------------------------------- #
def iter_tsv_records(path: str, id_col: str, class_col: str,
                     tax_col: Optional[str], name_col: Optional[str],
                     sep: str) -> Iterator[Dict[str, object]]:
    for chunk in tu.iter_tsv_chunks(path, sep=sep):
        cols = list(chunk.columns)
        cid = tu.detect_column(cols, [id_col, "family_id", "accession", "id", "name"])
        ccl = tu.detect_column(cols, [class_col, "classification", "class", "type", "label"])
        ctx = tu.detect_column(cols, [tax_col] if tax_col else []) or \
            tu.detect_column(cols, ["taxonomy", "clade", "lineage", "organism"])
        cnm = tu.detect_column(cols, [name_col] if name_col else []) or \
            tu.detect_column(cols, ["family_name", "name"])
        if cid is None or ccl is None:
            raise SystemExit(
                f"ERROR: could not find id/class columns in {path}. "
                f"Available: {cols}")
        for _, row in chunk.iterrows():
            raw = str(row[ccl])
            rm_type = raw.split("/", 1)[0] if "/" in raw else raw
            rm_subtype = raw.split("/", 1)[1] if "/" in raw else ""
            yield {
                "id": str(row[cid]),
                "name": str(row[cnm]) if cnm else "",
                "accession": str(row[cid]),
                "type": rm_type, "subtype": rm_subtype,
                "keywords": raw,
                "taxonomy": str(row[ctx]) if ctx else "",
                "sequence": "",
            }


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
class LineageWriters:
    """Lazily-opened per-lineage FASTA + TSV writers."""

    def __init__(self, outdir: str, prefix: str, write_fasta: bool):
        self.outdir = tu.ensure_dir(outdir)
        self.prefix = prefix
        self.write_fasta = write_fasta
        self._fa: Dict[str, TextIO] = {}
        self._tsv: Dict[str, TextIO] = {}
        self._tsv_w: Dict[str, csv.writer] = {}

    def _fasta(self, lineage: str) -> TextIO:
        if lineage not in self._fa:
            p = os.path.join(self.outdir, f"{self.prefix}_{lineage}.fasta")
            self._fa[lineage] = open(p, "w", encoding="utf-8")
        return self._fa[lineage]

    def _table(self, lineage: str) -> csv.writer:
        if lineage not in self._tsv_w:
            p = os.path.join(self.outdir, f"{self.prefix}_{lineage}.tsv")
            fh = open(p, "w", encoding="utf-8", newline="")
            w = csv.writer(fh, delimiter="\t")
            w.writerow(TSV_COLUMNS)
            self._tsv[lineage] = fh
            self._tsv_w[lineage] = w
        return self._tsv_w[lineage]

    def write(self, lineage: str, row: Dict[str, object], sequence: str) -> None:
        self._table(lineage).writerow([row[c] for c in TSV_COLUMNS])
        if self.write_fasta and sequence:
            fh = self._fasta(lineage)
            desc = f"{row['family_name']}|{row['label']}|{lineage}"
            fh.write(f">{row['family_id']} {desc}\n")
            for i in range(0, len(sequence), 60):
                fh.write(sequence[i:i + 60] + "\n")

    def close(self) -> None:
        for fh in list(self._fa.values()) + list(self._tsv.values()):
            fh.close()


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Partition a local Dfam release into Plant/Animal/Fungi "
                    "ground-truth sets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--embl", help="Dfam EMBL file (.embl or .embl.gz).")
    src.add_argument("--tsv", help="Pre-exported Dfam summary TSV.")
    p.add_argument("--fasta", help="Companion consensus FASTA (TSV mode; "
                                   "sequences fetched by id).")
    p.add_argument("--outdir", default="dfam_partitions", help="Output directory.")
    p.add_argument("--prefix", default="dfam", help="Output file-name prefix.")
    p.add_argument("--lineages", nargs="+", default=["plant", "animal", "fungi"],
                   choices=["plant", "animal", "fungi"],
                   help="Lineages to emit.")
    p.add_argument("--min-length", type=int, default=0,
                   help="Drop families whose consensus is shorter than this (bp).")
    p.add_argument("--keep-nonte", action="store_true",
                   help="Keep non-TE repeats (satellites, rRNA, simple) that "
                        "normalise to order 'Other' (dropped by default).")
    p.add_argument("--keep-other-lineage", action="store_true",
                   help="Also emit families not in plant/animal/fungi into an "
                        "'other' partition.")
    p.add_argument("--label-map", help="JSON of extra label synonyms (see "
                                        "te_label_utils.load_label_overrides).")
    p.add_argument("--no-fasta", action="store_true",
                   help="Write only the TSV tables, not per-lineage FASTA.")
    # TSV-mode column hints
    p.add_argument("--id-col", default="family_id")
    p.add_argument("--class-col", default="classification")
    p.add_argument("--tax-col", default=None)
    p.add_argument("--name-col", default=None)
    p.add_argument("--sep", default="\t")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)

    wanted = set(args.lineages)
    if args.keep_other_lineage:
        wanted.add("other")

    # Optional FASTA index for TSV mode (streaming random access).
    seq_index = None
    if args.tsv and args.fasta and not args.no_fasta:
        from Bio import SeqIO  # local import
        LOG.info("Indexing companion FASTA %s ...", args.fasta)
        seq_index = SeqIO.index(args.fasta, "fasta")

    if args.embl:
        record_iter = iter_embl_records(args.embl)
    else:
        record_iter = iter_tsv_records(
            args.tsv, args.id_col, args.class_col, args.tax_col,
            args.name_col, args.sep)

    writers = LineageWriters(args.outdir, args.prefix,
                             write_fasta=not args.no_fasta)
    counts: "Counter[Tuple[str, str]]" = Counter()
    stats = Counter()
    try:
        for rec in record_iter:
            stats["total"] += 1
            fam_id = str(rec.get("id") or rec.get("accession") or "").strip()
            if not fam_id:
                stats["no_id"] += 1
                continue

            lineage = tu.classify_lineage(str(rec.get("taxonomy", "")))
            if lineage is None:
                lineage = "other"
            if lineage not in wanted:
                stats["lineage_skipped"] += 1
                continue

            rm_type, rm_subtype, raw = resolve_class(rec)
            te_class = tu.normalize_label(raw, "class")
            order = tu.normalize_label(raw, "order")
            superfam = tu.normalize_label(raw, "superfamily")

            if order == "Other" and not args.keep_nonte:
                stats["nonte_dropped"] += 1
                continue

            sequence = str(rec.get("sequence", "") or "")
            if not sequence and seq_index is not None and fam_id in seq_index:
                sequence = str(seq_index[fam_id].seq)
            if args.min_length and len(sequence) < args.min_length:
                stats["too_short"] += 1
                continue

            row = {
                "family_id": fam_id,
                "family_name": str(rec.get("name") or rec.get("description") or ""),
                "dfam_accession": str(rec.get("accession") or fam_id),
                "rm_type": rm_type,
                "rm_subtype": rm_subtype,
                "label": raw,
                "label_class": te_class,
                "label_order": order,
                "label_superfamily": superfam,
                "lineage": lineage,
                "seq_length": len(sequence),
            }
            writers.write(lineage, row, sequence)
            counts[(lineage, order)] += 1
            stats["written"] += 1
            if stats["total"] % 50000 == 0:
                LOG.info("... parsed %d families (%d written)",
                         stats["total"], stats["written"])
    finally:
        writers.close()
        if seq_index is not None:
            seq_index.close()

    # Partition summary.
    summ = os.path.join(args.outdir, f"{args.prefix}_partition_summary.tsv")
    with open(summ, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["lineage", "order", "n_families"])
        for (lineage, order), n in sorted(counts.items()):
            w.writerow([lineage, order, n])

    LOG.info("Done. total=%d written=%d "
             "(no_id=%d lineage_skipped=%d nonte_dropped=%d too_short=%d)",
             stats["total"], stats["written"], stats["no_id"],
             stats["lineage_skipped"], stats["nonte_dropped"], stats["too_short"])
    LOG.info("Partitions + summary written to %s/", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
