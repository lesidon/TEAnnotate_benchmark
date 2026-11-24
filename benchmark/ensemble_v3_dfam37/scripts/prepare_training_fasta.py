#!/usr/bin/env python3
"""
prepare_training_fasta.py — convert Dfam / RepBase releases into the plain,
uncompressed, single-line-per-record nucleotide FASTA that leakage_filter.py
(MMseqs2) needs as a ``--training`` corpus.

The two training databases ship in awkward shapes:

  * Dfam    ``Dfam_curatedonly.embl[.gz]``            -- one big EMBL file.
  * RepBase ``RepBase29.10.fasta.tar.gz``             -- a tar of per-taxon
            ``*.ref`` files (angrep.ref, celrep.ref, diarep.ref, ...). Despite
            the ".fasta" in the tarball name these ``.ref`` files may be EMBL
            *or* FASTA depending on the edition.

This script does not care which: it auto-detects EMBL vs FASTA per member,
streams each record (peak memory = one record), keeps only nucleotide
sequences, and writes one combined FASTA per source. Record ids are prefixed
with the source name (``repbase__MARINER1_DR``) so target ids stay unique and
traceable across the merged corpus.

Inputs may be a file, a directory, a glob, or a ``.tar[.gz]`` archive, and may
be gzip-compressed. Multiple sources are given as ``NAME=PATH`` so the output
FASTA is named after NAME.

Example
-------
  python prepare_training_fasta.py \
      --source dfam=Dfam_curatedonly.embl.gz \
      --source repbase=RepBase29.10.fasta.tar.gz \
      --outdir training_fasta

  # -> training_fasta/dfam.fasta, training_fasta/repbase.fasta
  #    plus the ready-to-run leakage_filter.py --training line (printed).
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import sys
import tarfile
from glob import glob
from typing import Iterator, List, Optional, TextIO, Tuple

import te_label_utils as tu

LOG = tu.setup_logger("prepare_training_fasta")

_FASTA_EXTS = (".ref", ".embl", ".fa", ".fasta", ".fna", ".txt", ".dat", ".seq")
_NUCL = set("ACGTNURYSWKMBDHVacgtnuryswkmbdhv")

# EMBL flat-file line: a 2-letter tag in fixed columns (e.g. "CC   ", "ID   ",
# "SQ   ", "AC   ", "XX") or the "//" record terminator.
_EMBL_LINE = re.compile(r"^([A-Z]{2})( {3,}|\t)")


# --------------------------------------------------------------------------- #
# Streaming record readers (format auto-detected per handle)
# --------------------------------------------------------------------------- #
def _sniff_format(handle: TextIO, max_scan: int = 500) -> str:
    """Scan the head of a handle to decide 'fasta' | 'embl' | 'unknown'.

    Robust to leading comment/preamble blocks: Dfam's EMBL release opens with a
    ``CC`` license header (no ``ID`` line for hundreds of bytes), so we look for
    *any* EMBL tag line — not just ``ID`` — and stop at the first '>' (FASTA) or
    EMBL signature, whichever comes first.
    """
    pos = handle.tell()
    fmt = "unknown"
    scanned = 0
    for line in handle:
        s = line.rstrip("\n")
        if not s.strip():
            continue
        scanned += 1
        if s.startswith(">"):
            fmt = "fasta"
            break
        if s.startswith("//") or s.strip() == "XX" or _EMBL_LINE.match(s):
            fmt = "embl"
            break
        if s.startswith((";", "#")):
            continue  # comment lines seen in some .ref headers
        if scanned >= max_scan:
            break
    handle.seek(pos)
    return fmt


def _iter_fasta(handle: TextIO) -> Iterator[Tuple[str, str]]:
    rec_id: Optional[str] = None
    chunks: List[str] = []
    for line in handle:
        if line.startswith(">"):
            if rec_id is not None:
                yield rec_id, "".join(chunks)
            header = line[1:].strip()
            rec_id = header.split()[0] if header else ""
            chunks = []
        elif rec_id is not None:
            chunks.append(line.strip())
    if rec_id is not None:
        yield rec_id, "".join(chunks)


def _iter_embl(handle: TextIO) -> Iterator[Tuple[str, str]]:
    """Lean EMBL reader: (id, sequence). id = ID/NM/AC first token, in order."""
    rec_id: Optional[str] = None
    acc = ""
    seq_chunks: List[str] = []
    in_seq = False
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith("//"):
            if rec_id or seq_chunks:
                yield (rec_id or acc or "unnamed"), "".join(seq_chunks)
            rec_id, acc, seq_chunks, in_seq = None, "", [], False
            continue
        if in_seq:
            seq_chunks.append("".join(c for c in line if c.isalpha()))
            continue
        tag = line[:2]
        body = line[5:].strip() if len(line) > 5 else ""
        if tag == "ID":
            first = line[2:].strip().split()
            rec_id = first[0].rstrip(";") if first else None
        elif tag == "AC" and not acc:
            acc = body.split(";", 1)[0].strip()
        elif tag == "NM" and not rec_id:
            rec_id = body.strip() or None
        elif tag == "SQ":
            in_seq = True
    if rec_id or seq_chunks:
        yield (rec_id or acc or "unnamed"), "".join(seq_chunks)


def iter_records(handle: TextIO, fmt: str) -> Iterator[Tuple[str, str]]:
    if fmt == "auto":
        fmt = _sniff_format(handle)
    if fmt == "fasta":
        yield from _iter_fasta(handle)
    elif fmt == "embl":
        yield from _iter_embl(handle)
    else:
        raise ValueError("could not detect FASTA/EMBL format")


# --------------------------------------------------------------------------- #
# Member enumeration (file / dir / glob / tar), gzip-aware
# --------------------------------------------------------------------------- #
def _open_member(path: str) -> TextIO:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def iter_source_handles(path: str, pattern: Optional[str]
                        ) -> Iterator[Tuple[str, TextIO]]:
    """Yield (member_name, text_handle) for every file inside a source path."""
    low = path.lower()
    if low.endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(path, "r:*") as tar:
            for m in tar.getmembers():
                if not m.isfile():
                    continue
                if pattern and not _match_ext(m.name, pattern):
                    continue
                if not pattern and not m.name.lower().endswith(_FASTA_EXTS):
                    continue
                fh = tar.extractfile(m)
                if fh is None:
                    continue
                raw = fh.read()
                if m.name.endswith(".gz"):
                    raw = gzip.decompress(raw)
                yield m.name, io.StringIO(raw.decode("utf-8", "replace"))
        return

    if os.path.isdir(path):
        pats = [pattern] if pattern else [f"*{e}" for e in _FASTA_EXTS]
        seen = set()
        for pat in pats:
            for f in sorted(glob(os.path.join(path, pat))):
                if os.path.isfile(f) and f not in seen:
                    seen.add(f)
                    yield os.path.basename(f), _open_member(f)
        if not seen:  # fall back to every file in the directory
            for f in sorted(glob(os.path.join(path, "*"))):
                if os.path.isfile(f):
                    yield os.path.basename(f), _open_member(f)
        return

    if any(ch in path for ch in "*?[") or (pattern and not os.path.exists(path)):
        for f in sorted(glob(path)):
            if os.path.isfile(f):
                yield os.path.basename(f), _open_member(f)
        return

    if not os.path.isfile(path):
        raise SystemExit(f"ERROR: source not found: {path}")
    yield os.path.basename(path), _open_member(path)


def _match_ext(name: str, pattern: str) -> bool:
    import fnmatch
    return fnmatch.fnmatch(os.path.basename(name).lower(), pattern.lower())


# --------------------------------------------------------------------------- #
# Nucleotide check + writer
# --------------------------------------------------------------------------- #
def is_nucleotide(seq: str, min_frac: float = 0.7) -> bool:
    letters = [c for c in seq if c.isalpha()]
    if not letters:
        return False
    n = sum(1 for c in letters if c in _NUCL)
    return (n / len(letters)) >= min_frac


def convert_source(name: str, path: str, out_fa: str, fmt: str,
                   pattern: Optional[str], nucleotide_only: bool,
                   min_len: int) -> Tuple[int, int]:
    written = skipped = 0
    seen_ids: set = set()
    with open(out_fa, "w", encoding="utf-8") as out:
        for member, handle in iter_source_handles(path, pattern):
            try:
                for rec_id, seq in iter_records(handle, fmt):
                    seq = "".join(c for c in seq if c.isalpha())
                    if not seq or len(seq) < min_len:
                        skipped += 1
                        continue
                    if nucleotide_only and not is_nucleotide(seq):
                        skipped += 1
                        continue
                    rid = rec_id or "unnamed"
                    tag = f"{name}__{rid}"
                    if tag in seen_ids:  # guarantee unique target ids
                        suffix = 1
                        while f"{tag}_{suffix}" in seen_ids:
                            suffix += 1
                        tag = f"{tag}_{suffix}"
                    seen_ids.add(tag)
                    out.write(f">{tag}\n{seq}\n")
                    written += 1
            except ValueError as exc:
                LOG.warning("  %s/%s: %s (skipped)", name, member, exc)
            finally:
                try:
                    handle.close()
                except Exception:
                    pass
            LOG.debug("  processed member %s", member)
    return written, skipped


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Convert Dfam (EMBL) and RepBase (.ref) releases into plain "
                    "nucleotide FASTA for leakage_filter.py / MMseqs2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--source", required=True, nargs="+", metavar="NAME=PATH",
                   help="Training source(s). PATH may be a file, directory, "
                        "glob or .tar[.gz]; gzip is fine. e.g. "
                        "dfam=Dfam_curatedonly.embl.gz "
                        "repbase=RepBase29.10.fasta.tar.gz")
    p.add_argument("--outdir", default="training_fasta")
    p.add_argument("--format", default="auto", choices=["auto", "fasta", "embl"],
                   help="Force a record format (default: auto-detect per member).")
    p.add_argument("--pattern", default=None,
                   help="Glob (e.g. '*.ref') to restrict members inside a "
                        "directory/tar. Default: common sequence extensions.")
    p.add_argument("--min-length", type=int, default=1,
                   help="Drop sequences shorter than this (bp).")
    p.add_argument("--keep-nonnucl", action="store_true",
                   help="Keep non-nucleotide (e.g. protein) records too. By "
                        "default only nucleotide sequences are written.")
    p.add_argument("--verbose", action="store_true")
    return p


def parse_sources(items: List[str]) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for it in items:
        if "=" in it:
            name, path = it.split("=", 1)
        else:
            name = os.path.splitext(os.path.basename(it.rstrip("/\\")))[0]
            path = it
        out.append((name.strip(), path.strip()))
    return out


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.ensure_dir(args.outdir)

    produced: List[Tuple[str, str]] = []
    for name, path in parse_sources(args.source):
        out_fa = os.path.join(args.outdir, f"{name}.fasta")
        LOG.info("Converting source '%s' <- %s", name, path)
        written, skipped = convert_source(
            name, path, out_fa, args.format, args.pattern,
            nucleotide_only=not args.keep_nonnucl, min_len=args.min_length)
        LOG.info("  %-10s: %d sequences written, %d skipped -> %s",
                 name, written, skipped, out_fa)
        if written == 0:
            LOG.warning("  '%s' produced 0 sequences — check --format/--pattern.",
                        name)
        produced.append((name, out_fa))

    # Print a ready-to-run leakage_filter.py --training string.
    training = " ".join(f"{n}={p}" for n, p in produced)
    LOG.info("Done. Feed these to leakage_filter.py, e.g.:")
    LOG.info("  python leakage_filter.py --query <panteon_kingdom>.fasta \\")
    LOG.info("      --training %s \\", training)
    LOG.info("      --identity 0.7 --coverage 0.5 --threads 16 "
             "--outdir leakage_out --prefix <kingdom> --write-leaked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
