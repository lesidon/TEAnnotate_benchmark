#!/usr/bin/env python3
"""
panteon_partitioner.py — split the PanTEon Database FASTA into per-kingdom
(Plant / Animal / Fungi) FASTA + label tables for TE-classification benchmarking.

PanTEon headers carry a species binomial but NOT a kingdom, e.g.::

    >PDB00018676#CLASSI/LINE/CR1 @Bucorvus abyssinicus

so each sequence's organism type is resolved from its ``@Species`` via one of
four pluggable, offline-friendly backends (checked in this order per species):

  1. --metadata FILE      PanTEon_Database_metadata_v*.csv — the authoritative
                          ``Kingdom,...,Species,...`` table shipped with the DB.
                          This is the recommended, exact source of truth.
  2. --species-map FILE   explicit TSV/CSV  ``species<TAB>kingdom``  (or
                          ``genus<TAB>kingdom``); manual override, wins over (1).
  3. --taxdump DIR        NCBI taxdump (nodes.dmp+names.dmp) read via `taxopy`.
  4. --ete3               ete3 NCBITaxa (downloads its own db on first use).

Use (1) alone for a faithful split; (3)/(4) only help recover species that are
absent from the metadata (rare — the CSV is the DB's own species list).

The resolved lineage string is classified with te_label_utils.classify_lineage
(Viridiplantae->plant, Metazoa->animal, Fungi->fungi). Species that resolve to
none of the three go to an 'unassigned' partition and are listed in
``<prefix>.unresolved_species.txt`` so you can extend --species-map and re-run.
Every lookup is cached to ``<prefix>.species_kingdom.tsv`` (auditable + reusable
as a future --species-map).

For each kingdom it writes (streaming, O(1) memory per record):
  <prefix>_<kingdom>.fasta   sequences (original PanTEon headers preserved)
  <prefix>_<kingdom>.tsv      labels: family_id, species, label, order, superfamily, ...
and a global <prefix>_partition_summary.tsv.

Example
-------
  python panteon_partitioner.py \
      --fasta PanTEon_Database_v1.5.1_benchmark_edition.fasta \
      --taxdump /path/to/taxdump \
      --outdir panteon_partitions --prefix panteon
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

import te_label_utils as tu

LOG = tu.setup_logger("panteon_partitioner")

TSV_COLUMNS = [
    "family_id", "species", "kingdom", "label",
    "label_class", "label_order", "label_superfamily", "resolved_via",
    "seq_length",
]
KINGDOMS = ("plant", "animal", "fungi")


# --------------------------------------------------------------------------- #
# Streaming FASTA reader (header + sequence)
# --------------------------------------------------------------------------- #
def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    chunks: List[str] = []
    with tu.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


# --------------------------------------------------------------------------- #
# Species -> kingdom resolution
# --------------------------------------------------------------------------- #
class SpeciesResolver:
    """Resolve a species binomial to 'plant'|'animal'|'fungi'|None, with cache.

    Precedence per species: explicit map (species then genus) -> taxonomy
    backend (species then genus). Results are memoised so each unique species is
    looked up once.
    """

    def __init__(self, metadata: Optional[str], species_map: Optional[str],
                 taxdump: Optional[str], use_ete3: bool):
        self._map: Dict[str, str] = {}
        self._cache: Dict[str, Tuple[Optional[str], str]] = {}
        if metadata:
            self._load_metadata(metadata)
        if species_map:  # manual overrides win over the metadata table
            self._load_map(species_map)
        self._taxdb = None
        self._taxopy = None
        self._ncbi = None
        if taxdump:
            self._init_taxopy(taxdump)
        elif use_ete3:
            self._init_ete3()

    # ---- backends -------------------------------------------------------- #
    def _load_metadata(self, path: str) -> None:
        """Load the PanTEon metadata CSV: Species -> Kingdom (authoritative).

        Quoted, multi-column CSV (Kingdom,Phylum,...,Species,...) — parsed with
        the csv module so embedded commas/quotes are safe. Also derives a
        genus->kingdom map (ambiguous genera dropped) to catch species whose
        binomial in the FASTA differs slightly from the metadata row.
        """
        genus_king: Dict[str, str] = {}
        genus_conflict: set = set()
        n_species = 0
        with tu.open_text(path) as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames:
                raise SystemExit(f"ERROR: metadata CSV has no header: {path}")
            sp_col = tu.detect_column(reader.fieldnames, ["species", "organism", "taxon"])
            kg_col = tu.detect_column(reader.fieldnames, ["kingdom", "clade", "group"])
            if not sp_col or not kg_col:
                raise SystemExit(
                    f"ERROR: metadata CSV needs Species+Kingdom columns; "
                    f"found {reader.fieldnames}")
            for row in reader:
                species = (row.get(sp_col) or "").strip()
                king = _canon_kingdom(row.get(kg_col) or "")
                if not species or not king:
                    continue
                self._map.setdefault(species.lower(), king)
                n_species += 1
                genus = species.lower().split()[0]
                if genus in genus_king and genus_king[genus] != king:
                    genus_conflict.add(genus)
                else:
                    genus_king[genus] = king
        for genus, king in genus_king.items():
            if genus not in genus_conflict:
                self._map.setdefault(genus, king)
        LOG.info("Loaded %d species (+%d genera) from metadata %s",
                 n_species, len(genus_king) - len(genus_conflict), path)

    def _load_map(self, path: str) -> None:
        sep = "," if path.lower().endswith(".csv") else "\t"
        with tu.open_text(path) as fh:
            for i, line in enumerate(fh):
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = [p.strip() for p in line.split(sep)]
                if len(parts) < 2:
                    continue
                key, king = parts[0], parts[1].lower()
                if i == 0 and king not in KINGDOMS and \
                        key.lower() in ("species", "organism", "name", "taxon"):
                    continue  # header row
                king = _canon_kingdom(king)
                if king:
                    self._map[key.lower()] = king
        LOG.info("Loaded %d species/genus->kingdom entries from %s",
                 len(self._map), path)

    def _init_taxopy(self, taxdump: str) -> None:
        try:
            import taxopy  # noqa: F401
        except ImportError:
            raise SystemExit("ERROR: --taxdump needs the 'taxopy' package "
                             "(pip install taxopy).")
        nodes = os.path.join(taxdump, "nodes.dmp")
        names = os.path.join(taxdump, "names.dmp")
        if not (os.path.isfile(nodes) and os.path.isfile(names)):
            raise SystemExit(f"ERROR: nodes.dmp/names.dmp not found in {taxdump}")
        LOG.info("Loading NCBI taxdump via taxopy from %s ...", taxdump)
        self._taxopy = taxopy
        self._taxdb = taxopy.TaxDb(nodes_dmp=nodes, names_dmp=names)

    def _init_ete3(self) -> None:
        try:
            from ete3 import NCBITaxa
        except ImportError:
            raise SystemExit("ERROR: --ete3 needs the 'ete3' package "
                             "(pip install ete3).")
        LOG.info("Initialising ete3 NCBITaxa (may download its db on first run)")
        self._ncbi = NCBITaxa()

    def _lineage_blob(self, name: str) -> Optional[str]:
        """Return a space-joined lineage string for a taxon name, or None."""
        if self._taxdb is not None:
            try:
                taxids = self._taxopy.taxid_from_name(name, self._taxdb)
            except Exception:
                taxids = []
            if not taxids:
                return None
            try:
                taxon = self._taxopy.Taxon(taxids[0], self._taxdb)
                return " ".join(taxon.name_lineage)
            except Exception:
                return None
        if self._ncbi is not None:
            try:
                mapping = self._ncbi.get_name_translator([name])
                if not mapping.get(name):
                    return None
                lineage = self._ncbi.get_lineage(mapping[name][0])
                names = self._ncbi.get_taxid_translator(lineage)
                return " ".join(names[t] for t in lineage)
            except Exception:
                return None
        return None

    @property
    def has_backend(self) -> bool:
        return bool(self._map) or self._taxdb is not None or self._ncbi is not None

    # ---- public API ------------------------------------------------------ #
    def resolve(self, species: str) -> Tuple[Optional[str], str]:
        """Return (kingdom|None, resolved_via) for a species binomial."""
        key = (species or "").strip()
        if not key:
            return (None, "no_species")
        if key in self._cache:
            return self._cache[key]
        low = key.lower()
        genus = low.split()[0] if low.split() else low

        result: Tuple[Optional[str], str] = (None, "unresolved")
        if low in self._map:
            result = (self._map[low], "map:species")
        elif genus in self._map:
            result = (self._map[genus], "map:genus")
        else:
            blob = self._lineage_blob(key)
            if blob is None and genus != low:
                blob = self._lineage_blob(genus)
            if blob:
                king = tu.classify_lineage(blob)
                if king:
                    result = (king, "taxonomy")
        self._cache[key] = result
        return result

    def dump_cache(self, path: str) -> None:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["species", "kingdom", "resolved_via"])
            for sp, (king, via) in sorted(self._cache.items()):
                w.writerow([sp, king or "unassigned", via])


def _canon_kingdom(value: str) -> Optional[str]:
    v = (value or "").strip().lower()
    if v in ("plant", "plants", "plantae", "viridiplantae"):
        return "plant"
    if v in ("animal", "animals", "animalia", "metazoa"):
        return "animal"
    if v in ("fungi", "fungus", "fungal"):
        return "fungi"
    return None


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
class KingdomWriters:
    def __init__(self, outdir: str, prefix: str, write_fasta: bool):
        self.outdir = tu.ensure_dir(outdir)
        self.prefix = prefix
        self.write_fasta = write_fasta
        self._fa: Dict[str, TextIO] = {}
        self._tsv: Dict[str, TextIO] = {}
        self._tsv_w: Dict[str, "csv.writer"] = {}

    def _fasta(self, kingdom: str) -> TextIO:
        if kingdom not in self._fa:
            p = os.path.join(self.outdir, f"{self.prefix}_{kingdom}.fasta")
            self._fa[kingdom] = open(p, "w", encoding="utf-8")
        return self._fa[kingdom]

    def _table(self, kingdom: str) -> "csv.writer":
        if kingdom not in self._tsv_w:
            p = os.path.join(self.outdir, f"{self.prefix}_{kingdom}.tsv")
            fh = open(p, "w", encoding="utf-8", newline="")
            w = csv.writer(fh, delimiter="\t")
            w.writerow(TSV_COLUMNS)
            self._tsv[kingdom] = fh
            self._tsv_w[kingdom] = w
        return self._tsv_w[kingdom]

    def write(self, kingdom: str, header: str, sequence: str,
              row: List[object]) -> None:
        self._table(kingdom).writerow(row)
        if self.write_fasta and sequence:
            fh = self._fasta(kingdom)
            fh.write(f">{header}\n")
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
        description="Split the PanTEon Database FASTA into per-kingdom "
                    "(plant/animal/fungi) FASTA + label tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--fasta", required=True,
                   help="PanTEon Database FASTA (.fasta or .fasta.gz).")
    p.add_argument("--outdir", default="panteon_partitions")
    p.add_argument("--prefix", default="panteon")
    p.add_argument("--metadata",
                   help="PanTEon_Database_metadata_v*.csv (authoritative "
                        "Kingdom,...,Species,... table). Recommended resolver.")
    p.add_argument("--species-map",
                   help="TSV/CSV of species-or-genus<TAB>kingdom overrides "
                        "(wins over --metadata).")
    p.add_argument("--taxdump",
                   help="NCBI taxdump directory (nodes.dmp+names.dmp) for "
                        "taxopy-based resolution.")
    p.add_argument("--ete3", action="store_true",
                   help="Use ete3 NCBITaxa for resolution (auto-downloads db).")
    p.add_argument("--min-length", type=int, default=0,
                   help="Drop sequences shorter than this (bp).")
    p.add_argument("--keep-nonte", action="store_true",
                   help="Keep sequences whose label normalises to Unknown/Other "
                        "(dropped by default).")
    p.add_argument("--keep-unassigned", action="store_true",
                   help="Emit an 'unassigned' partition for species that resolve "
                        "to no kingdom (otherwise only counted + listed).")
    p.add_argument("--label-map", help="JSON of extra label synonyms.")
    p.add_argument("--no-fasta", action="store_true",
                   help="Write only the label TSVs, not per-kingdom FASTA.")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.verbose:
        LOG.setLevel("DEBUG")
    tu.load_label_overrides(args.label_map)
    tu.ensure_dir(args.outdir)

    resolver = SpeciesResolver(args.metadata, args.species_map,
                               args.taxdump, args.ete3)
    if not resolver.has_backend:
        raise SystemExit(
            "ERROR: no species->kingdom backend configured. Provide at least one "
            "of --metadata (recommended), --species-map, --taxdump or --ete3 "
            "(see --help).")

    writers = KingdomWriters(args.outdir, args.prefix,
                             write_fasta=not args.no_fasta)
    counts: "Counter[Tuple[str, str]]" = Counter()
    stats: "Counter[str]" = Counter()
    unresolved_species: "Counter[str]" = Counter()

    try:
        for header, sequence in iter_fasta(args.fasta):
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
            if args.min_length and len(sequence) < args.min_length:
                stats["too_short"] += 1
                continue

            kingdom, via = resolver.resolve(species)
            if kingdom is None:
                stats["unassigned"] += 1
                if species:
                    unresolved_species[species] += 1
                if not args.keep_unassigned:
                    continue
                kingdom = "unassigned"

            row = [fam_id, species, kingdom, raw_label, te_class, order,
                   superfam, via, len(sequence)]
            writers.write(kingdom, header, sequence, row)
            counts[(kingdom, order)] += 1
            stats["written"] += 1
            if stats["total"] % 50000 == 0:
                LOG.info("... parsed %d sequences (%d written)",
                         stats["total"], stats["written"])
    finally:
        writers.close()

    # Species->kingdom cache (auditable; reusable as a future --species-map).
    resolver.dump_cache(os.path.join(args.outdir, f"{args.prefix}.species_kingdom.tsv"))

    # Unresolved species list (extend --species-map and re-run to recover these).
    if unresolved_species:
        upath = os.path.join(args.outdir, f"{args.prefix}.unresolved_species.txt")
        with open(upath, "w", encoding="utf-8") as fh:
            for sp, n in unresolved_species.most_common():
                fh.write(f"{sp}\t{n}\n")
        LOG.warning("%d sequences from %d species could not be assigned to a "
                    "kingdom; see %s", stats["unassigned"],
                    len(unresolved_species), upath)

    # Partition summary.
    summ = os.path.join(args.outdir, f"{args.prefix}_partition_summary.tsv")
    with open(summ, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["kingdom", "order", "n_sequences"])
        for (kingdom, order), n in sorted(counts.items()):
            w.writerow([kingdom, order, n])

    LOG.info("Done. total=%d written=%d "
             "(no_id=%d nonte_dropped=%d too_short=%d unassigned=%d)",
             stats["total"], stats["written"], stats["no_id"],
             stats["nonte_dropped"], stats["too_short"], stats["unassigned"])
    LOG.info("Per-kingdom FASTA + label TSV written to %s/", args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
