#!/usr/bin/env python3
"""
Standalone samplesheet validator for nf-core/teannotate.

Performs checks beyond what nf-schema's JSON Schema validation covers:
  - Duplicate sample IDs
  - Recognised FASTA file extensions
  - Valid organism_type values

Outputs a normalised CSV with default values filled in.
Can be run directly by a user before launching the pipeline:

    python bin/check_samplesheet.py --input samplesheet.csv --output validated.csv
"""
import argparse
import csv
import sys
from pathlib import Path

VALID_ORGANISM_TYPES = {"plant", "animal", "fungus"}
FASTA_EXTENSIONS     = {".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"}


def validate(input_path: str, output_path: str) -> None:
    errors  = []
    samples = []
    seen_ids: set = set()

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []

        required = {"sample", "genome", "organism_type"}
        missing  = required - set(fieldnames)
        if missing:
            sys.exit(
                f"ERROR: samplesheet is missing required columns: "
                f"{', '.join(sorted(missing))}"
            )

        for lineno, row in enumerate(reader, start=2):
            sid      = row.get("sample",        "").strip()
            genome   = row.get("genome",         "").strip()
            org_type = row.get("organism_type",  "animal").strip().lower() or "animal"

            # ── Duplicate sample IDs ───────────────────────────────────────────
            if not sid:
                errors.append(f"Line {lineno}: 'sample' column is empty.")
            elif sid in seen_ids:
                errors.append(f"Line {lineno}: duplicate sample ID '{sid}'.")
            seen_ids.add(sid)

            # ── Genome path: extension check only (existence checked by nf-schema) ─
            if not genome:
                errors.append(f"Line {lineno} [{sid}]: 'genome' column is empty.")
            else:
                # Normalise: handle double extensions like .fa.gz
                name = genome.lower()
                if not any(name.endswith(ext) for ext in FASTA_EXTENSIONS):
                    errors.append(
                        f"Line {lineno} [{sid}]: 'genome' file '{genome}' does not "
                        f"have a recognised FASTA extension "
                        f"({', '.join(sorted(FASTA_EXTENSIONS))})."
                    )

            # ── organism_type ──────────────────────────────────────────────────
            if org_type not in VALID_ORGANISM_TYPES:
                errors.append(
                    f"Line {lineno} [{sid}]: organism_type='{org_type}' is invalid. "
                    f"Must be one of: {', '.join(sorted(VALID_ORGANISM_TYPES))}."
                )

            samples.append({
                "sample"       : sid,
                "genome"       : genome,
                "organism_type": org_type or "animal",
            })

    if errors:
        print("ERROR: samplesheet validation failed:", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample", "genome", "organism_type"],
        )
        writer.writeheader()
        writer.writerows(samples)

    print(
        f"Samplesheet OK: {len(samples)} sample(s) passed validation -> '{output_path}'",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Validate nf-core/teannotate samplesheet CSV."
    )
    parser.add_argument("--input",  required=True, help="Input samplesheet CSV")
    parser.add_argument("--output", required=True, help="Output validated CSV")
    args = parser.parse_args()
    validate(args.input, args.output)


if __name__ == "__main__":
    main()
