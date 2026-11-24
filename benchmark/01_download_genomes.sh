#!/usr/bin/env bash
# 01_download_genomes.sh — pull each assembly from NCBI by accession, verify the
# package checksums, log provenance, and store under genomes/<kingdom>/<species>/.
#
# Usage:
#   ./01_download_genomes.sh                 # all genomes in the manifest
#   ./01_download_genomes.sh oryza_sativa    # just one species_key
#
# Requires the NCBI 'datasets' CLI (conda install -c conda-forge ncbi-datasets-cli).
# Re-entrant: a species with a non-empty <species>.fna and provenance.json is skipped.
set -euo pipefail
source "$(dirname "$0")/env.sh"

ONLY="${1:-}"

command -v datasets >/dev/null 2>&1 || die \
    "NCBI 'datasets' CLI not found. Install: conda install -c conda-forge ncbi-datasets-cli"

download_one() {
    local kingdom="$1" species="$2" accession="$3" display="$4"
    local dest="${GENOME_DIR}/${kingdom}/${species}"
    local fna="${dest}/${species}.fna"
    local prov="${dest}/provenance.json"

    if nonempty "$fna" && nonempty "$prov"; then
        log "skip ${species}: already downloaded ($fna)"; return 0
    fi
    mkdir -p "$dest"

    local tmp; tmp="$(mktemp -d "${WORK_DIR}/dl_${species}.XXXXXX")"
    local zip="${tmp}/${species}.zip"
    # Binomial for the fallback = display name with any "(assembly)" suffix stripped.
    local binomial; binomial="$(echo "$display" | sed 's/ *(.*//; s/[[:space:]]*$//')"
    local source_desc used_taxon=0

    log "download ${species} (${accession}) ..."
    if datasets download genome accession "$accession" --include genome --filename "$zip" 2>/dev/null; then
        source_desc="NCBI datasets (accession ${accession})"
    else
        # NCBI suppresses/supersedes assemblies (e.g. GRCz11 -> GRCz12); fall back to
        # the current reference genome for the species so the run self-heals.
        log "WARN ${species}: accession ${accession} unavailable (suppressed/superseded?)"
        log "       falling back to the current NCBI reference for '${binomial}'"
        datasets download genome taxon "$binomial" --reference --include genome --filename "$zip" \
            || die "download failed for ${accession} AND taxon reference '${binomial}'"
        used_taxon=1
        source_desc="NCBI datasets (taxon '${binomial}' --reference; requested ${accession} was unavailable)"
    fi

    ( cd "$tmp" && unzip -q -o "$zip" )

    # Verify the package checksums NCBI ships (md5sum.txt at the archive root).
    if [[ -f "${tmp}/md5sum.txt" ]]; then
        ( cd "$tmp" && md5sum -c --quiet md5sum.txt ) \
            || die "checksum mismatch for ${species} — download corrupt"
        log "checksum OK (${species})"
    else
        log "WARN ${species}: no md5sum.txt in package; skipping checksum verify"
    fi

    # Discover the ACTUAL accession directory (unknown when the taxon fallback ran).
    local datadir="${tmp}/ncbi_dataset/data" actual
    actual="$(find "$datadir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null \
                | grep -E '^GC[AF]_' | head -1 || true)"
    [[ -n "$actual" ]] || actual="$accession"
    if [[ "$used_taxon" == "1" && "$actual" != "$accession" ]]; then
        log "       taxon fallback resolved to ${actual}"
    fi

    # Locate the assembled genome FASTA inside the datasets layout.
    local src
    src="$(find "${datadir}/${actual}" \
              \( -name '*.fna' -o -name '*.fasta' -o -name '*.fa' \) -print -quit 2>/dev/null || true)"
    [[ -n "$src" && -s "$src" ]] || die "no genome FASTA found for ${species} (${actual})"
    cp "$src" "$fna"

    # Provenance: the actual accession downloaded, assembly name, md5, date, source.
    local md5 asmname
    md5="$(md5sum "$fna" | awk '{print $1}')"
    asmname="$(datasets summary genome accession "$actual" 2>/dev/null \
                | grep -oE '"assembly_name":[^,]*' | head -1 | sed 's/.*: *"//; s/"$//' || true)"
    cat > "$prov" <<JSON
{
  "species_key": "${species}",
  "requested_accession": "${accession}",
  "accession": "${actual}",
  "assembly_name": "${asmname:-unknown}",
  "source": "${source_desc}",
  "genome_file": "${species}.fna",
  "md5": "${md5}",
  "downloaded": "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
}
JSON
    rm -rf "$tmp"
    log "done ${species}: ${fna} (accession ${actual}, md5 ${md5})"
}

read_manifest "$GENOMES_TSV" | while IFS=$'\t' read -r kingdom species display accession organism scoring curated; do
    [[ -n "$ONLY" && "$species" != "$ONLY" ]] && continue
    download_one "$kingdom" "$species" "$accession" "$display"
done

log "all downloads complete."
