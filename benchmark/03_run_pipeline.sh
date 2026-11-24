#!/usr/bin/env bash
# 03_run_pipeline.sh — run the FULL nf-core/teannotate pipeline on each genome and
# capture its final consensus library as the "pipeline" method.
#
#   libraries/<species>/pipeline/pipeline_library.fa
#
# The pipeline's final library is CURATION.out.library -> *.classified_library.fa
# (the file that feeds RepeatMasker/MASKING). We copy that out of the run's outdir.
#
# Usage:
#   NF_PIPELINE=~/nf-core-te-annotation ./03_run_pipeline.sh                 # all
#   NF_PIPELINE=~/nf-core-te-annotation ./03_run_pipeline.sh oryza_sativa    # one
#
# Heavy: each genome is a complete de-novo + curation + DL-ensemble run (hours to
# days on large genomes). -resume makes it re-entrant across restarts.
set -euo pipefail
source "$(dirname "$0")/env.sh"

ONLY="${1:-}"
[[ -n "$NF_PIPELINE" ]] || die "set NF_PIPELINE to the nf-core/teannotate checkout (contains main.nf)"
[[ -f "${NF_PIPELINE}/main.nf" ]] || die "no main.nf under NF_PIPELINE=${NF_PIPELINE}"
command -v nextflow >/dev/null 2>&1 || die "nextflow not on PATH"

run_one() {
    local kingdom="$1" species="$2" organism="$3"
    local genome; genome="$(genome_fasta "$kingdom" "$species")"
    local out="${LIB_DIR}/${species}/pipeline"
    local dest="${out}/pipeline_library.fa"

    nonempty "$genome" || { log "WARN ${species}: genome missing ($genome) — run 01 first"; return 0; }
    if nonempty "$dest"; then log "skip ${species}: pipeline library already captured"; return 0; fi

    local rundir="${WORK_DIR}/pipeline/${species}"
    local outdir="${rundir}/out"
    mkdir -p "$rundir" "$out"

    # One-row samplesheet: sample,genome,organism_type,gff  (gff optional/fungal).
    local sheet="${rundir}/samplesheet.csv"
    printf 'sample,genome,organism_type,gff\n%s,%s,%s,\n' \
        "$species" "$(readlink -f "$genome")" "$organism" > "$sheet"

    log "pipeline run ${species} (${organism}) ..."
    # Shared, pre-stageable caches so images/conda envs are reused across genomes and
    # the compute node needs no registry/conda internet (pre-fill on a login node).
    ( cd "$rundir" && nextflow run "${NF_PIPELINE}/main.nf" \
        -profile singularity \
        -work-dir "${rundir}/work" \
        --input "$sheet" \
        --outdir "$outdir" \
        --max_cpus 16 \
        --max_memory 256GB \
        -resume ) || die "nextflow failed for ${species}"

    # Capture the final classified/curated consensus library from the outdir.
    local lib
    lib="$(find "$outdir" -type f \( -name '*.classified_library.fa' -o -name '*.curated_library.fa' \) \
             -print -quit 2>/dev/null || true)"
    [[ -n "$lib" && -s "$lib" ]] || die "no *_library.fa found under ${outdir} for ${species}"
    cp "$lib" "$dest"
    log "pipeline library captured: $(grep -c '^>' "$dest") sequences -> $dest"
}

read_manifest "$GENOMES_TSV" | while IFS=$'\t' read -r kingdom species display accession organism scoring curated; do
    [[ -n "$ONLY" && "$species" != "$ONLY" ]] && continue
    run_one "$kingdom" "$species" "$organism"
done

log "pipeline stage complete."
