#!/usr/bin/env bash
# 00_prestage.sh — ONE-TIME, SERIAL pre-staging on a node WITH internet (e.g. an01).
#
# Populates the shared caches so the compute-node array (submit_benchmark.slurm)
# neither pulls images nor builds conda envs. This avoids BOTH failure modes seen:
#   * "read tcp ...:443: connection timed out" / "CONNECTION FAILED conda.anaconda.org"
#     — compute nodes have no reliable internet;
#   * "java.nio FileAlreadyExistsException: .../<image>.img"
#     — concurrent `nextflow run`s racing to pull the SAME image into the shared
#       cacheDir (Nextflow does not lock the image cache across separate runs).
#
# Run this to COMPLETION before `sbatch submit_benchmark.slurm`. Everything here is
# serial, so nothing races. Prereqs: 01 + 02 done, NF_PIPELINE set. Heavy — run it
# under tmux/screen.
set -euo pipefail
source "$(dirname "$0")/env.sh"

[[ -n "$NF_PIPELINE" ]] || die "set NF_PIPELINE to the nf-core/teannotate checkout (contains main.nf)"
[[ -f "${NF_PIPELINE}/bin/prepull_containers.sh" ]] || die "no bin/prepull_containers.sh under NF_PIPELINE"

# Clear any half-pulled image locks left by an interrupted/raced pull.
log "clearing stale pull locks in $NXF_SINGULARITY_CACHEDIR"
find "$NXF_SINGULARITY_CACHEDIR" -maxdepth 1 -name '*.pulling.*' -delete 2>/dev/null || true

# 1. Pull EVERY pipeline image serially into the shared Nextflow image cache. Once
#    cached, no `nextflow run` pulls it, so the array tasks can never race on it.
log "pre-pulling all pipeline images -> $NXF_SINGULARITY_CACHEDIR"
bash "${NF_PIPELINE}/bin/prepull_containers.sh" -c "$NXF_SINGULARITY_CACHEDIR"

# 2. Build the conda-only envs (EarlGrey / HELIANO / CD-HIT / Pfam / clone-fetch) into
#    the shared conda cache by warming up ONE genome per de-novo route, SERIALLY. The
#    fungus route covers EarlGrey/CD-HIT/Pfam; the plant route adds EDTA/HELIANO; the
#    animal route reuses EarlGrey. These warm-ups double as the real pipeline runs for
#    those two genomes (03 is re-entrant, so the array later skips them).
#    Set PRESTAGE_CONDA=0 to skip (e.g. if your compute nodes CAN reach conda).
if [[ "${PRESTAGE_CONDA:-1}" == "1" ]]; then
    log "conda warm-up (fungus route) — saccharomyces_cerevisiae"
    bash ./03_run_pipeline.sh saccharomyces_cerevisiae
    log "conda warm-up (plant route) — arabidopsis_thaliana"
    bash ./03_run_pipeline.sh arabidopsis_thaliana
fi

log "pre-stage complete. Caches:"
log "  images: $NXF_SINGULARITY_CACHEDIR"
log "  conda : $NXF_CONDA_CACHEDIR"
log "Now submit the compute array:  sbatch submit_benchmark.slurm"
