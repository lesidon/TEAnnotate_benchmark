#!/usr/bin/env bash
# prepull_containers.sh — pull the DL classifier images into the Singularity cache
# ONCE, with retries, so the pipeline reuses them instead of pulling mid-run.
#
# Why: the classifier images are served from GHCR over ORAS, whose HTTP/2 stream
# occasionally dies with "stream error ... PROTOCOL_ERROR; received from peer".
# A pull that fails inside a run is FATAL — Nextflow raises it during task staging,
# before the task starts, so the process-level `errorStrategy 'ignore'` cannot catch
# it and the whole pipeline aborts. Pre-pulling (serially, with retries) into a
# persistent cache sidesteps the problem: once the .img is cached, Nextflow's
# resolver (sifCached() in nextflow.config) sees it and never pulls again.
#
# Usage:
#   bin/prepull_containers.sh [-c CACHE_DIR] [-n RETRIES] [image ...]
#     -c CACHE_DIR   cache directory (default: $NXF_SINGULARITY_CACHEDIR,
#                    else $SINGULARITY_CACHEDIR, else ./singularity_cache)
#     -n RETRIES     attempts per image (default 5)
#   Extra positional args replace the default image list.
#
# Then run the pipeline pointed at the same directory:
#   export NXF_SINGULARITY_CACHEDIR=<CACHE_DIR>       # or --singularity_cache_dir <CACHE_DIR>
#   nextflow run . -profile singularity ...
set -euo pipefail

CACHE="${NXF_SINGULARITY_CACHEDIR:-${SINGULARITY_CACHEDIR:-$PWD/singularity_cache}}"
RETRIES=5

while getopts ":c:n:h" opt; do
    case "$opt" in
        c) CACHE="$OPTARG" ;;
        n) RETRIES="$OPTARG" ;;
        h) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "invalid option; run with -h" >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))

# Default image set = EVERY registry image the pipeline pulls (keep in sync with
# conf/modules.config + the *_container params). Pre-pulling all of them on a
# networked node means the compute nodes need no registry access at all — the fix
# for the intermittent "read tcp ...:443: connection timed out" blob-download
# failures (e.g. dfam/tetools). Override by passing images as args.
# NB: conda-only steps (EARLGREY, HELIANO, CD_HIT_EST, PFAM_FETCH, STARFISH, the
# *_CLONE/*_FETCH_MODEL git/curl envs) have no image — pre-build those conda envs
# into a shared $NXF_CONDA_CACHEDIR on a networked node as well.
IMAGES=("$@")
if [[ ${#IMAGES[@]} -eq 0 ]]; then
    IMAGES=(
        # DL classifier images (GHCR / ORAS)
        oras://ghcr.io/lesidon/teannotate-teclass2:1.0
        oras://ghcr.io/lesidon/teannotate-create:1.0
        oras://ghcr.io/lesidon/teannotate-deepte:1.0
        oras://ghcr.io/lesidon/teannotate-neuralte:1.0
        oras://ghcr.io/lesidon/teannotate-terrier:1.0
        # de novo / masking / QC images (docker registries)
        docker://dfam/tetools:1.89.2
        docker://kanghu/hite:3.3.3
        docker://quay.io/biocontainers/edta:2.3.0--hdfd78af_0
        docker://quay.io/biocontainers/seqkit:2.8.0--h9ee0642_0
        docker://quay.io/biocontainers/tetrimmer:1.7.2--hdfd78af_0
        docker://rocker/tidyverse:4.3.2
        docker://docker.io/python:3.11
    )
fi

ENGINE="$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || true)"
[[ -n "$ENGINE" ]] || { echo "ERROR: neither singularity nor apptainer on PATH" >&2; exit 1; }
mkdir -p "$CACHE"

# Replicate Nextflow's SingularityCache.simpleName(): drop protocol, ':'/'/' -> '-'.
# Nextflow then appends '.img' (confirmed by the *.img.pulling.* temp name it uses).
cache_name() { local s="${1#*://}"; s="${s//:/-}"; s="${s//\//-}"; printf '%s.img' "$s"; }

fail=0
for img in "${IMAGES[@]}"; do
    out="${CACHE}/$(cache_name "$img")"
    if [[ -s "$out" ]]; then
        echo "[cached]  $out"
        continue
    fi
    echo "[pull]    $img"
    n=0
    until "$ENGINE" pull --name "$out" "$img"; do
        n=$((n + 1))
        rm -f "$out"
        if [[ $n -ge $RETRIES ]]; then
            echo "[FAILED]  $img  (after ${RETRIES} attempts)" >&2
            fail=1
            break
        fi
        echo "[retry ${n}/${RETRIES}] $img in 15s ..." >&2
        sleep 15
    done
    [[ -s "$out" ]] && echo "[ok]      $out"
done

echo
if [[ $fail -eq 0 ]]; then
    echo "All images cached in: $CACHE"
else
    echo "Some images failed to pull (see above). Re-run to retry the missing ones." >&2
fi
echo "Use this cache in the run with:"
echo "  export NXF_SINGULARITY_CACHEDIR=\"$CACHE\"   # or --singularity_cache_dir \"$CACHE\""
exit $fail
