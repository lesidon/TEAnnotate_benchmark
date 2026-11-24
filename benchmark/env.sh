#!/usr/bin/env bash
# env.sh — shared configuration + helpers for the whole-pipeline TE benchmark.
# Source it at the top of every stage script:  source "$(dirname "$0")/env.sh"
#
# Tuned for Zhores (Singularity/Apptainer + Conda, Slurm, read-only $HOME). Every
# value can be overridden from the environment before sourcing, e.g.
#   NF_PIPELINE=~/nf-core-te-annotation THREADS=32 ./04_run_tools.sh oryza_sativa
#
# NOTE: do not `set -e` here — this file is *sourced*; a non-zero from an optional
# probe (module/command -v) must not kill the caller. Stage scripts set their own
# `set -euo pipefail`.

# ── Locations ──────────────────────────────────────────────────────────────
BENCH_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"   # -P: physical path
        # (resolve symlinks). Inputs are passed to containers via `readlink -f` (physical),
        # so BENCH_ROOT must be physical too or the `-B "$BENCH_ROOT"` bind won't cover them
        # (e.g. a /trinity symlinked home whose files really live on /beegfs).
export BENCH_ROOT
export CONFIG_DIR="${BENCH_ROOT}/config"
export GENOME_DIR="${GENOME_DIR:-${BENCH_ROOT}/genomes}"
export CURATED_DIR="${CURATED_DIR:-${BENCH_ROOT}/curated}"
export LIB_DIR="${LIB_DIR:-${BENCH_ROOT}/libraries}"
export MASK_DIR="${MASK_DIR:-${BENCH_ROOT}/masking}"
export RESULT_DIR="${RESULT_DIR:-${BENCH_ROOT}/results}"
export FIG_DIR="${FIG_DIR:-${RESULT_DIR}/figures}"
export LOG_DIR="${LOG_DIR:-${BENCH_ROOT}/logs}"
export WORK_DIR="${WORK_DIR:-${BENCH_ROOT}/work}"       # nextflow + scratch (put on beegfs)
export TOOLS_DIR="${TOOLS_DIR:-${BENCH_ROOT}/tools}"    # cloned EDTA (lib-test.py + curated libs)
mkdir -p "$GENOME_DIR" "$CURATED_DIR" "$LIB_DIR" "$MASK_DIR" \
         "$RESULT_DIR" "$FIG_DIR" "$LOG_DIR" "$WORK_DIR" "$TOOLS_DIR"

export GENOMES_TSV="${CONFIG_DIR}/genomes.tsv"
export METHODS_TSV="${CONFIG_DIR}/methods.tsv"

# ── The pipeline checkout (main.nf + nextflow.config) ──────────────────────
# Must point at the nf-core/teannotate repo. No sane default across machines, so
# leave empty and let 03_run_pipeline.sh error clearly if unset.
export NF_PIPELINE="${NF_PIPELINE:-}"

# ── Compute ────────────────────────────────────────────────────────────────
export THREADS="${THREADS:-16}"
# LAI (LTR_retriever) is per-genome context and expensive on large genomes.
# Set RUN_LAI=0 to skip it globally.
export RUN_LAI="${RUN_LAI:-1}"

# ── Containers / conda envs (pin tags to match the pipeline) ───────────────
export EDTA_SIF="${EDTA_SIF:-docker://quay.io/biocontainers/edta:2.3.0--hdfd78af_0}"
export HITE_SIF="${HITE_SIF:-docker://kanghu/hite:3.3.3}"
export RM_SIF="${RM_SIF:-docker://dfam/tetools:1.89.2}"      # RepeatMasker (matches MASKING stage)
export EARLGREY_ENV="${EARLGREY_ENV:-eg_bench}"             # conda env name (built from the pipeline yml)
# EarlGrey via prebuilt container — set to skip the (intractable) conda solve. Use a FULL
# (Dfam-preconfigured) image, NOT a '-nodfam' one: earlGreyLibConstruct hard-checks for a
# configured Dfam partition 0 (dfam39_full.0.h5 + rmlib.config in RepeatMasker's famdb dir)
# and, if absent, drops into an interactive "biocontainer?" prompt that hangs on
# non-interactive stdin (endless "Please answer yes or no."). The full image ships partition
# 0 so the check passes with no prompt, and works offline. run_earlgrey uses this when set.
#   export EARLGREY_SIF=docker://tobybaril/earlgrey:v6.3.4
export EARLGREY_SIF="${EARLGREY_SIF:-}"
export HELIANO_ENV="${HELIANO_ENV:-heliano_bench}"          # conda env name (bioconda::heliano)
# EDTA repo checkout — ships lib-test.py + database/{rice7.0.0.liban,maizeTE02052020}.
export EDTA_REPO="${EDTA_REPO:-${TOOLS_DIR}/EDTA}"
export EDTA_GIT="${EDTA_GIT:-https://github.com/oushujun/EDTA.git}"

# Standardized RepeatMasker params — identical for the reference and every test
# library so the .out files (and lib-test.py) are comparable. Mirrors EDTA's
# benchmarking RM call: NCBI engine, TE-only (no simple/low-complexity/rna).
export RM_OPTS="${RM_OPTS:- -e ncbi -no_is -norna -nolow -div 40 }"

# ── Container engine + cache ───────────────────────────────────────────────
# Shared Singularity .img cache, REUSED across runs. On HPCs where the compute (and often
# login) nodes have flaky registry access — mid-run pulls die with "TLS handshake timeout"
# / "connection timed out" fetching quay.io blobs — pre-pull EVERY image into this dir on a
# networked node so runs never pull:
#     source ./env.sh && bash "$NF_PIPELINE/bin/prepull_containers.sh"   # uses $NXF_SINGULARITY_CACHEDIR
# Set explicitly (not via :-) so a stale exported value can't send runs at a wrong cache;
# a failed mid-run pull is FATAL (raised at staging, uncatchable by errorStrategy).
export NXF_SINGULARITY_CACHEDIR="${WORK_DIR}/singularity_cache"
# Shared conda-env cache, REUSED across runs so the heavy EarlGrey env (and the others) is
# built ONCE, not per genome. Requires `jq` on PATH when envs build (a bioconductor data
# package's post-link script needs it): `conda install -n base -c conda-forge jq`.
# A FAILED/interrupted build leaves a poisoned partial env that blocks the next run
# ("Could not remove or rename ..."), so if a build breaks, delete that one env dir under
# this cache before retrying. Set explicitly so a stale value can't point at a wrong cache.
export NXF_CONDA_CACHEDIR="${WORK_DIR}/conda_cache"

CONTAINER_ENGINE="${CONTAINER_ENGINE:-$(command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null || echo singularity)}"
export CONTAINER_ENGINE
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR}}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SINGULARITY_CACHEDIR}}"
mkdir -p "$SINGULARITY_CACHEDIR"

# Extra bind mounts (colon/space separated -B args) — e.g. beegfs scratch.
export SING_BINDS="${SING_BINDS:-}"

# Best-effort module loads (names vary; never fatal on a sourced file).
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true
module load nextflow    2>/dev/null || true

# Nextflow 26.04+ makes the STRICT config parser the default, which forbids `def` function
# declarations in nextflow.config — but the pipeline defines helper functions there
# (resolveContainer/resolveConda/singularityCacheName/…). Force the legacy (v1) parser so
# the config parses. Offline-safe: no new Nextflow download (unlike pinning NXF_VER).
export NXF_SYNTAX_PARSER=v1

# ── Helpers ────────────────────────────────────────────────────────────────
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2; exit 1; }

# container_run <image> <cmd...> — writable-tmpfs so bundled tools can write, and
# bind BENCH_ROOT (+ any SING_BINDS) so inputs/outputs are visible.
container_run() {
    local img="$1"; shift
    local binds=( -B "$BENCH_ROOT" )
    # Also bind the *physical* directory of any existing absolute-path argument. Inputs
    # are passed as `readlink -f` paths, which can resolve onto a filesystem (e.g. /beegfs)
    # outside the symlinked BENCH_ROOT; without this the file is invisible in the container
    # and EDTA/HiTE fail with a misleading "input not found". Dedup keeps the arg list short.
    local a d seen=":"
    for a in "$@"; do
        if [[ "$a" == /* && -e "$a" ]]; then
            d="$(cd "$(dirname "$a")" 2>/dev/null && pwd -P)" || continue
            [[ "$seen" == *":${d}:"* ]] || { binds+=( -B "$d" ); seen="${seen}${d}:"; }
        fi
    done
    # shellcheck disable=SC2086
    "$CONTAINER_ENGINE" exec --writable-tmpfs "${binds[@]}" ${SING_BINDS} "$img" "$@"
}

# conda_run <env> <cmd...> — run inside a named conda env without a login shell.
conda_run() {
    local envname="$1"; shift
    # Accept either a named env (-n) or an env *prefix path* (-p). A value containing '/'
    # (or an existing directory) is treated as a prefix — this lets EARLGREY_ENV / HELIANO_ENV
    # point at a Nextflow-built env from the pipeline's work/conda cache, e.g.
    #   export EARLGREY_ENV=/.../work/pipeline/<sp>/work/conda/env-<hash>
    # so we reuse an already-solved env instead of rebuilding a slow one.
    local sel=(-n "$envname")
    [[ "$envname" == */* || -d "$envname" ]] && sel=(-p "$envname")
    if command -v conda >/dev/null 2>&1; then
        conda run --no-capture-output "${sel[@]}" "$@"
    elif command -v micromamba >/dev/null 2>&1; then
        micromamba run "${sel[@]}" "$@"
    else
        die "neither conda nor micromamba on PATH (needed for env '$envname')"
    fi
}

# read_manifest <tsv> — emit non-comment, non-blank rows (fields tab-separated).
read_manifest() { grep -vE '^[[:space:]]*(#|$)' "$1"; }

# genome_row <species_key> — echo the manifest row for one species (tab-separated).
genome_row() { read_manifest "$GENOMES_TSV" | awk -F'\t' -v k="$1" '$2==k{print; exit}'; }

# genome_fasta <kingdom> <species_key> — canonical path of the downloaded genome.
genome_fasta() { echo "${GENOME_DIR}/$1/$2/$2.fna"; }

# nonempty <file> — true if file exists and is non-empty (re-entrancy guard).
nonempty() { [[ -s "$1" ]]; }
