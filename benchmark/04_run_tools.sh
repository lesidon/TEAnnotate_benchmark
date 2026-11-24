#!/usr/bin/env bash
# 04_run_tools.sh — build each INDIVIDUAL de-novo library, one tool per env, into
#   libraries/<species>/<method>/<method>_library.fa
# for method in {edta, earlgrey, hite, heliano}. (The full pipeline is stage 03.)
#
# Each invocation mirrors the corresponding pipeline module (same CLI, same output
# glob, same HOME/MPLCONFIGDIR/writable-tmpfs/famdb workarounds). Re-entrant: a
# non-empty output library is skipped.
#
# Usage:
#   ./04_run_tools.sh                          # every genome, every method
#   ./04_run_tools.sh oryza_sativa             # one genome, every method
#   ./04_run_tools.sh oryza_sativa hite        # one genome, one method
#   METHODS="edta hite" ./04_run_tools.sh      # restrict the method set globally
set -euo pipefail
source "$HOME/tmp/teannotate_benchmark/whole_pipeline_benchmark/env.sh"

ONLY_SP="${1:-}"
ONLY_M="${2:-}"
METHODS="${METHODS:-edta earlgrey hite heliano}"
[[ -n "$ONLY_M" ]] && METHODS="$ONLY_M"

# ── conda-env bootstrap (built once, then cached) ──────────────────────────
ensure_earlgrey_env() {
    conda_run "$EARLGREY_ENV" true 2>/dev/null && return 0
    [[ -n "$NF_PIPELINE" && -f "${NF_PIPELINE}/assets/environments/earl-grey.yml" ]] \
        || die "earlgrey env '${EARLGREY_ENV}' missing and no earl-grey.yml under NF_PIPELINE to build it"
    log "creating conda env ${EARLGREY_ENV} from earl-grey.yml (one-time)"
    conda env create -n "$EARLGREY_ENV" -f "${NF_PIPELINE}/assets/environments/earl-grey.yml"
}
ensure_heliano_env() {
    conda_run "$HELIANO_ENV" true 2>/dev/null && return 0
    log "creating conda env ${HELIANO_ENV} (bioconda::heliano) (one-time)"
    conda create -y -n "$HELIANO_ENV" -c conda-forge -c bioconda heliano
}

# HOME/cache redirection for tools that write dotfiles (read-only $HOME on Zhores).
tool_env_prefix() {  # emit exports to run before a tool, into $1 = workdir
    local wd="$1"
    export HOME="$wd"
    export MPLCONFIGDIR="${wd}/.mplconfig"
    export XDG_CACHE_HOME="${wd}/.cache"
    mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
}

# ── per-method runners ─────────────────────────────────────────────────────
run_edta() {  # run_edta <genome> <workdir> <dest>
    local genome="$1" wd="$2" dest="$3"
    ( cd "$wd"; tool_env_prefix "$wd"
      container_run "$EDTA_SIF" EDTA.pl \
          --genome "$genome" --species others --step all --anno 0 --force 1 \
          --threads "$THREADS" --overwrite 0 )
    local lib; lib="$(find "$wd" -name '*.mod.EDTA.TElib.fa' -print -quit 2>/dev/null || true)"
    [[ -n "$lib" && -s "$lib" ]] && cp "$lib" "$dest"
}

run_earlgrey_container() {  # run_earlgrey_container <genome> <workdir> <dest> <species>
    local genome="$1" wd="$2" dest="$3" species="$4"
    local out="${wd}/earlgrey_out"; mkdir -p "$out"
    # Write an in-container driver (a real file in the bound work dir — avoids shell-quoting
    # hell). Host-side vars are expanded now; \$… stays literal for in-container evaluation.
    cat > "${wd}/_eg_run.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
# EarlGrey's tools may live in a conda env the image doesn't auto-activate under
# 'singularity exec' (no Docker ENTRYPOINT). Put them on PATH if not already.
if ! command -v earlGreyLibConstruct >/dev/null 2>&1; then
    for p in /opt/conda/etc/profile.d/conda.sh /usr/local/etc/profile.d/conda.sh; do
        [[ -f "\$p" ]] && source "\$p" || true
    done
    conda activate earlgrey 2>/dev/null || conda activate earlGrey 2>/dev/null || true
fi
command -v earlGreyLibConstruct >/dev/null 2>&1 || { echo "earlGreyLibConstruct not on PATH in \$EARLGREY_SIF" >&2; exit 3; }
# earlGreyLibConstruct gates de novo behind TWO interactive Dfam checks; both must be made
# non-interactive (--writable-tmpfs makes the image FS writable for these edits):
#  (1) dfamCheck() exits 2 when the famdb dir has EXACTLY 2 files (partition 0 + rmlib.config).
#      Touch its .earlgrey.config.complete sentinel so the count is >=3 and the check passes.
for f in \$(find /opt /usr/local /usr/share -type d -name famdb 2>/dev/null); do
    touch "\$f/.earlgrey.config.complete" 2>/dev/null || true
done
#  (2) A separate block (in Checks()) prompts ("Please answer yes or no." forever on
#      non-interactive stdin) when RepeatMasker's LEGACY Dfam.h5 is an unconfigured
#      "Placeholder". De novo does not use it (RepeatMasker here is famdb-configured), so
#      point Dfam.h5 at the real partition-0 h5 via a SYMLINK so 'grep -q Placeholder' is
#      false and the prompt is skipped. A symlink (not a copy) avoids the tiny --writable-tmpfs
#      overlay running out of space on the ~77 MB partition file.
for lib in \$(find /opt /usr/local /usr/share -type d -name Libraries -path '*RepeatMasker*' 2>/dev/null); do
    dfamh5="\$lib/Dfam.h5"
    [ -f "\$dfamh5" ] && grep -q Placeholder "\$dfamh5" 2>/dev/null || continue
    part0=\$(find "\$lib/famdb" -maxdepth 1 -name 'dfam*_full.0.h5' 2>/dev/null | head -1)
    [ -n "\$part0" ] && ln -sf "\$part0" "\$dfamh5" 2>/dev/null || true
    echo "Dfam.h5 gate (\$lib): \$(grep -q Placeholder "\$dfamh5" 2>/dev/null && echo STILL-PLACEHOLDER || echo cleared)" >&2
done
earlGreyLibConstruct -g "$genome" -s "$species" -o "$out" -t "$THREADS"
EOF
    chmod +x "${wd}/_eg_run.sh"
    ( cd "$wd"; tool_env_prefix "$wd"
      container_run "$EARLGREY_SIF" bash "${wd}/_eg_run.sh" )
    local lib
    lib="$(find "$out" -path '*summaryFiles*' \( -name '*-families.fa' -o -name '*-families.fa.strained' \) -print -quit 2>/dev/null || true)"
    [[ -z "$lib" ]] && lib="$(find "$out" \( -name '*-families.fa.strained' -o -name '*-families.fa' \) -print -quit 2>/dev/null || true)"
    [[ -n "$lib" && -s "$lib" ]] && cp "$lib" "$dest"
}

run_earlgrey() {  # run_earlgrey <genome> <workdir> <dest> <species>
    local genome="$1" wd="$2" dest="$3" species="$4"
    # Prefer the prebuilt container when EARLGREY_SIF is set (skips the slow conda solve).
    if [[ -n "${EARLGREY_SIF:-}" ]]; then
        run_earlgrey_container "$genome" "$wd" "$dest" "$species"; return
    fi
    ensure_earlgrey_env
    ( cd "$wd"; tool_env_prefix "$wd"
      # Resolve earlGreyLibConstruct + create the famdb sentinel so EarlGrey skips
      # the multi-GB Dfam download (de novo needs no clade library). Mirrors the module.
      EG_LIB="$(conda_run "$EARLGREY_ENV" bash -lc 'command -v earlGreyLibConstruct || find "${CONDA_PREFIX:-/opt}" -name earlGreyLibConstruct -type f 2>/dev/null | head -1')"
      [[ -n "$EG_LIB" ]] || die "earlGreyLibConstruct not found in ${EARLGREY_ENV}"
      RM_BIN="$(conda_run "$EARLGREY_ENV" bash -lc 'command -v RepeatMasker' || true)"
      if [[ -n "$RM_BIN" ]]; then
          FAMDB="$(dirname "$RM_BIN")/../share/RepeatMasker/Libraries/famdb"
          [[ -d "$FAMDB" ]] && touch "${FAMDB}/.earlgrey.config.complete" || true
      fi
      conda_run "$EARLGREY_ENV" bash "$EG_LIB" -g "$genome" -s "$species" -o "${wd}/earlgrey_out" -t "$THREADS" )
    local lib
    lib="$(find "${wd}/earlgrey_out" -path '*summaryFiles*' \( -name '*-families.fa' -o -name '*-families.fa.strained' \) -print -quit 2>/dev/null || true)"
    [[ -z "$lib" ]] && lib="$(find "${wd}/earlgrey_out" \( -name '*-families.fa.strained' -o -name '*-families.fa' \) -print -quit 2>/dev/null || true)"
    [[ -n "$lib" && -s "$lib" ]] && cp "$lib" "$dest"
}

run_hite() {  # run_hite <genome> <workdir> <dest>
    local genome="$1" wd="$2" dest="$3"
    ( cd "$wd"; tool_env_prefix "$wd"
      container_run "$HITE_SIF" python /HiTE/main.py \
          --genome "$genome" --out_dir "${wd}/hite_out" --thread "$THREADS" --annotate 0 )
    local lib; lib="$(find "${wd}/hite_out" -name 'confident_TE.cons.fa' -print -quit 2>/dev/null || true)"
    [[ -n "$lib" && -s "$lib" ]] && cp "$lib" "$dest"
}

run_heliano() {  # run_heliano <genome> <workdir> <dest>
    local genome="$1" wd="$2" dest="$3"
    ensure_heliano_env
    ( cd "$wd"; tool_env_prefix "$wd"
      conda_run "$HELIANO_ENV" heliano -g "$genome" -o "${wd}/heliano_out" -n "$THREADS" )
    local rep; rep="$(find "${wd}/heliano_out" -name 'RC.representative.fa' -print -quit 2>/dev/null || true)"
    if [[ -n "$rep" && -s "$rep" ]]; then
        # Relabel headers name#RC/Helitron (RepeatMasker convention), as the module does.
        awk '/^>/{ split($0,a," "); n=a[1]; if (n !~ /#/) n=n"#RC/Helitron"; print n; next } {print}' \
            "$rep" > "$dest"
    fi
}

run_method() {  # run_method <kingdom> <species> <method>
    local kingdom="$1" species="$2" method="$3"
    local genome; genome="$(readlink -f "$(genome_fasta "$kingdom" "$species")")"
    local dest="${LIB_DIR}/${species}/${method}/${method}_library.fa"
    local wd="${WORK_DIR}/tools/${species}/${method}"

    nonempty "$genome" || { log "WARN ${species}: genome missing — run 01 first"; return 0; }
    if nonempty "$dest"; then log "skip ${species}/${method}: library already built"; return 0; fi
    mkdir -p "$(dirname "$dest")" "$wd"

    log "build ${species}/${method} ..."
    case "$method" in
        edta)     run_edta     "$genome" "$wd" "$dest" ;;
        earlgrey) run_earlgrey "$genome" "$wd" "$dest" "$species" ;;
        hite)     run_hite     "$genome" "$wd" "$dest" ;;
        heliano)  run_heliano  "$genome" "$wd" "$dest" ;;
        *) log "WARN: unknown method '${method}'"; return 0 ;;
    esac

    if nonempty "$dest"; then
        log "done ${species}/${method}: $(grep -c '^>' "$dest") sequences"
    else
        log "WARN ${species}/${method}: produced no library (left absent; excluded from scoring)"
    fi
}

read_manifest "$GENOMES_TSV" | while IFS=$'\t' read -r kingdom species display accession organism scoring curated; do
    [[ -n "$ONLY_SP" && "$species" != "$ONLY_SP" ]] && continue
    for m in $METHODS; do run_method "$kingdom" "$species" "$m"; done
done

log "individual-tool stage complete."
