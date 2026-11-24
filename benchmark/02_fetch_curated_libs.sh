#!/usr/bin/env bash
# 02_fetch_curated_libs.sh — obtain the curated gold-standard TE libraries for the
# genomes that have one, into curated/<species>/<species>.curated.lib.fa.
#
#   rice/maize/arabidopsis -> EDTA repo (ship with EDTA)
#   dmel  -> RepBase drorep.ref — subscription (REPBASE_USER/PASS); proxy fallback if absent
#   human/danio/septoria/zymoseptoria -> Dfam clade libraries (fetch_dfam_lib;
#                                        9606 / 7955 / 215465 / 1047171)
#
# Any genome whose curated lib can't be obtained is left absent, so 05 scores it on the
# proxy branch instead. Re-entrant: an existing curated lib is kept.
#
# Also clones the EDTA repo (for lib-test.py used by 05). Re-entrant.
#
# Usage: ./02_fetch_curated_libs.sh
set -euo pipefail
source "$(dirname "$0")/env.sh"

# ── EDTA checkout (lib-test.py + curated rice/maize libraries) ─────────────
if [[ ! -f "${EDTA_REPO}/lib-test.py" ]]; then
    log "cloning EDTA -> ${EDTA_REPO}"
    rm -rf "$EDTA_REPO"
    git clone --depth 1 "$EDTA_GIT" "$EDTA_REPO" || die "git clone EDTA failed"
else
    log "EDTA checkout present (${EDTA_REPO})"
fi

# True if $1 is FASTA — first non-blank char is '>' (decompressing .gz transparently).
# Guards against an HTML error page or a JSON envelope being saved as a "library".
is_fasta() {
    [[ -s "$1" ]] || return 1
    local c; c="$(zcat -f "$1" 2>/dev/null | grep -m1 -v '^[[:space:]]*$' | cut -c1)"
    [[ "$c" == ">" ]]
}

place() {  # place <species_key> <source_file>
    local species="$1" src="$2"
    local dest="${CURATED_DIR}/${species}/${species}.curated.lib.fa"
    mkdir -p "$(dirname "$dest")"
    if is_fasta "$dest"; then log "skip ${species}: curated lib already in place"; return 0; fi
    [[ -s "$src" ]] || { log "WARN ${species}: curated source missing ($src)"; return 1; }
    cp "$src" "$dest"
    log "curated lib for ${species}: $(grep -c '^>' "$dest") sequences -> $dest"
}

# fetch_dfam_lib <species_key> <clade_taxid> — download a clade's curated TE library from
# Dfam (clade + descendants) and place it as the species' curated lib. Dfam's FASTA API
# omits RepeatMasker class tags, so dfam_to_rmlib.py re-attaches them (joining each family
# to its Dfam repeat_type/subtype) — otherwise lib-test.py can only score at the total
# level. On any failure the curated lib is left absent, so 05 falls back to the proxy
# branch. Override the URL per species with DFAM_URL_<species_key>=<raw or download=true url>.
fetch_dfam_lib() {
    local species="$1" clade="$2"
    local dest="${CURATED_DIR}/${species}/${species}.curated.lib.fa"
    if is_fasta "$dest"; then log "skip ${species}: curated lib already in place"; return 0; fi
    mkdir -p "$(dirname "$dest")"
    local raw="${CURATED_DIR}/${species}/dfam_${species}.raw.fasta"
    local ovar="DFAM_URL_${species}"
    local url="${!ovar:-https://www.dfam.org/api/families?format=fasta&clade=${clade}&clade_relatives=descendants&download=true}"

    # Download until the envelope FULLY decodes. The download=true response is one large
    # base64(gzip(FASTA)) blob; on a flaky link it truncates -> the base64 is cut off and
    # decoding dies with "Incorrect padding". A truncated file is still non-empty, so we
    # must NOT reuse it: validate with `dfam_to_rmlib.py --check` and re-download on failure
    # (bump attempts with DFAM_DL_RETRIES). After that, fall through to the proxy branch.
    local tries=0 max="${DFAM_DL_RETRIES:-5}"
    while ! python3 "${BENCH_ROOT}/dfam_to_rmlib.py" --check --fasta "$raw" 2>/dev/null; do
        if (( tries >= max )); then
            log "WARN ${species}: Dfam download never completed cleanly after ${max} tries"
            log "      (truncated on this link) -> proxy fallback. Retry on a stabler node,"
            log "      or raise DFAM_DL_RETRIES."
            rm -f "$raw"; return 0
        fi
        tries=$((tries + 1))
        log "downloading Dfam TE library for ${species} (clade ${clade}, descendants) [try ${tries}/${max}]"
        rm -f "$raw"
        if command -v wget >/dev/null 2>&1; then
            wget -q --tries=3 --timeout=120 -O "$raw" "$url" || true
        else
            curl -fsSL --retry 3 --retry-delay 5 --max-time 600 -o "$raw" "$url" || true
        fi
    done

    # Download validated complete; reclassify (dfam_to_rmlib.py attaches RepeatMasker class
    # tags from Dfam repeat_type/subtype) and place the library.
    if python3 "${BENCH_ROOT}/dfam_to_rmlib.py" --fasta "$raw" --out "$dest" --clade "$clade" \
            && is_fasta "$dest"; then
        log "Dfam ${species} curated lib: $(grep -c '^>' "$dest") sequences (classified)"
    else
        log "WARN ${species}: download OK but classification/rewrite failed -> proxy fallback."
        rm -f "$dest"
    fi
}

# rice + maize + arabidopsis curated libraries all ship with EDTA.
place oryza_sativa        "${EDTA_REPO}/database/rice7.0.0.liban"            || true
place zea_mays            "${EDTA_REPO}/database/maizeTE02052020"           || true
place arabidopsis_thaliana "${EDTA_REPO}/database/athrep.updated.nonredun.fasta" || true

# Drosophila — RepBase Drosophila reference (drorep.ref), a well-curated TE library.
# RepBase's 'protected' area needs a subscription login (HTTP Basic auth): set
#   export REPBASE_USER=<user> REPBASE_PASS=<pass>
# Override the source with DMEL_STD_URL=<raw FASTA url>, or drop a file at $DMEL_LOCAL.
# RepBase FASTA headers are 'name<TAB>class<TAB>species' (no RepeatMasker '#Class'),
# so repbase_to_rmlib.py rewrites them so lib-test.py can score per TE category.
REPBASE_DMEL_URL="${REPBASE_DMEL_URL:-https://www.girinst.org/server/RepBase/protected/RepBase31.06.fasta/drorep.ref}"
DMEL_LOCAL="${DMEL_LOCAL:-${CURATED_DIR}/drosophila_melanogaster/dmel_std6.28.fasta}"
dmel_url="${DMEL_STD_URL:-$REPBASE_DMEL_URL}"
dmel_dest="${CURATED_DIR}/drosophila_melanogaster/drosophila_melanogaster.curated.lib.fa"
dmel_raw="${CURATED_DIR}/drosophila_melanogaster/dmel.raw.ref"

if is_fasta "$dmel_dest"; then
    log "skip drosophila_melanogaster: curated lib already in place"
else
    mkdir -p "$(dirname "$dmel_dest")"
    # 1. Obtain the raw library: a local file wins; else download (with auth if set).
    if [[ -s "$DMEL_LOCAL" ]]; then
        cp "$DMEL_LOCAL" "$dmel_raw"
    elif ! is_fasta "$dmel_raw"; then
        log "downloading Drosophila RepBase library: ${dmel_url}"
        auth=()
        if [[ -n "${REPBASE_USER:-}" ]]; then
            if command -v wget >/dev/null 2>&1; then
                auth=(--user "$REPBASE_USER" --password "${REPBASE_PASS:-}")
            else
                auth=(-u "${REPBASE_USER}:${REPBASE_PASS:-}")
            fi
        fi
        # ${auth[@]+"${auth[@]}"} expands to nothing for an empty array — the plain
        # "${auth[@]}" trips `set -u` ("unbound variable") on bash < 4.4 (e.g. RHEL 7).
        if command -v wget >/dev/null 2>&1; then
            wget -q ${auth[@]+"${auth[@]}"} -O "$dmel_raw" "$dmel_url" || true
        else
            curl -fsSL ${auth[@]+"${auth[@]}"} -o "$dmel_raw" "$dmel_url" || true
        fi
    fi
    # 2. Validate + rewrite RepBase headers -> RepeatMasker name#Class.
    if is_fasta "$dmel_raw"; then
        if python3 "${BENCH_ROOT}/repbase_to_rmlib.py" --fasta "$dmel_raw" --out "$dmel_dest" \
                && is_fasta "$dmel_dest"; then
            log "Drosophila curated lib: $(grep -c '^>' "$dmel_dest") sequences (RepBase, classified)"
        else
            log "WARN drosophila_melanogaster: header conversion failed; using raw library"
            cp "$dmel_raw" "$dmel_dest"
        fi
    else
        rm -f "$dmel_dest"
        cat >&2 <<MSG
[02] NOTE: could not obtain a FASTA Drosophila library — Drosophila -> PROXY fallback.
     RepBase is subscription-protected. Provide credentials and re-run:
       export REPBASE_USER=<user> REPBASE_PASS=<pass>   ./02_fetch_curated_libs.sh
     Or override the source (must resolve to raw FASTA, not an HTML page):
       export DMEL_STD_URL=<raw url>                     ./02_fetch_curated_libs.sh
     Or place the FASTA at:
       ${DMEL_LOCAL}
MSG
    fi
fi

# Dfam curated TE libraries (clade + descendants):
#   homo_sapiens         9606    (~thousands of families)
#   danio_rerio          7955    (~1740 families — strong zebrafish coverage)
#   septoria_linicola    215465  (~43 families — modest but curated; flax-pasmo fungus,
#                                 INRAE MPMI 2022; a thin reference, so tools that call
#                                 more TEs than those 43 cover look falsely over-sensitive.
#                                 Its gold is ~94% LTR with ZERO LINE/SINE bp, so it cannot
#                                 test anything class-specific outside LTR.)
#   zymoseptoria_tritici 1047171 (324 families, ZymTri_2023_*, manually curated per Goubert
#                                 et al. 2022. ~36% carry a repeat_type — DNA/TIR ~23%,
#                                 LTR ~8%, LINE ~4%, RC/Helitron ~1% — so unlike Septoria
#                                 this genome has real per-class gold in every major class,
#                                 which is what the fungus route needs. The other ~64% are
#                                 type "Unknown" in Dfam and can only be scored at the
#                                 total level.)
#
# WARNING on the Zymoseptoria clade id: use the SPECIES taxid 1047171, not the strain
# taxid 336722 (Zymoseptoria tritici IPO323). The strain id is valid in NCBI taxonomy and
# the API accepts it happily — but Dfam returns total_count=0 for it, which would fetch an
# empty library and silently demote the genome to the proxy branch.
fetch_dfam_lib homo_sapiens         9606
fetch_dfam_lib danio_rerio          7955
fetch_dfam_lib septoria_linicola    215465
fetch_dfam_lib zymoseptoria_tritici 1047171

log "curated-library stage complete."
