
# $HOME on this HPC is read-only from inside the container; redirect HOME and
# caches to the (writable) task work dir for RepeatMasker/matplotlib/etc.
export HOME="$PWD"
export MPLCONFIGDIR="$PWD/.mplconfig"
export XDG_CACHE_HOME="$PWD/.cache"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

# TIR-Learner 3.0.1 as bundled with EDTA 2.3.0 does not run under the image's
# pandas 3: it crashes (x[0] on a labelled row -> KeyError: 0; and chromosome IDs
# with dots, e.g. NC_003070.9, never map back from their '<id>_split_NofM' chunks
# -> KeyError). EDTA then skips HelitronScanner too. Run a patched copy: EDTA finds
# TIR-Learner with `command -v TIR-Learner`, so a wrapper first on PATH suffices.
# The edits apply only where the buggy lines are present.
TL_WRAPPER=$(command -v TIR-Learner || true)
TL_SRC=$(grep -oE '/[^ ]*TIR-Learner[^ /]*\.py' "${TL_WRAPPER:-/dev/null}" 2>/dev/null | head -1 || true)
if [[ -n "${TL_SRC}" && -f "${TL_SRC}" ]]; then
    mkdir -p tirlearner/bin
    cp -r "$(dirname "${TL_SRC}")" tirlearner/TIR-Learner
    sed -i 's/^\(\s*\)family = x\[0\]$/\1family = x.iloc[0]/' tirlearner/TIR-Learner/bin/check_TIR_TSD.py
    sed -i 's/id_pattern = r"^(\\w+)_split_/id_pattern = r"^(\\S+?)_split_/' tirlearner/TIR-Learner/bin/run_GRF.py
    printf '#!/bin/bash\npython3 %s "$@"\n' "$PWD/tirlearner/TIR-Learner/$(basename "${TL_SRC}")" > tirlearner/bin/TIR-Learner
    chmod +x tirlearner/bin/TIR-Learner
    export PATH="$PWD/tirlearner/bin:$PATH"
else
    echo "WARN [__SP__]: TIR-Learner source not found; running it unpatched" >&2
fi

# EDTA's redundancy filter cleanup_nested.pl (final library; stage-1 TIR/Helitron
# library) runs threads with a race: two near-identical copies blasted at the same
# time each delete the other, so multi-copy (mostly LTR) families vanish at random
# (Zymoseptoria tritici: a 13 kb Gypsy and a 6 kb Copia, ~1 Mb of LTR). Run EDTA
# from a work-dir copy (database linked, not copied) whose cleanup_nested.pl applies
# each query's edits under a lock and skips a query that another thread already
# changed or removed. EDTA's scripts find their helpers via FindBin, so every stage
# uses the copy. (the bioconda EDTA.pl is a wrapper: perl /usr/local/share/EDTA/EDTA.pl)
EDTA_BIN=$(command -v EDTA.pl)
EDTA_DIR=$(grep -oE '/[^ ]*/EDTA/' "${EDTA_BIN}" 2>/dev/null | head -1 || true)
[[ -d "${EDTA_DIR}database" ]] || EDTA_DIR="$(dirname "$(readlink -f "${EDTA_BIN}")")/"
EDTA_CMD=(EDTA.pl)
if [[ -f "${EDTA_DIR}EDTA.pl" && -f "${EDTA_DIR}bin/cleanup_nested.pl" ]]; then
    rm -rf edta_src && mkdir edta_src
    for f in "${EDTA_DIR}"*; do
        case "$(basename "$f")" in
            database|test) ln -s "$f" edta_src/ ;;
            *)             cp -r "$f" edta_src/ ;;
        esac
    done
    chmod -R u+w edta_src
    sed -i '/# removing the regions from the subject that are inserted into the query/a lock(%seq); next if !exists $seq{$id} or $touched_seq{$id};' edta_src/bin/cleanup_nested.pl
    if grep -q '^lock(%seq); next if' edta_src/bin/cleanup_nested.pl; then
        EDTA_CMD=(perl "$PWD/edta_src/EDTA.pl")
    else
        echo "WARN [__SP__]: cleanup_nested.pl not patched (unexpected EDTA version); running EDTA unpatched" >&2
    fi
else
    echo "WARN [__SP__]: EDTA scripts not found under '${EDTA_DIR}'; running EDTA unpatched" >&2
fi

# --force 1: continue when a sub-module legitimately finds no TEs instead of
# aborting. EDTA then copies its bundled RICE library for that module, which must
# never reach the output: that is checked right after the run.
"${EDTA_CMD[@]}" \
    --genome __SP__.fna \
    --species others \
    --step all \
    --anno 0 \
    --force 1 \
    --threads __CPUS__ \
    --overwrite 0

# Fail if EDTA substituted its rice library for any module (--force 1 fallback):
# the raw file is then a verbatim copy of <EDTA>/database/rice7.0.0.liban.<TYPE>.
EDTA_DB="${EDTA_DIR}database"
[[ -d "${EDTA_DB}" ]] || echo "WARN [__SP__]: EDTA database dir not found; rice-fallback check skipped" >&2
RICE_USED=""
for t in LTR:LTR.raw.fa LINE:LINE.raw.fa SINE:SINE.raw.fa TIR:TIR.intact.raw.fa Helitron:Helitron.intact.raw.fa; do
    RAW=$(find . -maxdepth 2 -path '*.EDTA.raw/*' -name "*.mod.${t#*:}" -print -quit)
    RICE="${EDTA_DB}/rice7.0.0.liban.${t%%:*}"
    if [[ -n "${RAW}" && -s "${RICE}" ]] && cmp -s "${RAW}" "${RICE}"; then
        RICE_USED="${RICE_USED} ${t%%:*}"
    fi
done
if [[ -n "${RICE_USED}" ]]; then
    echo "ERROR [__SP__]: EDTA found no${RICE_USED} candidates and substituted its bundled RICE library (--force 1)." >&2
    echo "  Check the EDTA log (.command.log) for the failing module; the library would contain rice TEs." >&2
    exit 1
fi

# Final non-redundant, classified library: <genome>.mod.EDTA.TElib.fa
LIB=$(find . -name '*.mod.EDTA.TElib.fa' -print -quit 2>/dev/null || true)
if [[ -n "${LIB}" && -s "${LIB}" ]]; then
    cp "${LIB}" __SP__.edta_library.fa
else
    echo "WARN [__SP__]: EDTA produced no TElib.fa" >&2
    touch __SP__.edta_library.fa
fi

EDTA_VER=$(ls "${CONDA_PREFIX:-}"/conda-meta/edta-*.json /usr/local/conda-meta/edta-*.json 2>/dev/null | sed -E 's@.*/edta-([0-9][^-]*)-.*@\1@' | head -1 || true)
[[ -z "${EDTA_VER}" ]] && EDTA_VER=2.3.0
cat <<-END_VERSIONS > versions.yml
"EDTA":
    edta: ${EDTA_VER}
END_VERSIONS
