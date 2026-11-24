#!/usr/bin/env bash
# After a -stub run: every sample has a final library with sequences in it, the run has
# a version and an asset manifest, no real asset was downloaded into the caches (stub
# outputs live in 'stub-' cache directories), and, given the user library, each of its
# records reached every sample's final library (as <id> or, renamed, <id>_user*).
#   tests/stub/check_outputs.sh <outdir> <samplesheet> [<launch dir>] [<user library>]
set -euo pipefail
OUT=$1; SHEET=$2; LAUNCH=${3:-.}; USER_LIB=${4:-}
fail=0

for id in $(tail -n +2 "$SHEET" | cut -d, -f1); do
    lib="$OUT/$id/05_classification/$id.te_library.fa"
    n=$(grep -c '^>' "$lib" 2>/dev/null || true)
    if [[ "${n:-0}" -gt 0 ]]; then echo "ok    $lib ($n sequences)"; else echo "EMPTY OR MISSING $lib"; fail=1; continue; fi
    if [[ -n "$USER_LIB" ]]; then
        for u in $(grep '^>' "$USER_LIB" | sed 's/^>//; s/[# ].*//'); do
            if grep -qE "^>${u}(_user[0-9]*)?(#|\$| )" "$lib"; then echo "ok    user record $u in $id"
            else echo "USER RECORD $u MISSING from $lib"; fail=1; fi
        done
    fi
done

for f in software_versions.yml asset_manifest.tsv; do
    v="$OUT/pipeline_info/$f"
    if [[ -s "$v" ]]; then echo "ok    $v"; else echo "MISSING $v"; fail=1; fi
done

while IFS= read -r m; do
    case "$m" in
        */stub-*) echo "ok    $m" ;;
        *)        echo "REAL ASSET IN A STUB RUN: $m"; fail=1 ;;
    esac
done < <(find "$LAUNCH" -maxdepth 6 -path '*_cache/*' -name source.tsv 2>/dev/null)

exit $fail
