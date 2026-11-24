#!/usr/bin/env bash
# 05_mask_and_score.sh — mask each genome with every method's library (one
# standardized RepeatMasker run), then score by branch:
#
#   GOLD  (rice/maize/dmel, curated lib present): build a reference .out from the
#         curated library, then lib-test.py per method -> results/gold/<sp>.<method>.metrics.tsv
#   PROXY (everything else): proxy_metrics.py over all method .out files ->
#         results/proxy/<sp>.{coverage,concordance}.tsv  (+ LAI context)
#
# Usage:
#   ./05_mask_and_score.sh                 # all genomes
#   ./05_mask_and_score.sh oryza_sativa    # one species_key
#
# Re-entrant: existing non-empty .out / metrics files are reused.
set -euo pipefail
source "$(dirname "$0")/env.sh"

ONLY="${1:-}"
# RepeatMasker -pa launches N jobs of ~4 cores each; divide the thread budget.
RM_PA=$(( THREADS / 4 )); (( RM_PA < 1 )) && RM_PA=1 || true

# repeatmask <genome> <lib.fa> <out.out> — one standardized RepeatMasker run.
repeatmask() {
    local genome="$1" lib="$2" outfile="$3"
    nonempty "$outfile" && { log "  reuse mask -> $(basename "$outfile")"; return 0; }
    local wd; wd="$(mktemp -d "${WORK_DIR}/rm.XXXXXX")"
    # shellcheck disable=SC2086
    container_run "$RM_SIF" RepeatMasker -pa "$RM_PA" $RM_OPTS -lib "$lib" -dir "$wd" "$genome" \
        || { log "  RepeatMasker failed for $(basename "$lib")"; rm -rf "$wd"; return 1; }
    local produced; produced="$(find "$wd" -name '*.out' -print -quit 2>/dev/null || true)"
    [[ -n "$produced" && -s "$produced" ]] || { log "  no .out produced for $(basename "$lib")"; rm -rf "$wd"; return 1; }
    mkdir -p "$(dirname "$outfile")"; cp "$produced" "$outfile"; rm -rf "$wd"
}

# compute_lai <genome> <outfile.txt> — whole-genome LAI via LTR_retriever (EDTA img).
compute_lai() {
    local genome="$1" outfile="$2"
    [[ "$RUN_LAI" == "1" ]] || { echo "NA" > "$outfile"; return 0; }
    nonempty "$outfile" && { log "  reuse LAI"; return 0; }
    local wd; wd="$(mktemp -d "${WORK_DIR}/lai.XXXXXX")"
    ln -sf "$genome" "${wd}/genome.fa"
    ( cd "$wd"
      container_run "$EDTA_SIF" bash -c '
        set -e; cd "'"$wd"'"
        gt suffixerator -db genome.fa -indexname genome.fa -tis -suf -lcp -des -ssp -sds -dna
        gt ltrharvest -index genome.fa -similar 90 -vic 10 -seed 20 -seqids yes \
           -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 \
           > genome.fa.harvest.scn
        LTR_retriever -genome genome.fa -inharvest genome.fa.harvest.scn -threads '"$THREADS"'
      ' ) || { log "  LAI computation failed (non-fatal)"; echo "NA" > "$outfile"; rm -rf "$wd"; return 0; }
    local lai; lai="$(grep -w 'whole_genome' "${wd}"/genome.fa*.LAI 2>/dev/null | awk '{print $NF}' | head -1 || true)"
    [[ -n "$lai" ]] && echo "$lai" > "$outfile" || echo "NA" > "$outfile"
    rm -rf "$wd"
}

score_one() {
    local kingdom="$1" species="$2" scoring="$3"
    local genome; genome="$(readlink -f "$(genome_fasta "$kingdom" "$species")")"
    nonempty "$genome" || { log "WARN ${species}: genome missing — run 01"; return 0; }

    local maskdir="${MASK_DIR}/${species}"; mkdir -p "$maskdir"
    log "== ${species} (${scoring}) =="

    # Mask with every method library that exists.
    local -a methods=() outs=()
    local libf method mout
    for libf in "${LIB_DIR}/${species}"/*/*_library.fa; do
        [[ -s "$libf" ]] || continue
        method="$(basename "$(dirname "$libf")")"
        mout="${maskdir}/${species}.${method}.out"
        log "  mask ${method}"
        repeatmask "$genome" "$libf" "$mout" || continue
        methods+=("$method"); outs+=("$mout")
    done
    [[ ${#methods[@]} -gt 0 ]] || { log "  no method libraries to score (run 03/04)"; return 0; }

    # Does this genome have a curated library present? (gold only if so.)
    local curated="${CURATED_DIR}/${species}/${species}.curated.lib.fa"
    if [[ "$scoring" == "gold" && -s "$curated" ]]; then
        # ---- GOLD: lib-test.py against a reference .out built from the curated lib.
        local ref="${maskdir}/${species}.reference.out"
        log "  mask reference (curated)"
        repeatmask "$genome" "$curated" "$ref" || { log "  reference mask failed; skipping gold"; return 0; }
        mkdir -p "${RESULT_DIR}/gold"
        local i
        for i in "${!methods[@]}"; do
            method="${methods[$i]}"; mout="${outs[$i]}"
            local report="${mout}.repeatmasker.lib.report"
            local metrics="${RESULT_DIR}/gold/${species}.${method}.metrics.tsv"
            nonempty "$metrics" && { log "  reuse gold metrics ${method}"; continue; }
            log "  lib-test ${method}"
            ( cd "$maskdir" && container_run "$EDTA_SIF" python "${EDTA_REPO}/lib-test.py" \
                --genome "$genome" --reference "$ref" --test "$mout" ) \
                || { log "  lib-test.py failed for ${method}"; continue; }
            # The report is <test>.repeatmasker.lib.report; fall back to a search if
            # lib-test.py wrote it relative to CWD instead of next to the test path.
            if [[ ! -s "$report" ]]; then
                report="$(find "$maskdir" -name "*.${method}.out.repeatmasker.lib.report" -print -quit 2>/dev/null || true)"
                [[ -z "$report" ]] && report="$(find "$maskdir" -name '*.repeatmasker.lib.report' -newer "$mout" -print -quit 2>/dev/null || true)"
            fi
            [[ -n "$report" && -s "$report" ]] || { log "  no report for ${method}"; continue; }
            python "${BENCH_ROOT}/parse_libtest.py" \
                --report "$report" --species "$species" --kingdom "$kingdom" \
                --method "$method" --out "$metrics"
        done
    else
        # ---- PROXY: coverage + pairwise concordance (+ LAI context).
        [[ "$scoring" == "gold" ]] && log "  NOTE: curated lib absent -> proxy fallback for ${species}"
        local lai="${maskdir}/${species}.lai.txt"
        compute_lai "$genome" "$lai"
        mkdir -p "${RESULT_DIR}/proxy"
        local -a outargs=()
        for i in "${!methods[@]}"; do outargs+=("${methods[$i]}=${outs[$i]}"); done
        python "${BENCH_ROOT}/proxy_metrics.py" \
            --genome "$genome" --species "$species" --kingdom "$kingdom" \
            --out "${outargs[@]}" --lai-file "$lai" \
            --outdir "${RESULT_DIR}/proxy"
    fi
    log "  ${species} scored."
}

read_manifest "$GENOMES_TSV" | while IFS=$'\t' read -r kingdom species display accession organism scoring curated; do
    [[ -n "$ONLY" && "$species" != "$ONLY" ]] && continue
    score_one "$kingdom" "$species" "$scoring"
done

log "scoring stage complete. Aggregate with:  python aggregate.py --outdir results"
