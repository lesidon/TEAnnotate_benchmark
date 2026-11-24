process INPUT_CHECK {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    // A gzipped genome is decompressed here once (<name> without .gz): every
    // downstream tool (RepeatMasker, EDTA, awk-based plots, ...) needs plain FASTA.
    // A plain genome passes through unchanged.
    output:
    tuple val(meta), path("${fasta.name.endsWith('.gz') ? fasta.baseName : fasta.name}", includeInputs: true), emit: genome
    path  'versions.yml'        , emit: versions

    script:
    def gz = fasta.name.endsWith('.gz')
    def plain = gz ? fasta.baseName : fasta.name
    """
    # ── 1. Non-empty check ────────────────────────────────────────────────────
    [[ -s "${fasta}" ]] || {
        echo "ERROR [${meta.id}]: genome file '${fasta}' is empty." >&2
        exit 1
    }

    # ── 2. FASTA format: first non-blank line must start with '>' ─────────────
    # grep -m 1 exits after the first match, avoiding the grep|head SIGPIPE
    # that triggers exit 141 under set -o pipefail.
    # zcat still receives SIGPIPE when grep -m 1 exits; pipefail is suspended
    # for that one pipeline only.
    set +o pipefail
    if ${gz}; then
        FIRST=\$(zcat "${fasta}" | grep -m 1 -v '^[[:space:]]*\$')
    else
        FIRST=\$(grep -m 1 -v '^[[:space:]]*\$' "${fasta}")
    fi
    set -o pipefail

    [[ "\${FIRST:0:1}" == ">" ]] || {
        echo "ERROR [${meta.id}]: '${fasta}' does not appear to be FASTA." >&2
        echo "  First non-blank line: \${FIRST}" >&2
        exit 1
    }

    # ── 3. At least one sequence header ───────────────────────────────────────
    if ${gz}; then
        SEQ_COUNT=\$(zcat "${fasta}" | grep -c '^>' || true)
    else
        SEQ_COUNT=\$(grep -c '^>' "${fasta}" || true)
    fi

    [[ \$SEQ_COUNT -ge 1 ]] || {
        echo "ERROR [${meta.id}]: no FASTA headers found in '${fasta}'." >&2
        exit 1
    }

    # ── 4. Plain FASTA for every downstream step ──────────────────────────────
    if ${gz}; then
        gzip -t "${fasta}" || { echo "ERROR [${meta.id}]: '${fasta}' is a corrupt gzip stream." >&2; exit 1; }
        gzip -cd "${fasta}" | tr -d '\\r' > "${plain}"
    elif grep -q \$'\\r' "${fasta}"; then
        echo "ERROR [${meta.id}]: '${fasta}' has Windows (CRLF) line endings, which break" >&2
        echo "  sequence-length counting in several tools. Convert it first, e.g. dos2unix." >&2
        exit 1
    fi

    # ── 5. Unique sequence IDs (EDTA, RepeatMasker and the evidence gate need them) ─
    set +o pipefail          # head may close the pipe early
    DUPS=\$(awk '/^>/ { print substr(\$1, 2) }' "${plain}" | sort | uniq -d | head -5 | tr '\\n' ' ')
    set -o pipefail
    [[ -z "\${DUPS}" ]] || {
        echo "ERROR [${meta.id}]: duplicated sequence IDs in '${fasta}': \${DUPS}" >&2
        exit 1
    }

    echo "INFO [${meta.id}]: FASTA OK — \${SEQ_COUNT} sequence(s)" >&2

    set +o pipefail
    BASH_VER=\$(bash --version | head -1 | sed 's/GNU bash, version //;s/ .*//')
    set -o pipefail
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$BASH_VER
    END_VERSIONS
    """

    stub:
    def plain = fasta.name.endsWith('.gz') ? fasta.baseName : fasta.name
    """
    if [[ "${fasta.name}" == *.gz ]]; then gzip -cd "${fasta}" > "${plain}"; fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub"
    END_VERSIONS
    """
}
