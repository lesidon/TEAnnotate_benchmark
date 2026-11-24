process INPUT_CHECK {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), emit: genome
    path  'versions.yml'        , emit: versions

    script:
    def gz = fasta.name.endsWith('.gz')
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
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub"
    END_VERSIONS
    """
}
