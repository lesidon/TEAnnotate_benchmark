// Dedicated Helitron detection with HELIANO (Li et al. 2024).
// HELIANO is a fast, accurate Helitron-like element (HLE) detector — better than
// EDTA's older bundled HelitronScanner — so it complements EDTA in the plant
// route. Conda-only (bioconda::heliano; no container image published). The
// representative Helitron sequences are relabelled #RC/Helitron for the library.
process HELIANO {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}.heliano_library.fa"), emit: library
    path 'versions.yml'                                   , emit: versions

    script:
    """
    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    mkdir -p "\$MPLCONFIGDIR"

    heliano \\
        -g ${genome} \\
        -o heliano_out \\
        -n ${task.cpus}

    # Representative HLE sequences → library, labelled name#RC/Helitron (RepeatMasker
    # convention), dropping any description so the header is a clean TE library entry.
    REP=\$(find heliano_out -name 'RC.representative.fa' -print -quit 2>/dev/null || true)
    if [[ -n "\${REP}" && -s "\${REP}" ]]; then
        awk '/^>/{ split(\$0,a," "); n=a[1]; if (n !~ /#/) n=n"#RC/Helitron"; print n; next } {print}' \\
            "\${REP}" > ${meta.id}.heliano_library.fa
    else
        echo "WARN [${meta.id}]: HELIANO produced no RC.representative.fa" >&2
        touch ${meta.id}.heliano_library.fa
    fi

    HEL_VER=\$(ls "\${CONDA_PREFIX:-}"/conda-meta/heliano-*.json 2>/dev/null | sed -E 's@.*/heliano-([0-9][^-]*)-.*@\\1@' | head -1 || true)
    [[ -z "\${HEL_VER}" ]] && HEL_VER=1.3.1
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        heliano: \${HEL_VER}
    END_VERSIONS
    """

    stub:
    """
    cat > ${meta.id}.heliano_library.fa <<-EOF
    >HLE_stub#RC/Helitron
    GGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCT
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        heliano: "stub"
    END_VERSIONS
    """
}
