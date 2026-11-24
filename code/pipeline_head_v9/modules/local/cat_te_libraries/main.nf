// Concatenate the per-tool TE consensus libraries into one raw combined library, and
// record which tool each sequence came from (<id>.provenance.tsv: name<TAB>tool) for
// the evidence gate's cross-tool support. Deduplication happens downstream
// (CD-HIT-EST at the 80-80-80 rule). Inputs keep their names (<id>.<tool>_library.fa),
// which carry the tool.
process CAT_TE_LIBRARIES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(libraries)

    output:
    tuple val(meta), path("${meta.id}.raw_library.fa"), emit: library
    tuple val(meta), path("${meta.id}.provenance.tsv"), emit: provenance
    path 'versions.yml'                               , emit: versions

    script:
    """
    : > "${meta.id}.raw_library.fa"
    : > "${meta.id}.provenance.tsv"
    for f in ${libraries}; do
        tool=\$(basename "\$f" .fa)
        tool=\${tool#${meta.id}.}
        tool=\${tool%_library}
        cat "\$f" >> "${meta.id}.raw_library.fa"
        awk -v t="\$tool" '/^>/ { n = substr(\$1, 2); sub(/#.*/, "", n); print n "\\t" t }' "\$f" \\
            >> "${meta.id}.provenance.tsv"
    done

    set +o pipefail
    CAT_VER=\$(cat --version 2>&1 | head -1 | sed 's/cat (GNU coreutils) //')
    set -o pipefail

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$CAT_VER
    END_VERSIONS
    """

    stub:
    """
    cat ${libraries} > "${meta.id}.raw_library.fa"
    touch "${meta.id}.provenance.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: "stub"
    END_VERSIONS
    """
}
