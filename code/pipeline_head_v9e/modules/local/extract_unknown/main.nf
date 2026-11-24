// Stage 3 (re-classification) — subset the curated library to its "Unknown"
// fraction. Only these consensi are sent to the DL ensemble (TEclass2/DeepTE/
// CREATE); already-classified consensi keep their curated labels. See
// bin/extract_unknown.py for the Unknown-matching rules.
process EXTRACT_UNKNOWN {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(library)

    output:
    tuple val(meta), path("${meta.id}.unknown.fa"), emit: unknown
    path 'versions.yml'                           , emit: versions

    script:
    """
    extract_unknown.py --in ${library} --out ${meta.id}.unknown.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${library} ${meta.id}.unknown.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
