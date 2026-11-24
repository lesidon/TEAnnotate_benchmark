// The final TE library of a sample, i.e. the one RepeatMasker masks with, published
// under one name whichever branch produced it: the ensemble-relabelled library, the
// curated library when re-classification is skipped or every classifier failed, or
// the deduplicated library under --skip_tetrimmer.
process FINAL_LIBRARY {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(library, stageAs: 'input_library.fa')

    output:
    tuple val(meta), path("${meta.id}.te_library.fa"), emit: library
    path 'versions.yml'                              , emit: versions

    script:
    """
    cp -L input_library.fa ${meta.id}.te_library.fa
    N=\$(grep -c '^>' ${meta.id}.te_library.fa || true)
    U=\$(grep -ciE '^>[^ ]*#unknown' ${meta.id}.te_library.fa || true)
    echo "INFO [${meta.id}]: final TE library: \${N} consensi, \${U} still Unknown" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(cp --version 2>&1 | head -1 | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    """
    cp -L input_library.fa ${meta.id}.te_library.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: "stub"
    END_VERSIONS
    """
}
