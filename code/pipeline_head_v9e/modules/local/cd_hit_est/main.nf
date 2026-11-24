// Deduplicate/merge a TE consensus library with CD-HIT-EST at the 80-80-80 rule
// (params.cdhit_identity / params.cdhit_coverage).
// Stage 2 (LIBRARY_DEDUP): the de novo library plus an optional user library.
// Stage 3 (CURATION, as CD_HIT_EST_CURATED with ext.prefix '<id>.curated'): TEtrimmer's
// un-merged output, replacing TEtrimmer's own length-based merge (see DEDUP_REFINE).
process CD_HIT_EST {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(libraries, stageAs: 'lib_??.fa')

    output:
    tuple val(meta), path("*.dedup_library.fa")       , emit: library
    // representatives + cluster file + all inputs, for DEDUP_REFINE
    tuple val(meta), path("*.dedup_library.fa"), path("*.dedup_library.fa.clstr"),
          path("*.combined_library.fa")               , emit: clustered
    path "*.dedup_library.fa.clstr"                   , emit: clusters
    path 'versions.yml'                               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // cd-hit-est word size (-n) must match the identity threshold (-c). The
    // 0.80 default maps to -n 5; bump it for stricter thresholds. -p 1 writes each
    // member's alignment coordinates on its representative into the .clstr, which
    // DEDUP_REFINE uses to cut foreign-class segments out of chimeric consensi.
    // -G 0 -g 1 -b 500 make this the standard 80-80-80 family rule (Goubert et al.
    // 2022): identity over the local alignment, not the whole shorter sequence (-G 1,
    // the default, splits families whose consensi differ in indels or ends), each
    // sequence goes to its best cluster, and the alignment band allows long indels.
    def ident = params.cdhit_identity as double
    def word  = ident >= 0.90 ? 8 : ident >= 0.88 ? 7 : ident >= 0.85 ? 6 : 5
    def mem_mb = task.memory ? task.memory.toMega() : 4000
    """
    cat ${libraries} > ${prefix}.combined_library.fa

    cd-hit-est \\
        -i ${prefix}.combined_library.fa \\
        -o ${prefix}.dedup_library.fa \\
        -c ${params.cdhit_identity} \\
        -aS ${params.cdhit_coverage} \\
        -G 0 -g 1 -b 500 \\
        -n ${word} \\
        -d 0 \\
        -p 1 \\
        -T ${task.cpus} \\
        -M ${mem_mb}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$(cd-hit-est -h 2>&1 | head -1 | sed -E 's/.*version ([0-9.]+).*/\\1/' || echo unknown)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${libraries} > ${prefix}.combined_library.fa
    cp ${prefix}.combined_library.fa ${prefix}.dedup_library.fa
    touch ${prefix}.dedup_library.fa.clstr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: "stub"
    END_VERSIONS
    """
}
