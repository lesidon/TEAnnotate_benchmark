// Deduplicate/merge a TE consensus library with CD-HIT-EST at the 80-80-80 rule
// (params.cdhit_identity / params.cdhit_coverage).
// Stage 2 (LIBRARY_DEDUP): the de novo library plus an optional user library.
// Stage 3 (CURATION, as CD_HIT_EST_CURATED with ext.prefix '<id>.curated'): TEtrimmer's
// un-merged output, replacing TEtrimmer's own length-based merge (see DEDUP_REFINE).
// <prefix>.lineage.tsv records every sequence: id, origin (denovo | user), parent (-)
// and its ID in the input library. It is the only source of "which sequences are the
// user's" downstream (DEDUP_REFINE, TE_EVIDENCE).
process CD_HIT_EST {
    tag "$meta.id"
    label 'process_medium'

    input:
    // user_library: the user's records (Stage 2: --user_te_library; CURATION: those
    // records after TEtrimmer, kept as given), or []
    tuple val(meta), path(libraries, stageAs: 'lib_??.fa'), path(user_library, stageAs: 'user_library.fa')

    output:
    tuple val(meta), path("*.dedup_library.fa")       , emit: library
    // representatives + cluster file + all inputs + lineage, for DEDUP_REFINE
    tuple val(meta), path("*.dedup_library.fa"), path("*.dedup_library.fa.clstr"),
          path("*.combined_library.fa"), path("*.lineage.tsv"), emit: clustered
    tuple val(meta), path("*.lineage.tsv")            , emit: lineage
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
    def user = user_library ? 'user_library.fa' : ''
    """
    # De novo IDs (before '#') must already be unique: downstream steps key sequences,
    # evidence and provenance by them. A user ID that is already taken becomes <id>_user
    # (then <id>_user2, ...). Sources are told apart by file, so an empty de novo
    # library still leaves the user's sequences recognised as the user's.
    printf 'id\\torigin\\tparent\\tsource_id\\n' > ${prefix}.lineage.tsv
    awk -v user="${user}" -v lin="${prefix}.lineage.tsv" '
        /^>/ { id = substr(\$1, 2); b = id; lab = ""
               if ((i = index(id, "#")) > 0) { b = substr(id, 1, i - 1); lab = substr(id, i) }
               n = b
               if (FILENAME == user) {
                   k = 1
                   while (n in seen) { n = b "_user" (k > 1 ? k : ""); k++ }
                   if (n != b) printf "WARN: user library ID %s is already taken (de novo, or earlier in the user library); renamed %s\\n", b, n > "/dev/stderr"
                   origin = "user"
               } else {
                   if (n in seen) { printf "ERROR: ID %s occurs twice in the de novo library (%s); IDs before \\"#\\" must be unique\\n", b, FILENAME > "/dev/stderr"; exit 1 }
                   origin = "denovo"
               }
               seen[n] = 1
               print n "\\t" origin "\\t-\\t" b >> lin
               print ">" n lab substr(\$0, length(\$1) + 1); next }
        { print }' ${libraries} ${user} > ${prefix}.combined_library.fa

    if ! grep -q '^>' ${prefix}.combined_library.fa; then
        echo "WARN [${meta.id}]: no sequences to deduplicate; empty library" >&2
        : > ${prefix}.dedup_library.fa
        : > ${prefix}.dedup_library.fa.clstr
    else
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
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$(cd-hit-est -h 2>&1 | head -1 | sed -E 's/.*version ([0-9.]+).*/\\1/' || echo unknown)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def user = user_library ? 'user_library.fa' : ''
    """
    cat ${libraries} ${user} > ${prefix}.combined_library.fa
    cp ${prefix}.combined_library.fa ${prefix}.dedup_library.fa
    touch ${prefix}.dedup_library.fa.clstr
    printf 'id\\torigin\\tparent\\tsource_id\\n' > ${prefix}.lineage.tsv
    awk -v user="${user}" '/^>/ { b = substr(\$1, 2); sub(/#.*/, "", b)
                                  print b "\\t" (FILENAME == user ? "user" : "denovo") "\\t-\\t" b }' \\
        ${libraries} ${user} >> ${prefix}.lineage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: "stub"
    END_VERSIONS
    """
}
