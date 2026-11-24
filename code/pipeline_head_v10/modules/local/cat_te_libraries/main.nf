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
    # IDs must be unique across the tools' libraries: downstream steps key sequences,
    # provenance and evidence by the ID before '#'. A repeated ID becomes <id>_<tool>
    # (then <id>_<tool>2, ...), in the library and the provenance table alike. The tools'
    # libraries are taken in name order, not in the order the tools finished, so the
    # same inputs always give the same IDs.
    : > seen_ids.txt
    for f in \$(printf '%s\\n' ${libraries} | sort); do
        tool=\$(basename "\$f" .fa)
        tool=\${tool#${meta.id}.}
        tool=\${tool%_library}
        awk -v t="\$tool" -v lib="${meta.id}.raw_library.fa" -v prov="${meta.id}.provenance.tsv" '
            FILENAME == "seen_ids.txt" { seen[\$1] = 1; next }
            /^>/ { id = substr(\$1, 2); b = id; lab = ""
                   if ((i = index(id, "#")) > 0) { b = substr(id, 1, i - 1); lab = substr(id, i) }
                   n = b; k = 1
                   while (n in seen) { n = b "_" t (k > 1 ? k : ""); k++ }
                   if (n != b) { printf "WARN: duplicate ID %s in the %s library renamed %s\\n", b, t, n > "/dev/stderr" }
                   seen[n] = 1; print n >> "seen_ids.txt"
                   print n "\\t" t >> prov
                   print ">" n lab substr(\$0, length(\$1) + 1) >> lib
                   next }
            { print >> lib }' seen_ids.txt "\$f"
    done
    rm -f seen_ids.txt

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
    cat \$(printf '%s\\n' ${libraries} | sort) > "${meta.id}.raw_library.fa"
    touch "${meta.id}.provenance.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: "stub"
    END_VERSIONS
    """
}
