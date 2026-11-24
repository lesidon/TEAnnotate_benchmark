// Stage 3 — post-process TEtrimmer's un-merged output (bin/tetrimmer_postprocess.py)
// before the class-aware merge (CD_HIT_EST_CURATED + DEDUP_REFINE_CURATED), which
// replaces TEtrimmer's own length-based merge. Drops Need_check consensi and re-checks
// labels: TEtrimmer keeps an input's label on every consensus it derives from it, even
// when the consensus was rebuilt from another element (e.g. a chimeric LTR-int comes
// back as the LINE inside it), so such a consensus takes the label of the input it
// actually matches. With TE_EVIDENCE's table, inputs TEtrimmer dropped return
// unchanged when they carry a TE domain or another tool's support: TEtrimmer keeps a
// skipped input only with >=2 full-length copies, which discards real low-copy and
// degraded families; inputs the gate kept unassessed (failed domain search) return too.
// With the Stage 2 lineage, the user's records leave as given (<id>.user_records.fa),
// consensi TEtrimmer rebuilt from them are dropped, and the caller merges them back as
// user records. <id>.tetrimmer_lineage.tsv gives each consensus its input as parent.
// Needs blastn, hence the TEtrimmer image.
process TETRIMMER_POSTPROCESS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(tetrimmer_dir), path(input_library), path(evidence), path(lineage)

    output:
    tuple val(meta), path("${meta.id}.tetrimmer_consensus.fa")    , emit: library
    tuple val(meta), path("${meta.id}.user_records.fa")           , emit: user
    tuple val(meta), path("${meta.id}.tetrimmer_lineage.tsv")     , emit: lineage
    tuple val(meta), path("${meta.id}.tetrimmer_postprocess.tsv") , emit: report
    path 'versions.yml'                                           , emit: versions

    script:
    def reinject = (evidence && params.tetrimmer_reinject) ? "--evidence ${evidence} --reinject-tiers 1,2,T" : ''
    def lin      = lineage ? "--lineage ${lineage}" : ''
    """
    tetrimmer_postprocess.py \\
        --tetrimmer-dir ${tetrimmer_dir} \\
        --input ${input_library} \\
        --out ${meta.id}.tetrimmer_consensus.fa \\
        --out-user ${meta.id}.user_records.fa \\
        --out-lineage ${meta.id}.tetrimmer_lineage.tsv \\
        --report ${meta.id}.tetrimmer_postprocess.tsv \\
        --threads ${task.cpus} ${reinject} ${lin}

    [[ -s ${meta.id}.tetrimmer_consensus.fa ]] || \\
        echo "WARN [${meta.id}]: no TEtrimmer consensus left after Need_check filtering" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        blast: \$(blastn -version 2>/dev/null | head -1 | sed 's/blastn: //' || echo none)
    END_VERSIONS
    """

    stub:
    // Sequences flow through (user records by the lineage), so a stub run can check
    // that families, and the user's records in particular, reach the final library.
    def lin = lineage ? "${lineage}" : ''
    """
    awk -v lin="${lin}" '
        FILENAME == lin { if (FNR > 1 && \$2 == "user") u[\$1] = 1; next }
        /^>/ { b = substr(\$1, 2); sub(/#.*/, "", b); out = (b in u) ? "${meta.id}.user_records.fa" : "${meta.id}.tetrimmer_consensus.fa" }
        { print > out }' ${lin} ${input_library}
    touch ${meta.id}.tetrimmer_consensus.fa ${meta.id}.user_records.fa
    printf 'id\\torigin\\tparent\\tsource_id\\n' > ${meta.id}.tetrimmer_lineage.tsv
    printf 'consensus\\tinput\\ttier\\taction\\tlabel_in\\tlabel_out\\tdetail\\n' > ${meta.id}.tetrimmer_postprocess.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
