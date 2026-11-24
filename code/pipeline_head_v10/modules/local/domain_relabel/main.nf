// Domain-label check for TEtrimmer-curated consensi (bin/te_evidence.py, no gating).
// A curated consensus is one rebuilt element, so when its TE protein domains (TEsorter,
// REXdb) point to another class than its label, the domain class is used: e.g. HELIANO
// "Helitrons" that carry an LTR/Gypsy or MuDR transposase domain. Unknown consensi with
// a domain take its class too. Uncurated consensi are left alone: there a chimera's
// domain describes only part of the sequence (run only in the TEtrimmer arm).
// Only real domain hits override a label, never TEsorter's homology-only second pass;
// the user's records (origin "user" in the lineage) keep their labels.
process DOMAIN_RELABEL {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(library), path(lineage)   // lineage: curated lineage, or []

    output:
    tuple val(meta), path("${meta.id}.domain_relabelled.fa"), emit: library
    tuple val(meta), path("${meta.id}.domain_relabel.tsv")  , emit: report
    path 'versions.yml'                                      , emit: versions

    script:
    def lin = lineage ? "--lineage ${lineage}" : ''
    """
    # A failed TEsorter run may leave a partial table; it is discarded, so a failure
    # changes no label (the same policy as TE_EVIDENCE), and the report says so.
    STATUS=""
    if ! TEsorter ${library} -db ${params.evidence_tesorter_db} -p ${task.cpus} \\
            -pre dr -tmp ./tesorter_tmp > tesorter.log 2>&1; then
        echo "WARN [${meta.id}]: TEsorter failed (see tesorter.log); domain labels not assessed, labels unchanged" >&2
        : > dr.cls.tsv
        STATUS="# TEsorter failed: domain labels were not assessed and every label was kept"
    fi
    touch dr.cls.tsv empty_hits.tsv

    te_evidence.py \\
        --library ${library} \\
        --raw ${library} ${lin} \\
        --tesorter-cls dr.cls.tsv \\
        --support-hits empty_hits.tsv \\
        --genome-hits empty_hits.tsv \\
        --mode off --relabel-unknown --relabel-conflicts \\
        --out-library ${meta.id}.domain_relabelled.fa \\
        --out-nonte /dev/null \\
        --out-removed /dev/null \\
        --evidence domain_table.tsv
    { [[ -n "\${STATUS}" ]] && echo "\${STATUS}"
      awk -F'\\t' 'NR == 1 || \$3 != \$4' domain_table.tsv | cut -f1,3,4,6; } > ${meta.id}.domain_relabel.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tesorter: \$(TEsorter -v 2>&1 | sed 's/TEsorter //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    cp ${library} ${meta.id}.domain_relabelled.fa
    printf 'name\\tlabel_in\\tlabel_out\\tdomain\\n' > ${meta.id}.domain_relabel.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tesorter: "stub"
    END_VERSIONS
    """
}
