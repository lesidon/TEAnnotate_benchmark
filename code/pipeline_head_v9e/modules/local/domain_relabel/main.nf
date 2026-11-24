// Domain-label check for TEtrimmer-curated consensi (bin/te_evidence.py, no gating).
// A curated consensus is one rebuilt element, so when its TE protein domains (TEsorter,
// REXdb) point to another class than its label, the domain class is used: e.g. HELIANO
// "Helitrons" that carry an LTR/Gypsy or MuDR transposase domain. Unknown consensi with
// a domain take its class too. Uncurated consensi are left alone: there a chimera's
// domain describes only part of the sequence (run only in the TEtrimmer arm).
process DOMAIN_RELABEL {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(library)

    output:
    tuple val(meta), path("${meta.id}.domain_relabelled.fa"), emit: library
    tuple val(meta), path("${meta.id}.domain_relabel.tsv")  , emit: report
    path 'versions.yml'                                      , emit: versions

    script:
    """
    TEsorter ${library} -db ${params.evidence_tesorter_db} -p ${task.cpus} \\
        -pre dr -tmp ./tesorter_tmp > tesorter.log 2>&1 \\
        || echo "WARN [${meta.id}]: TEsorter failed (see tesorter.log); labels unchanged" >&2
    touch dr.cls.tsv empty_hits.tsv

    te_evidence.py \\
        --library ${library} \\
        --raw ${library} \\
        --tesorter-cls dr.cls.tsv \\
        --support-hits empty_hits.tsv \\
        --genome-hits empty_hits.tsv \\
        --mode off --relabel-unknown --relabel-conflicts \\
        --out-library ${meta.id}.domain_relabelled.fa \\
        --out-nonte /dev/null \\
        --out-removed /dev/null \\
        --evidence domain_table.tsv
    awk -F'\\t' 'NR == 1 || \$3 != \$4' domain_table.tsv | cut -f1,3,4,6 > ${meta.id}.domain_relabel.tsv

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
