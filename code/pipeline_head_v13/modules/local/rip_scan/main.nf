// Stage 3 (masking, fungal genomes only) — RIP (repeat-induced point mutation) scan.
//
// Runs bin/rip_scan.py on the genome and the RepeatMasker .align of the final masking
// run: RIP-type C->T mutations per TE copy against its consensus, RIP indices per copy
// and per family, and large RIP-affected genomic regions (LRARs) in sliding windows.
// A report only; nothing downstream consumes it. Thresholds can be changed through
// ext.args (see the script's --help).
process RIP_SCAN {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(genome), path(align)

    output:
    tuple val(meta), path("${meta.id}.rip_*"), emit: report
    path 'versions.yml'                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    rip_scan.py --genome ${genome} --align ${align} --prefix ${meta.id} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    for f in rip_copies.tsv rip_families.tsv rip_summary.tsv rip_lrar.bed rip_cri.bedgraph; do
        touch ${meta.id}.\$f
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
