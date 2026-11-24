process GENOME_STATS {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/seqkit:2.8.0--h9ee0642_0'
    conda     'bioconda::seqkit=2.8.0'

    input:
    tuple val(meta), path(fasta)

    output:
    path "${meta.id}.seqkit_stats.tsv", emit: stats_tsv
    path 'versions.yml'               , emit: versions

    script:
    """
    seqkit stats \\
        --all \\
        --tabular \\
        --threads ${task.cpus} \\
        "${fasta}" > "${meta.id}.seqkit_stats.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\tQ1\\tQ2\\tQ3\\tsum_gap\\tN50\\tQ20(%%)\\tQ30(%%)\\tGC(%%)\\n' \\
        > "${meta.id}.seqkit_stats.tsv"
    printf '%s\\tFASTA\\tDNA\\t1\\t10000000\\t10000000\\t10000000.0\\t10000000\\t10000000\\t10000000\\t10000000\\t0\\t10000000\\t0.00\\t0.00\\t48.50\\n' \\
        "${fasta}" >> "${meta.id}.seqkit_stats.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: "stub"
    END_VERSIONS
    """
}
