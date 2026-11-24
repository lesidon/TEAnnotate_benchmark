process ASSESS_ASSEMBLY_QUALITY {
    tag "$meta.id"
    label 'process_single'

    // Only Python stdlib is needed (csv, json, argparse) — any slim Python 3 image works.
    container 'docker.io/python:3.11'
    conda     'conda-forge::python=3.11'

    input:
    tuple val(meta), path(fasta), path(stats_tsv)

    output:
    tuple val(meta), path(fasta), path("${meta.id}.stats.json"), emit: genome_with_stats
    path 'versions.yml'                                        , emit: versions

    script:
    def quality_override = (params.assembly_quality != 'auto') ? params.assembly_quality : ''
    """
    assess_assembly_quality.py \\
        --stats        "${stats_tsv}" \\
        --override     "${quality_override}" \\
        --n50-floor-mb ${params.chromosome_n50_mb} \\
        --n50-fraction ${params.chromosome_n50_fraction} \\
        --out          "${meta.id}.stats.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat > "${meta.id}.stats.json" <<-EOF
    {
      "assembly_quality": "draft",
      "quality_rule": "below_thresholds",
      "genome_size_mb": 10.0,
      "num_seqs": 1,
      "n50_mb": 10.0,
      "n50_to_genome_ratio": 0.05,
      "max_seq_mb": 10.0,
      "gc_percent": 48.5
    }
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
