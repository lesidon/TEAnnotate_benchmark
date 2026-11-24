process REPEATMODELER2 {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/repeatmodeler:2.0.5--pl5321hdfd78af_0'
    conda     'bioconda::repeatmodeler=2.0.5'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.rm2_library.fa"), emit: fasta
    path "${meta.id}-families.stk"                    , emit: stockholm
    path "${meta.id}.repeatmodeler.log"               , emit: log
    path 'versions.yml'                               , emit: versions

    script:
    // -LTRStruct enables LTR retrotransposon detection; requires chromosome-scale assembly
    def ltr_struct = (meta.assembly_quality in ['chromosome', 'T2T']) ? '-LTRStruct' : ''
    """
    BuildDatabase \\
        -name "${meta.id}" \\
        -engine ncbi \\
        "${fasta}"

    RepeatModeler \\
        -database "${meta.id}" \\
        -threads ${task.cpus} \\
        ${ltr_struct} \\
        > "${meta.id}.repeatmodeler.log" 2>&1

    # RepeatModeler may find nothing on very small / low-complexity genomes
    if [[ -f "${meta.id}-families.fa" ]]; then
        cp "${meta.id}-families.fa" "${meta.id}.rm2_library.fa"
    else
        echo "WARN [${meta.id}]: RepeatModeler2 found no repeat families" >&2
        touch "${meta.id}.rm2_library.fa"
    fi

    [[ -f "${meta.id}-families.stk" ]] || touch "${meta.id}-families.stk"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version 2>&1 | sed 's/RepeatModeler version //')
    END_VERSIONS
    """

    stub:
    """
    cat > "${meta.id}.rm2_library.fa" <<-EOF
    >${meta.id}_stub_RM2#Unknown/Unknown
    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
    GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
    EOF

    touch "${meta.id}-families.stk"
    touch "${meta.id}.repeatmodeler.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: "stub"
    END_VERSIONS
    """
}
