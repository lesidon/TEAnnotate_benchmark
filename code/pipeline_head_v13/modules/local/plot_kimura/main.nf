// Stage 4 (visualization) — TE divergence (Kimura) landscape.
//
// Runs bin/plot_kimura.R on the coverage matrix from CALC_DIVERGENCE. The R script is
// staged (and published) alongside the figure so the user can restyle and re-run it.
// The genome length for the % of genome y axis is the 'total length' row of the
// RepeatMasker .tbl (as in REPEAT_QC); without it the axis falls back to Mb.
// R container (rocker/tidyverse) preferred, conda-forge R fallback.
process PLOT_KIMURA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(matrix), path(tbl)
    path rscript

    output:
    tuple val(meta), path("${meta.id}.kimura_landscape*"), emit: plot, optional: true
    path "${matrix}"                                       , emit: data
    path "${rscript}"                                      , emit: rscript
    path 'versions.yml'                                    , emit: versions

    script:
    """
    GENOME_BP=\$(awk '/total length:/ { for (i = 1; i <= NF; i++) if (\$i ~ /^[0-9]+\$/) { print \$i; exit } }' ${tbl} 2>/dev/null || true)
    Rscript ${rscript} ${matrix} ${meta.id}.kimura_landscape "${params.plot_formats}" ${params.kimura_layout} "\${GENOME_BP:-0}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.kimura_landscape.png ${meta.id}.kimura_landscape_faceted.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: "stub"
    END_VERSIONS
    """
}
