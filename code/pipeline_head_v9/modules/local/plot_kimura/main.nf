// Stage 4 (visualization) — TE divergence (Kimura) landscape.
//
// Runs bin/plot_kimura.R on the coverage matrix from CALC_DIVERGENCE. The R script is
// staged (and published) alongside the figure so the user can restyle and re-run it.
// R container (rocker/tidyverse) preferred, conda-forge R fallback.
process PLOT_KIMURA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(matrix)
    path rscript

    output:
    tuple val(meta), path("${meta.id}.kimura_landscape.*"), emit: plot, optional: true
    path "${matrix}"                                       , emit: data
    path "${rscript}"                                      , emit: rscript
    path 'versions.yml'                                    , emit: versions

    script:
    def faceted = params.kimura_faceted ? '1' : '0'
    """
    Rscript ${rscript} ${matrix} ${meta.id}.kimura_landscape "${params.plot_formats}" ${faceted}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.kimura_landscape.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: "stub"
    END_VERSIONS
    """
}
