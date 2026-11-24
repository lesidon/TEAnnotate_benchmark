// Stage 4 (visualization) — TE density karyoplot (chromosome painting).
//
// Bins the combined repeat BED into fixed windows (--karyoplot_window_kb) per scaffold
// and draws a Nature-styled density heatmap of the longest scaffolds. Scaffold lengths
// come straight from the genome FASTA (awk). bin/plot_karyotype.R is staged + published
// so the user can restyle (colours, scaffold count, window) and re-run.
process PLOT_KARYOTYPE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(genome), path(bed)
    path rscript

    output:
    tuple val(meta), path("${meta.id}.karyoplot.*"), emit: plot, optional: true
    path "${meta.id}.karyoplot_density.tsv"        , emit: density
    path "${meta.id}.karyoplot_lengths.tsv"        , emit: lengths
    path "${rscript}"                              , emit: rscript
    path 'versions.yml'                            , emit: versions

    script:
    """
    # Scaffold lengths from the FASTA (first header token = sequence name in the .out).
    awk '/^>/{ if(name!=""){print name"\\t"len} name=substr(\$1,2); len=0; next }
              { len+=length(\$0) }
         END{ if(name!=""){print name"\\t"len} }' ${genome} > lengths.body
    printf 'chrom\\tlength\\n'  > ${meta.id}.karyoplot_lengths.tsv
    cat lengths.body          >> ${meta.id}.karyoplot_lengths.tsv

    # Per-window masked fraction from the combined BED (capped at 1 for overlaps).
    WIN=\$(( ${params.karyoplot_window_kb} * 1000 ))
    awk -v w=\$WIN 'BEGIN{OFS="\\t"}
        { s=\$2; e=\$3; b=int(s/w); we=int((e-1)/w);
          for(i=b;i<=we;i++){ ws=i*w; wend=ws+w; os=(s>ws?s:ws); oe=(e<wend?e:wend);
                              if(oe>os){cov[\$1 SUBSEP i]+=oe-os} } }
        END{ for(k in cov){ split(k,a,SUBSEP); ws=a[2]*w; f=cov[k]/w; if(f>1)f=1;
                            print a[1], ws, ws+w, f } }' ${bed} > dens.body
    printf 'chrom\\tstart\\tend\\tfrac\\n' > ${meta.id}.karyoplot_density.tsv
    sort -k1,1 -k2,2n dens.body           >> ${meta.id}.karyoplot_density.tsv

    Rscript ${rscript} ${meta.id}.karyoplot_density.tsv ${meta.id}.karyoplot_lengths.tsv \\
        ${meta.id}.karyoplot "${params.plot_formats}" 30

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tlength\\n'            > ${meta.id}.karyoplot_lengths.tsv
    printf 'chrom\\tstart\\tend\\tfrac\\n' > ${meta.id}.karyoplot_density.tsv
    touch ${meta.id}.karyoplot.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: "stub"
    END_VERSIONS
    """
}
