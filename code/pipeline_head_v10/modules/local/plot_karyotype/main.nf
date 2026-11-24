// Stage 4 (visualization) — TE density karyoplot (chromosome painting).
//
// Bins the interspersed repeats of the combined BED (simple / low-complexity / satellite
// repeats left out, as in te.gff3) into fixed windows (--karyoplot_window_kb) per
// scaffold and draws a Nature-styled density heatmap of the longest scaffolds. Each
// window's value is the union of its masked bases (overlapping annotations count once)
// over its own length (the last window of a scaffold is shorter); windows without
// repeats are 0. Scaffold lengths come straight from the genome FASTA (awk).
// bin/plot_karyotype.R is staged + published so the user can restyle (colours,
// scaffold count, window) and re-run.
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

    # Interspersed repeats only (BED col 4 = repeat#class), merged into non-overlapping
    # intervals per scaffold, so a base masked twice counts once.
    awk -F'\\t' 'BEGIN { OFS = "\\t" }
        { split(\$4, a, "#"); if (a[2] !~ /^(Simple_repeat|Low_complexity|Satellite)/) print \$1, \$2, \$3 }' \\
        ${bed} | sort -k1,1 -k2,2n -k3,3n \\
        | awk 'BEGIN { OFS = "\\t" }
               \$1 != c || \$2 > e { if (c != "") print c, s, e; c = \$1; s = \$2; e = \$3; next }
               \$3 > e { e = \$3 }
               END { if (c != "") print c, s, e }' > merged.bed

    # Masked fraction per window: covered bases over the window's length, the last
    # window clipped to the scaffold end; every window of every scaffold is written.
    WIN=\$(( ${params.karyoplot_window_kb} * 1000 ))
    awk -v w=\$WIN 'BEGIN { OFS = "\\t" }
        FNR == NR { len[\$1] = \$2; order[++n] = \$1; next }
        { s = \$2; e = \$3; if (\$1 in len && e > len[\$1]) e = len[\$1]
          for (i = int(s / w); i * w < e; i++) {
              ws = i * w; we = ws + w; os = (s > ws ? s : ws); oe = (e < we ? e : we)
              if (oe > os) cov[\$1 SUBSEP i] += oe - os } }
        END { for (j = 1; j <= n; j++) { c = order[j]
                  for (i = 0; i * w < len[c]; i++) {
                      ws = i * w; we = (ws + w < len[c] ? ws + w : len[c]); k = c SUBSEP i
                      print c, ws, we, (k in cov ? cov[k] : 0) / (we - ws) } } }' \\
        lengths.body merged.bed > dens.body
    printf 'chrom\\tstart\\tend\\tfrac\\n' > ${meta.id}.karyoplot_density.tsv
    cat dens.body                         >> ${meta.id}.karyoplot_density.tsv

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
