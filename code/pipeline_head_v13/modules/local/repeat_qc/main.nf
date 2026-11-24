// Stage 4 (visualization) — repeat detection & masking quality control.
//
// Derives QC metrics from the RepeatMasker .tbl (genome length, % masked), the combined
// repeat BED (per-class bp -> % of genome), and the classified library (family count,
// unclassified fraction). A high unclassified fraction is flagged WARN (weak
// classification signal). bin/repeat_qc.R draws the Nature-styled composition figure and
// echoes the QC table; the script + tables are published so the user can extend them.
process REPEAT_QC {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(tbl), path(bed), path(library)
    path rscript

    output:
    tuple val(meta), path("${meta.id}.repeat_qc_composition.*"), emit: plot, optional: true
    path "${meta.id}.repeat_qc_summary.tsv"                    , emit: summary
    path "${meta.id}.repeat_qc_composition.tsv"               , emit: composition
    path "${rscript}"                                         , emit: rscript
    path 'versions.yml'                                       , emit: versions

    script:
    """
    # ── Genome length + % masked from the RepeatMasker .tbl ──────────────────────
    # awk, not grep: a .tbl without these rows (or none at all) must give "unknown",
    # not abort the task under pipefail.
    TOTAL=\$(awk '/total length:/ { for (i = 1; i <= NF; i++) if (\$i ~ /^[0-9]+\$/) { print \$i; exit } }' ${tbl} 2>/dev/null || true)
    [[ -z "\${TOTAL}" ]] && TOTAL=NA
    MASKEDPCT=\$(awk -F'[()%]' '/bases masked:/{gsub(/ /,"",\$2); print \$2; exit}' ${tbl} 2>/dev/null || true)
    [[ -z "\${MASKEDPCT}" ]] && MASKEDPCT=NA

    # ── Per top-level class bp from the combined BED (col4 = repeat#class) ────────
    # The union of each class's intervals: a base covered by two hits of one class
    # counts once (a base covered by two classes counts in both).
    awk -F'\\t' 'BEGIN { OFS = "\\t" }
        { n = split(\$4, p, "#"); cls = (n >= 2) ? p[2] : p[1]; split(cls, c, "/")
          top = c[1]; if (top == "") top = "Unknown"; print top, \$1, \$2, \$3 }' ${bed} \\
        | sort -t \$'\\t' -k1,1 -k2,2 -k3,3n \\
        | awk -F'\\t' '\$1 != t || \$2 != c || \$3 > e { if (t != "") bp[t] += e - s; t = \$1; c = \$2; s = \$3; e = \$4; next }
                      \$4 > e { e = \$4 }
                      END { if (t != "") bp[t] += e - s; for (k in bp) printf "%s\\t%d\\n", k, bp[k] }' > class_bp.tsv

    printf 'class\\tbp\\tpct\\n' > ${meta.id}.repeat_qc_composition.tsv
    awk -v tot="\${TOTAL}" 'BEGIN{FS=OFS="\\t"} { p = (tot != "NA" && tot > 0) ? sprintf("%.4f", \$2 / tot * 100) : "NA";
        printf "%s\\t%d\\t%s\\n", \$1, \$2, p }' class_bp.tsv \\
        | sort -t \$'\\t' -k2,2nr >> ${meta.id}.repeat_qc_composition.tsv

    # ── Library family / unclassified counts ─────────────────────────────────────
    # grep -c already prints 0 (and exits 1) when nothing matches, so DON'T add
    # '|| echo 0' — that would append a second 0 and break the unquoted awk -v below.
    FAM=\$(grep -c '^>' ${library} 2>/dev/null || true); FAM=\${FAM:-0}
    UNK=\$(grep '^>' ${library} 2>/dev/null | grep -c '#Unknown' || true); UNK=\${UNK:-0}
    NCLS=\$(( \$(wc -l < ${meta.id}.repeat_qc_composition.tsv) - 1 ))
    if [[ "\${FAM}" -gt 0 ]]; then
        UNKPCT=\$(awk -v u="\${UNK}" -v f="\${FAM}" 'BEGIN{printf "%.2f", u/f*100}')
    else
        UNKPCT=NA
    fi
    UNKSTATUS=PASS
    awk -v u="\${UNKPCT}" 'BEGIN{ exit (u!="NA" && u+0>25) ? 0 : 1 }' && UNKSTATUS=WARN

    {
        printf 'metric\\tvalue\\tstatus\\n'
        printf 'genome_length_bp\\t%s\\tNA\\n'      "\${TOTAL}"
        printf 'bases_masked_pct\\t%s\\tNA\\n'      "\${MASKEDPCT}"
        printf 'repeat_classes\\t%s\\tNA\\n'        "\${NCLS}"
        printf 'repeat_families\\t%s\\tNA\\n'       "\${FAM}"
        printf 'unclassified_families\\t%s\\tNA\\n' "\${UNK}"
        printf 'unclassified_pct\\t%s\\t%s\\n'      "\${UNKPCT}" "\${UNKSTATUS}"
    } > ${meta.id}.repeat_qc_summary.tsv

    Rscript ${rscript} ${meta.id}.repeat_qc_composition.tsv ${meta.id}.repeat_qc_summary.tsv \\
        ${meta.id}.repeat_qc "${params.plot_formats}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'metric\\tvalue\\tstatus\\n'  > ${meta.id}.repeat_qc_summary.tsv
    printf 'class\\tbp\\tpct\\n'         > ${meta.id}.repeat_qc_composition.tsv
    touch ${meta.id}.repeat_qc_composition.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-ggplot2: "stub"
    END_VERSIONS
    """
}
