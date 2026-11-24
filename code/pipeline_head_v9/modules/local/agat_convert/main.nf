// Normalize a gene annotation (GTF or GFF) to a standard GFF3 for Starfish.
// Starfish's gene finder parses `mRNA` features with GFF3-style ID/Parent
// attributes. BRAKER/AUGUSTUS/StringTie emit GTF, whose transcript-level feature
// is `transcript` with `gene_id`/`transcript_id` attributes — Starfish cannot
// parse it ("could not parse any mRNA features"). AGAT's gxf2gxf standardizes any
// GTF/GFF input into spec-compliant GFF3 (transcript -> mRNA, ID=/Parent=), which
// Starfish consumes. Running it on already-GFF3 input is harmless (it re-emits
// a cleaned GFF3), so the conversion is unconditional for annotated samples.
process GXF_TO_GFF3 {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(fasta), path(gff)

    output:
    tuple val(meta), path(fasta), path("${meta.id}.normalized.gff3"), emit: annotation
    path 'versions.yml'                                            , emit: versions

    script:
    """
    agat_convert_sp_gxf2gxf.pl \\
        --gff ${gff} \\
        --output ${meta.id}.agat.gff3

    # Starfish's gene finder lifts over `mRNA` features (default --targetFeat).
    # BRAKER/AUGUSTUS GTF uses `transcript` at the isoform level and AGAT keeps
    # that primary tag, so Starfish parses zero mRNA features. Rename the
    # third-column `transcript` -> `mRNA` (genes/CDS/exon untouched) so Starfish's
    # default matches. If AGAT already emitted `mRNA`, this is a no-op.
    awk 'BEGIN{FS=OFS="\\t"} /^#/{print; next} \$3=="transcript"{\$3="mRNA"} {print}' \\
        ${meta.id}.agat.gff3 > ${meta.id}.normalized.gff3

    AGAT_VER=\$(ls "\${CONDA_PREFIX:-}"/conda-meta/agat-*.json /usr/local/conda-meta/agat-*.json 2>/dev/null | sed -E 's@.*/agat-([0-9][^-]*)-.*@\\1@' | head -1 || true)
    [[ -z "\${AGAT_VER}" ]] && AGAT_VER=\$(agat_convert_sp_gxf2gxf.pl --help 2>&1 | grep -oiE 'v[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 | tr -d 'v' || true)
    [[ -z "\${AGAT_VER}" ]] && AGAT_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \${AGAT_VER}
    END_VERSIONS
    """

    stub:
    """
    cp ${gff} ${meta.id}.normalized.gff3
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: "stub"
    END_VERSIONS
    """
}
