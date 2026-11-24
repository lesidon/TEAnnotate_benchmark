// Class-aware refinement of a CD-HIT-EST library (bin/refine_dedup.py).
// Cluster step: relabels Unknown representatives from their classified members and
// rescues members of another class as their own families (CD-HIT keeps each cluster's
// longest sequence with its own label, so they would otherwise be lost).
// Alignment step (--dedup_split_chimeras): an all-against-all blastn finds consensi
// that carry a segment of another TE class (e.g. an LTR-int spanning a LINE inserted
// next to it) and cuts that segment out. Needs blastn, hence the TEtrimmer image.
// Representative step (with the genome, --dedup_rep_by_copies): CD-HIT keeps each
// cluster's longest member, which may be an over-extended or chimeric consensus that
// absorbed an exact one. The member of the representative's class with the most
// full-length genomic copies (>= 75% identity over >= 90% of it) becomes the
// representative instead; ties keep CD-HIT's choice.
// User step (with a user library): a cluster holding one of the user's sequences keeps
// that sequence (the longest, if several) as its representative, so the user's curated
// record survives deduplication and stays trusted by the evidence gate.
// Lineage: CD_HIT_EST's record of every sequence (origin denovo | user) in, the same
// plus one row per chimera piece (parent = the host) out, for TE_EVIDENCE.
// Used after Stage 2 CD_HIT_EST and, as DEDUP_REFINE_CURATED, after TEtrimmer.
process DEDUP_REFINE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(dedup), path(clstr), path(combined), path(lineage), path(genome)

    output:
    tuple val(meta), path("*.refined_library.fa")  , emit: library
    tuple val(meta), path("*.dedup_refine.tsv")    , emit: report
    tuple val(meta), path("*.refined_lineage.tsv") , emit: lineage
    path 'versions.yml'                            , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def split  = params.dedup_split_chimeras ? '' : '--no-split'
    def copies = (genome && params.dedup_rep_by_copies) ? "--genome ${genome}" : ''
    """
    refine_dedup.py \\
        --combined ${combined} \\
        --dedup ${dedup} \\
        --clstr ${clstr} \\
        --lineage ${lineage} \\
        --out ${prefix}.refined_library.fa \\
        --out-lineage ${prefix}.refined_lineage.tsv \\
        --report ${prefix}.dedup_refine.tsv \\
        --threads ${task.cpus} \\
        ${split} ${copies}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        blast: \$(blastn -version 2>/dev/null | head -1 | sed 's/blastn: //' || echo none)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${dedup} ${prefix}.refined_library.fa
    cp ${lineage} ${prefix}.refined_lineage.tsv
    printf 'step\\taction\\tsequence\\tresult\\tdetail\\n' > ${prefix}.dedup_refine.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
