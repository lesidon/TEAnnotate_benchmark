// Evidence gate for the deduplicated library (bin/te_evidence.py).
// A union of discovery tools is sensitive but carries every tool's false positives;
// in the benchmark these sit in a minority of consensi that mask almost no reference
// TE. Each consensus is scored on independent evidence:
//   * TE protein domains: TEsorter (REXdb clade profiles) and, with Pfam-A, the
//     TE-specific Pfam families on six-frame ORF fragments (bin/te_pfam.py); Pfam's
//     broad family models also catch divergent TEs of unexplored lineages,
//   * the same element in another tool's library (blastn against the raw union),
//   * genomic copies (blastn against the genome): >=2 full-length (>=75% identity over
//     >=90% of the consensus, TEtrimmer's definition) or >=5 covering half of it. Old
//     families decay by mutation and deletion and keep partial copies only; RIP makes
//     this extreme in fungi, and in unexplored genomes (e.g. sponges) such families
//     have neither a known domain nor a second tool to vouch for them,
// and the gate (--evidence_gate lite|strict) drops non-TE consensi (Satellite, rRNA,
// ...; written to *.non_te_repeats.fa) and consensi without evidence. Unknown
// consensi with a TE domain take its class (--evidence_relabel_unknown). Sequences of
// the user's curated library (--user_te_library) are never gated: the lineage from
// LIBRARY_DEDUP says which they are (origin "user", including pieces cut from them).
process TE_EVIDENCE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(library), path(raw), path(genome), path(provenance), path(lineage)
    path(pfam_db)      // dir with Pfam-A.hmm, or [] to skip the Pfam domain search

    output:
    tuple val(meta), path("${meta.id}.gated_library.fa")  , emit: library
    tuple val(meta), path("${meta.id}.evidence.tsv")      , emit: evidence
    tuple val(meta), path("${meta.id}.non_te_repeats.fa") , emit: non_te
    tuple val(meta), path("${meta.id}.gate_removed.fa")   , emit: removed
    path 'versions.yml'                                   , emit: versions

    script:
    def prov    = provenance ? "--provenance ${provenance}" : ''
    def lin     = lineage ? "--lineage ${lineage}" : ''
    def relabel = params.evidence_relabel_unknown ? '--relabel-unknown' : ''
    def pfam_hmm = (pfam_db && params.evidence_pfam) ? "${pfam_db}/Pfam-A.hmm" : ''
    """
    # An empty library (no de novo family and no user library) has nothing to score.
    if ! grep -q '^>' ${library}; then
        echo "WARN [${meta.id}]: empty library; nothing to score" >&2
        : > ${meta.id}.gated_library.fa; : > ${meta.id}.non_te_repeats.fa; : > ${meta.id}.gate_removed.fa
        printf 'name\\ttool\\tlabel_in\\tlabel_out\\tlength\\tdomain\\tsupport_tools\\tfull_copies\\ttier\\taction\\thalf_copies\\tpfam\\thomology\\n' > ${meta.id}.evidence.tsv
        printf '"%s":\\n    te_evidence: "empty library"\\n' "${task.process}" > versions.yml
        exit 0
    fi

    # TE protein domains. TEsorter writes nothing for a library without domain hits;
    # an empty table then simply gives no tier-1 evidence. A FAILED run is different:
    # domain evidence is then missing, not negative, so its (possibly partial) table is
    # discarded and the gate keeps consensi it would otherwise remove for lack of evidence.
    UNASSESSED=""
    if ! TEsorter ${library} -db ${params.evidence_tesorter_db} -p ${task.cpus} \\
            -pre ev -tmp ./tesorter_tmp > tesorter.log 2>&1; then
        echo "WARN [${meta.id}]: TEsorter failed (see tesorter.log); domain evidence unavailable," >&2
        echo "  so no consensus is removed for lack of evidence (action keep_unassessed)" >&2
        : > ev.cls.tsv
        UNASSESSED="--domains-unassessed"
    fi
    touch ev.cls.tsv

    # TE-specific Pfam domains on six-frame ORF fragments (>= 80 aa, stop to stop).
    PFAM_ARG=""
    if [[ -n "${pfam_hmm}" && -s "${pfam_hmm}" ]]; then
        te_pfam.py subset ${pfam_hmm} te_pfam.hmm
        te_pfam.py orfs ${library} orfs.faa
        if [[ -s orfs.faa ]]; then
            hmmsearch --cpu ${task.cpus} --noali -E 1e-3 --domtblout pfam.domtbl te_pfam.hmm orfs.faa > /dev/null
            PFAM_ARG="--pfam-hits pfam.domtbl"
        else
            echo "INFO [${meta.id}]: no ORF fragment >= 80 aa; no Pfam domain search" >&2
        fi
    else
        echo "WARN [${meta.id}]: no Pfam-A.hmm; domain evidence from TEsorter only" >&2
    fi

    # The same element in another tool's library (support is counted per tool). With no
    # de novo sequence (only the user's library) there is nothing to search against.
    if grep -q '^>' ${raw}; then
        makeblastdb -in ${raw} -dbtype nucl -out rawdb > /dev/null
        blastn -task dc-megablast -query ${library} -db rawdb \\
            -outfmt '6 qseqid sseqid pident length qstart qend qlen' \\
            -evalue 1e-10 -perc_identity ${params.evidence_support_identity} -dust no \\
            -max_target_seqs 1000 -num_threads ${task.cpus} > support_hits.tsv
    else
        echo "INFO [${meta.id}]: empty de novo library; no cross-tool support search" >&2
        : > support_hits.tsv
    fi

    # Genomic copies covering >=50% of the consensus at >=70% identity; te_evidence.py
    # counts full-length ones (>=75% over >=90%) and half-length ones from the same hits,
    # one per genomic locus (overlapping HSPs on one subject sequence count once).
    makeblastdb -in ${genome} -dbtype nucl -out genomedb > /dev/null
    blastn -task dc-megablast -query ${library} -db genomedb \\
        -outfmt '6 qseqid sseqid pident length qstart qend qlen sstart send' \\
        -evalue 1e-20 -perc_identity 70 -qcov_hsp_perc 50 \\
        -max_target_seqs 100 -max_hsps 10 -num_threads ${task.cpus} > genome_hits.tsv

    te_evidence.py \\
        --library ${library} \\
        --raw ${raw} ${prov} ${lin} \\
        --tesorter-cls ev.cls.tsv \${PFAM_ARG} \${UNASSESSED} \\
        --support-hits support_hits.tsv \\
        --genome-hits genome_hits.tsv \\
        --mode ${params.evidence_gate} ${relabel} \\
        --support-id ${params.evidence_support_identity} \\
        --support-cov ${params.evidence_support_coverage} \\
        --min-full-copies ${params.evidence_min_full_copies} \\
        --min-half-copies ${params.evidence_min_half_copies} \\
        --out-library ${meta.id}.gated_library.fa \\
        --out-nonte ${meta.id}.non_te_repeats.fa \\
        --out-removed ${meta.id}.gate_removed.fa \\
        --evidence ${meta.id}.evidence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tesorter: \$(TEsorter -v 2>&1 | sed 's/TEsorter //' | head -1)
        hmmer: \$(hmmsearch -h | grep -oE 'HMMER [0-9.]+' | head -1 | sed 's/HMMER //')
        blast: \$(blastn -version 2>/dev/null | head -1 | sed 's/blastn: //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${library} ${meta.id}.gated_library.fa
    printf 'name\\ttool\\tlabel_in\\tlabel_out\\tlength\\tdomain\\tsupport_tools\\tfull_copies\\ttier\\taction\\thalf_copies\\tpfam\\thomology\\n' > ${meta.id}.evidence.tsv
    touch ${meta.id}.non_te_repeats.fa ${meta.id}.gate_removed.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tesorter: "stub"
    END_VERSIONS
    """
}
