include { CD_HIT_EST   } from '../../../modules/local/cd_hit_est/main'
include { DEDUP_REFINE } from '../../../modules/local/dedup_refine/main'

// Stage 2 — Library Deduplication & Merge.
// Combines the per-sample raw de novo library with an optional user-supplied TE
// library (--user_te_library) and collapses redundancy with CD-HIT-EST at the
// 80-80-80 rule, emitting one non-redundant library per sample. DEDUP_REFINE then
// repairs clusters that mix TE classes (--dedup_class_aware, on by default), so a
// family collapsed into a representative of another class is not lost, and, given
// the genome, makes the member with the most full-length genomic copies each
// cluster's representative instead of the longest (--dedup_rep_by_copies).
// The lineage (id, origin denovo | user, parent, source ID) records which sequences
// are the user's and where each chimera piece came from, for TE_EVIDENCE.
workflow LIBRARY_DEDUP {

    take:
    ch_raw_library  // channel: [ val(meta), path(raw_library.fa) ]
    ch_genome_by_id // channel: [ val(id), path(genome.fa) ]; may be empty

    main:
    ch_versions = Channel.empty()

    // Optional user library is a single shared file applied to every sample.
    def user_lib = params.user_te_library ? file(params.user_te_library, checkIfExists: true) : []

    CD_HIT_EST(ch_raw_library.map { meta, raw -> [ meta, [ raw ], user_lib ] })
    ch_versions = ch_versions.mix(CD_HIT_EST.out.versions.first())

    if (params.dedup_class_aware) {
        CD_HIT_EST.out.clustered
            .map { meta, dedup, clstr, combined, lineage -> [ meta.id, meta, dedup, clstr, combined, lineage ] }
            .join(ch_genome_by_id, by: 0, remainder: true)
            .filter { it[1] != null }
            .map { id, meta, dedup, clstr, combined, lineage, genome -> [ meta, dedup, clstr, combined, lineage, genome ?: [] ] }
            .set { ch_refine_input }
        DEDUP_REFINE(ch_refine_input)
        ch_versions = ch_versions.mix(DEDUP_REFINE.out.versions.first())
        ch_library  = DEDUP_REFINE.out.library
        ch_lineage  = DEDUP_REFINE.out.lineage
    } else {
        ch_library  = CD_HIT_EST.out.library
        ch_lineage  = CD_HIT_EST.out.lineage
    }

    emit:
    library  = ch_library  // channel: [ val(meta), path(*.{refined,dedup}_library.fa) ]
    lineage  = ch_lineage  // channel: [ val(meta), path(*lineage.tsv) ]
    versions = ch_versions
}
