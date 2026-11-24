include { INPUT_CHECK              } from '../../../modules/local/input_check/main'
include { GENOME_STATS             } from '../../../modules/local/genome_stats/main'

workflow INPUT_VALIDATION {

    take:
    ch_input  // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    // ── Step 1: FASTA format sanity check ─────────────────────────────────────
    INPUT_CHECK(ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions.first())

    // ── Step 2: Genome statistics with seqkit (a QC report; no step depends on it) ─
    GENOME_STATS(INPUT_CHECK.out.genome)
    ch_versions = ch_versions.mix(GENOME_STATS.out.versions.first())

    emit:
    genome   = INPUT_CHECK.out.genome        // channel: [ val(meta), path(fasta) ]
    stats    = GENOME_STATS.out.stats_tsv    // channel: path(*.seqkit_stats.tsv)
    versions = ch_versions                   // channel: path(versions.yml)
}
