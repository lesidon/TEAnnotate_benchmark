include { INPUT_CHECK              } from '../../../modules/local/input_check/main'
include { GENOME_STATS             } from '../../../modules/local/genome_stats/main'
include { ASSESS_ASSEMBLY_QUALITY  } from '../../../modules/local/assess_assembly_quality/main'

workflow INPUT_VALIDATION {

    take:
    ch_input  // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    // ── Step 1: FASTA format sanity check ─────────────────────────────────────
    INPUT_CHECK(ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions.first())

    // ── Step 2: Compute genome statistics with seqkit ─────────────────────────
    GENOME_STATS(INPUT_CHECK.out.genome)
    ch_versions = ch_versions.mix(GENOME_STATS.out.versions.first())

    // ── Step 3: Derive assembly quality tier from seqkit TSV (Python stdlib) ──
    ASSESS_ASSEMBLY_QUALITY(GENOME_STATS.out.genome_with_tsv)
    ch_versions = ch_versions.mix(ASSESS_ASSEMBLY_QUALITY.out.versions.first())

    // ── Step 4: Parse JSON stats → update meta map ────────────────────────────
    ch_validated = ASSESS_ASSEMBLY_QUALITY.out.genome_with_stats
        .map { meta, fasta, stats_json ->
            def stats   = new groovy.json.JsonSlurper().parseText(stats_json.text)
            def quality = stats.assembly_quality as String

            def updated_meta = meta + [
                assembly_quality : quality,
                genome_size_mb   : stats.genome_size_mb   as Float,
                num_seqs         : stats.num_seqs         as Long,
                n50_mb           : stats.n50_mb           as Float,
                gc_percent       : stats.gc_percent       as Float,
            ]

            log.info(
                "[${meta.id}] assembly_quality=${quality} " +
                "| size=${stats.genome_size_mb} Mb " +
                "| N50=${stats.n50_mb} Mb " +
                "| seqs=${stats.num_seqs} " +
                "| GC=${stats.gc_percent}%"
            )

            [ updated_meta, fasta ]
        }

    emit:
    genome   = ch_validated                  // channel: [ val(meta+stats), path(fasta) ]
    stats    = GENOME_STATS.out.stats_tsv    // channel: path(*.seqkit_stats.tsv)
    versions = ch_versions                   // channel: path(versions.yml)
}
