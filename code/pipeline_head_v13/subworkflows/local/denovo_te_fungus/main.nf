include { EARLGREY          } from '../../../modules/local/earlgrey/main'
include { HITE              } from '../../../modules/local/hite/main'
include { STARFISH_ANNOTATE } from '../../../modules/local/starfish/main'
include { CAT_TE_LIBRARIES  } from '../../../modules/local/cat_te_libraries/main'

// De novo TE discovery for fungal genomes. Per-genome consensus libraries from
// these engines are merged before Stage 2 dedup:
//   * Earl Grey (library-only) — RepeatModeler2 + BEAT curation. Always runs.
//   * HiTE (library-only, --annotate 0) — full-length LTR/TIR/Helitron/non-LTR.
//     Runs by default; disable with --skip_hite.
// Starfish (Starship captains; off with --skip_starfish) is a report only: its
// captain table and raw outputs are published to 02_denovo_te/starfish/ and never
// enter the library, so they are not curated, classified or used for masking.
workflow DENOVO_TE_FUNGUS {

    take:
    ch_genome       // channel: [ val(meta), path(fasta) ]
    ch_ltrdetector  // value channel: path(LtrDetector dir) for HiTE; [] with --skip_hite

    main:
    ch_versions = Channel.empty()

    // Earl Grey runs on the genome alone, for every fungal sample.
    EARLGREY(ch_genome)
    ch_versions = ch_versions.mix(EARLGREY.out.versions.first())

    ch_libraries = EARLGREY.out.library

    // HiTE runs by default on the genome alone (organism-agnostic).
    if (!params.skip_hite) {
        HITE(ch_genome, ch_ltrdetector)
        ch_versions  = ch_versions.mix(HITE.out.versions.first())
        ch_libraries = ch_libraries.mix(HITE.out.library)
    }

    // Starfish runs by default; with --skip_starfish the process is never
    // instantiated and its Conda env is never built. Its outputs are only published
    // (report), never mixed into ch_libraries.
    if (!params.skip_starfish) {
        STARFISH_ANNOTATE(ch_genome)
        ch_versions = ch_versions.mix(STARFISH_ANNOTATE.out.versions.first())
    }

    // Collect all per-sample libraries [meta, lib] into [meta, [libs]] for merging.
    // Earl Grey is always present; HiTE contributes when active.
    ch_libraries
        .groupTuple()
        .set { ch_to_merge }

    CAT_TE_LIBRARIES(ch_to_merge)
    ch_versions = ch_versions.mix(CAT_TE_LIBRARIES.out.versions.first())

    emit:
    raw_library = CAT_TE_LIBRARIES.out.library     // channel: [ val(meta), path(*.raw_library.fa) ]
    provenance  = CAT_TE_LIBRARIES.out.provenance  // channel: [ val(meta), path(*.provenance.tsv) ]
    versions    = ch_versions
}
