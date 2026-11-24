include { EARLGREY         } from '../../../modules/local/earlgrey/main'
include { HITE             } from '../../../modules/local/hite/main'
include { CAT_TE_LIBRARIES } from '../../../modules/local/cat_te_libraries/main'

// De novo TE discovery for the animal route (also the default when organism_type
// is unspecified).
// Two de novo engines run per genome and their consensus libraries are merged:
//   * Earl Grey (library-only) — RepeatModeler2 + BEAT curation.
//   * HiTE (library-only, --annotate 0) — full-length LTR/TIR/Helitron/non-LTR.
// HiTE runs by default; disable with --skip_hite. CAT_TE_LIBRARIES concatenates
// whatever libraries are present (dedup happens downstream at CD-HIT-EST).
workflow DENOVO_TE_ANIMAL {

    take:
    ch_input        // channel: [ val(meta), path(fasta) ]
    ch_ltrdetector  // value channel: path(LtrDetector dir) for HiTE; [] with --skip_hite

    main:
    ch_versions = Channel.empty()

    EARLGREY(ch_input)
    ch_versions = ch_versions.mix(EARLGREY.out.versions.first())

    ch_libraries = EARLGREY.out.library

    if (!params.skip_hite) {
        HITE(ch_input, ch_ltrdetector)
        ch_versions  = ch_versions.mix(HITE.out.versions.first())
        ch_libraries = ch_libraries.mix(HITE.out.library)
    }

    // Collect all per-sample libraries [meta, lib] into [meta, [libs]] for merging.
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
