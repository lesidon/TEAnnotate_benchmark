include { EDTA             } from '../../../modules/local/edta/main'
include { HELIANO          } from '../../../modules/local/heliano/main'
include { HITE             } from '../../../modules/local/hite/main'
include { CAT_TE_LIBRARIES } from '../../../modules/local/cat_te_libraries/main'

// De novo TE discovery for plant genomes. Per-genome consensus libraries from
// every engine are merged before Stage 2 dedup:
//   * EDTA (library-only, --anno 0, --species others) — the plant-standard builder
//     (LTR + TIR-Learner + HelitronScanner + GRF), in place of Earl Grey.
//   * HELIANO — dedicated, more accurate Helitron detector than EDTA's bundled one.
//   * HiTE (library-only) — full-length LTR/TIR/Helitron/non-LTR; runs by default
//     across all organisms (disable with --skip_hite).
// CAT_TE_LIBRARIES concatenates whatever libraries are present; CD-HIT-EST collapses
// the (expected) overlap downstream.
workflow DENOVO_TE_PLANT {

    take:
    ch_input        // channel: [ val(meta), path(fasta) ]
    ch_ltrdetector  // value channel: path(LtrDetector dir) for HiTE; [] with --skip_hite

    main:
    ch_versions = Channel.empty()

    EDTA(ch_input)
    ch_versions  = ch_versions.mix(EDTA.out.versions.first())
    ch_libraries = EDTA.out.library

    HELIANO(ch_input)
    ch_versions  = ch_versions.mix(HELIANO.out.versions.first())
    ch_libraries = ch_libraries.mix(HELIANO.out.library)

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
