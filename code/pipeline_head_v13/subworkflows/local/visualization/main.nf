include { CALC_DIVERGENCE } from '../../../modules/local/calc_divergence/main'
include { PLOT_KIMURA     } from '../../../modules/local/plot_kimura/main'
include { PLOT_KARYOTYPE  } from '../../../modules/local/plot_karyotype/main'
include { REPEAT_QC       } from '../../../modules/local/repeat_qc/main'

// Stage 4 — Visualization & QC (Nature-styled). Three deliverables, each shipping the
// exact R script that drew it (published next to the figure so the user can restyle):
//   KIMURA_LANDSCAPE  TE divergence landscape from the RepeatMasker .align
//   KARYOPLOT         per-scaffold TE-density chromosome painting
//   REPEAT_QC         repeat detection & masking QC (composition + PASS/WARN metrics)
// The R plotting scripts live in bin/ and are staged into each process.
workflow VISUALIZATION {

    take:
    ch_align    // channel: [ val(meta), path(repeats.align) ]
    ch_out      // channel: [ val(meta), path(repeats.out) ]
    ch_genome   // channel: [ val(meta), path(genome.fa) ]
    ch_bed      // channel: [ val(meta), path(repeats.bed) ]
    ch_tbl      // channel: [ val(meta), path(repeats.tbl) ]
    ch_library  // channel: [ val(meta), path(classified_library.fa) ]

    main:
    ch_versions = Channel.empty()

    rscript_kimura = file("${projectDir}/bin/plot_kimura.R")
    rscript_karyo  = file("${projectDir}/bin/plot_karyotype.R")
    rscript_qc     = file("${projectDir}/bin/repeat_qc.R")

    // ── Kimura divergence landscape (join .out + .align by sample id) ──────────
    ch_out
        .map { meta, o -> [ meta.id, meta, o ] }
        .join(ch_align.map { meta, a -> [ meta.id, a ] }, by: 0)
        .map { id, meta, o, a -> [ meta, o, a ] }
        .set { ch_div_in }
    CALC_DIVERGENCE(ch_div_in)
    ch_versions = ch_versions.mix(CALC_DIVERGENCE.out.versions.first())

    // the .tbl carries the genome length for the % of genome axis
    CALC_DIVERGENCE.out.matrix
        .map { meta, m -> [ meta.id, meta, m ] }
        .join(ch_tbl.map { meta, tbl -> [ meta.id, tbl ] }, by: 0)
        .map { id, meta, m, tbl -> [ meta, m, tbl ] }
        .set { ch_kimura_in }
    PLOT_KIMURA(ch_kimura_in, rscript_kimura)
    ch_versions = ch_versions.mix(PLOT_KIMURA.out.versions.first())

    // ── Karyoplot: join genome + combined BED by sample id ────────────────────
    ch_bed
        .map { meta, bed -> [ meta.id, meta, bed ] }
        .join(ch_genome.map { meta, g -> [ meta.id, g ] }, by: 0)
        .map { id, meta, bed, g -> [ meta, g, bed ] }
        .set { ch_karyo_in }
    PLOT_KARYOTYPE(ch_karyo_in, rscript_karyo)
    ch_versions = ch_versions.mix(PLOT_KARYOTYPE.out.versions.first())

    // ── QC: join .tbl + combined BED + classified library by sample id ────────
    ch_tbl
        .map { meta, tbl -> [ meta.id, meta, tbl ] }
        .join(ch_bed.map     { meta, bed -> [ meta.id, bed ] }, by: 0)
        .join(ch_library.map { meta, lib -> [ meta.id, lib ] }, by: 0)
        .map { id, meta, tbl, bed, lib -> [ meta, tbl, bed, lib ] }
        .set { ch_qc_in }
    REPEAT_QC(ch_qc_in, rscript_qc)
    ch_versions = ch_versions.mix(REPEAT_QC.out.versions.first())

    emit:
    kimura   = PLOT_KIMURA.out.plot
    karyo    = PLOT_KARYOTYPE.out.plot
    qc       = REPEAT_QC.out.plot
    versions = ch_versions
}
