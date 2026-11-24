include { REPEATMASKER } from '../../../modules/local/repeatmasker/main'

// Stage 3 — Genome masking (single-pass RepeatMasker, soft-masked by default).
//
// One RepeatMasker run over the ORIGINAL genome with the curated + classified library
// and simple/low-complexity detection left ON, so the summary .tbl reports BOTH the
// interspersed TEs AND the tandem/simple/satellite repeats. The run soft-masks every
// masked base (-xsmall) and writes the alignment file (-a) that the visualization
// stage turns into a Kimura divergence landscape.
//
// Deliverables (all published to <outdir>/<id>/masking/):
//   *.softmasked.fa   soft-masked genome
//   *.repeats.tbl     RepeatMasker summary — simple + interspersed
//   *.repeats.align   per-element alignments (Kimura landscape input)
//   *.simple.gff3     simple / low-complexity / satellite repeats only
//   *.te.gff3         interspersed TEs only
//   *.repeats.gff3    combined annotation
workflow MASKING {

    take:
    ch_genome_library  // channel: [ val(meta), path(genome.fa), path(te_library.fa) ]

    main:
    ch_versions = Channel.empty()

    REPEATMASKER(ch_genome_library)
    ch_versions = ch_versions.mix(REPEATMASKER.out.versions.first())

    emit:
    genome_masked = REPEATMASKER.out.masked      // channel: [ val(meta), path(*.softmasked.fa) ]
    repeats_gff   = REPEATMASKER.out.gff          // channel: [ val(meta), path(*.repeats.gff3) ]
    simple_gff    = REPEATMASKER.out.simple_gff   // channel: [ val(meta), path(*.simple.gff3) ]
    te_gff        = REPEATMASKER.out.te_gff        // channel: [ val(meta), path(*.te.gff3) ]
    repeats_bed   = REPEATMASKER.out.bed          // channel: [ val(meta), path(*.repeats.bed) ]
    repeats_out   = REPEATMASKER.out.out          // channel: [ val(meta), path(*.repeats.out) ]
    repeats_tbl   = REPEATMASKER.out.tbl          // channel: [ val(meta), path(*.repeats.tbl) ]
    repeats_align = REPEATMASKER.out.align        // channel: [ val(meta), path(*.repeats.align) ]
    versions      = ch_versions
}
