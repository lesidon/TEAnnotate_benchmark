include { REPEATMASKER } from '../../../modules/local/repeatmasker/main'
include { RIP_SCAN     } from '../../../modules/local/rip_scan/main'

// Stage 3 — Genome masking (single-pass RepeatMasker, soft-masked by default).
//
// One RepeatMasker run over the ORIGINAL genome with the curated + classified library
// and simple/low-complexity detection left ON, so the summary .tbl reports BOTH the
// interspersed TEs AND the tandem/simple/satellite repeats. The run soft-masks every
// masked base (-xsmall) and writes the alignment file (-a) that the visualization
// stage turns into a Kimura divergence landscape.
//
// Deliverables (all published to <outdir>/<id>/06_masking/):
//   *.softmasked.fa   soft-masked genome
//   *.repeats.tbl     RepeatMasker summary — simple + interspersed
//   *.repeats.align   per-element alignments (Kimura landscape input)
//   *.simple.gff3     simple / low-complexity / satellite repeats only
//   *.te.gff3         interspersed TEs only
//   *.repeats.gff3    combined annotation
//   *.rip_*           fungal genomes only (unless --skip_rip): RIP mutations per TE
//                     copy and family, RIP indices, large RIP-affected regions
workflow MASKING {

    take:
    ch_genome_library  // channel: [ val(meta), path(genome.fa), path(te_library.fa) ]

    main:
    ch_versions = Channel.empty()

    REPEATMASKER(ch_genome_library)
    ch_versions = ch_versions.mix(REPEATMASKER.out.versions.first())

    // RIP is a fungal genome defence, so only the fungal route is scanned.
    ch_rip = Channel.empty()
    if (!params.skip_rip) {
        ch_genome_library
            .filter { meta, genome, lib -> meta.organism_type == 'fungus' }
            .map { meta, genome, lib -> [ meta.id, meta, genome ] }
            .join(REPEATMASKER.out.align.map { meta, align -> [ meta.id, align ] }, by: 0)
            .map { id, meta, genome, align -> [ meta, genome, align ] }
            .set { ch_rip_in }
        RIP_SCAN(ch_rip_in)
        ch_versions = ch_versions.mix(RIP_SCAN.out.versions.first())
        ch_rip      = RIP_SCAN.out.report
    }

    emit:
    genome_masked = REPEATMASKER.out.masked      // channel: [ val(meta), path(*.softmasked.fa) ]
    repeats_gff   = REPEATMASKER.out.gff          // channel: [ val(meta), path(*.repeats.gff3) ]
    simple_gff    = REPEATMASKER.out.simple_gff   // channel: [ val(meta), path(*.simple.gff3) ]
    te_gff        = REPEATMASKER.out.te_gff        // channel: [ val(meta), path(*.te.gff3) ]
    repeats_bed   = REPEATMASKER.out.bed          // channel: [ val(meta), path(*.repeats.bed) ]
    repeats_out   = REPEATMASKER.out.out          // channel: [ val(meta), path(*.repeats.out) ]
    repeats_tbl   = REPEATMASKER.out.tbl          // channel: [ val(meta), path(*.repeats.tbl) ]
    repeats_align = REPEATMASKER.out.align        // channel: [ val(meta), path(*.repeats.align) ]
    rip           = ch_rip                        // channel: [ val(meta), [ path(*.rip_*) ] ] (fungi)
    versions      = ch_versions
}
