include { EARLGREY          } from '../../../modules/local/earlgrey/main'
include { HITE              } from '../../../modules/local/hite/main'
include { GXF_TO_GFF3       } from '../../../modules/local/agat_convert/main'
include { STARFISH_ANNOTATE } from '../../../modules/local/starfish/main'
include { CAT_TE_LIBRARIES  } from '../../../modules/local/cat_te_libraries/main'

// De novo TE discovery for fungal genomes. Per-genome consensus libraries from
// every active engine are merged before Stage 2 dedup:
//   * Earl Grey (library-only) — RepeatModeler2 + BEAT curation. Always runs.
//   * HiTE (library-only, --annotate 0) — full-length LTR/TIR/Helitron/non-LTR.
//     Runs by default; disable with --skip_hite.
//   * Starfish — Starship (giant fungal transposon) detection. OPT-IN
//     (--run_starfish); its heavy Conda solver can be OOM-killed (exit 137),
//     aborting the run. Only runs for genomes that ship a GFF gene annotation.
workflow DENOVO_TE_FUNGUS {

    take:
    ch_input        // channel: [ val(meta), path(fasta), path(gff) | [] ]
    ch_ltrdetector  // value channel: path(LtrDetector dir) for HiTE; [] with --skip_hite

    main:
    ch_versions = Channel.empty()

    // Earl Grey runs on the genome alone, for every fungal sample.
    ch_genome = ch_input.map { meta, fasta, gff -> [ meta, fasta ] }

    EARLGREY(ch_genome)
    ch_versions = ch_versions.mix(EARLGREY.out.versions.first())

    ch_libraries = EARLGREY.out.library

    // HiTE runs by default on the genome alone (organism-agnostic).
    if (!params.skip_hite) {
        HITE(ch_genome, ch_ltrdetector)
        ch_versions  = ch_versions.mix(HITE.out.versions.first())
        ch_libraries = ch_libraries.mix(HITE.out.library)
    }

    // Starfish is opt-in. When disabled, the process is never instantiated, so
    // its Conda env is never built and cannot take down the run.
    ch_starfish_lib = Channel.empty()
    if (params.run_starfish) {
        // Gate Starfish on the presence of a GFF gene annotation.
        ch_input
            .branch { meta, fasta, gff ->
                annotated  : gff          // file present → truthy
                unannotated: true
            }
            .set { ch_starfish }

        ch_starfish.unannotated
            .view { meta, fasta, gff ->
                "WARN | sample=${meta.id}: no gene annotation (gff) provided — " +
                "skipping Starfish (Starship) annotation"
            }

        // Starfish's gene finder needs spec-compliant GFF3. Convert ONLY when the
        // annotation is not already .gff3 (BRAKER/AUGUSTUS GTF, plain .gff);
        // files already in .gff3 skip AGAT and go straight to Starfish.
        ch_starfish.annotated
            .branch { meta, fasta, gff ->
                convert: !gff.toString().toLowerCase().endsWith('.gff3')
                ready  : true
            }
            .set { ch_starfish_ann }

        GXF_TO_GFF3(ch_starfish_ann.convert)
        ch_versions = ch_versions.mix(GXF_TO_GFF3.out.versions.first())

        // Union of freshly-normalized (GXF_TO_GFF3) and already-GFF3 annotations.
        ch_starfish_in = GXF_TO_GFF3.out.annotation.mix(ch_starfish_ann.ready)

        STARFISH_ANNOTATE(ch_starfish_in)
        ch_versions = ch_versions.mix(STARFISH_ANNOTATE.out.versions.first())
        ch_starfish_lib = STARFISH_ANNOTATE.out.library
    }
    ch_libraries = ch_libraries.mix(ch_starfish_lib)

    // Collect all per-sample libraries [meta, lib] into [meta, [libs]] for merging.
    // Earl Grey is always present; HiTE / Starfish contribute when active, and
    // genomes without a Starfish output simply group with the libraries they have.
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
