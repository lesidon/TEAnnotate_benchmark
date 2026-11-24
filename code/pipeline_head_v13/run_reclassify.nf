#!/usr/bin/env nextflow
// run_reclassify.nf — run the pipeline's library stages on a ready-made library,
// without re-running de novo discovery:
//   [--dedup]  Stage 2 LIBRARY_DEDUP (CD-HIT-EST + class-aware DEDUP_REFINE) on a raw
//              library, e.g. a finished run's 02_denovo_te/<id>.raw_library.fa
//   gate       TE_EVIDENCE (unless --evidence_gate off): keeps consensi with a TE
//              domain, another tool's support or >= 2 full-length copies; needs
//              --genome. Tools are read from --provenance (name<TAB>tool, written by
//              CAT_TE_LIBRARIES) or, without it, from each tool's naming scheme.
//   always     Stage 3 CURATION: TEtrimmer (unless --skip_tetrimmer), then the
//              kingdom-aware ensemble re-classifies the Unknown fraction and
//              ENSEMBLE_VOTE relabels it, exactly as in the full pipeline.
// With --skip_tetrimmer this is the "no TEtrimmer" ablation.
//
// Usage (run from the pipeline root so nextflow.config is picked up):
//   nextflow run run_reclassify.nf -profile singularity \
//       --input out/<sample>/02_denovo_te/<sample>.raw_library.fa --dedup \
//       --genome <genome.fa> --id <sample> --organism_type plant|animal|fungus \
//       [--provenance out/<sample>/02_denovo_te/<sample>.provenance.tsv] \
//       [--skip_tetrimmer] --outdir reclassify_out
//
// Output: <outdir>/<id>/05_classification/<id>.te_library.fa, the final library on
// every branch (re-classified, curated, or deduplicated when steps are skipped or every
// classifier fails), as in the full pipeline; with --dedup also <outdir>/<id>/03_library/.
// Published by the same modules.config rules as the full pipeline. The genome passes
// the same INPUT_CHECK (gzip -> plain FASTA; CRLF and duplicate IDs rejected).

nextflow.enable.dsl = 2

include { INPUT_CHECK   } from './modules/local/input_check/main'
include { LIBRARY_DEDUP } from './subworkflows/local/library_dedup/main'
include { TE_EVIDENCE   } from './modules/local/te_evidence/main'
include { PFAM_FETCH    } from './modules/local/pfam_fetch/main'
include { CURATION      } from './subworkflows/local/curation/main'
include { FINAL_LIBRARY } from './modules/local/final_library/main'
include { checkLibraryIds } from './modules/local/utils/library_ids'
include { suppliedAssets  } from './modules/local/utils/asset_key'

workflow {

    if (!params.input) { error "Provide --input <library.fasta>" }
    def need_genome = !params.skip_tetrimmer || params.evidence_gate != 'off'
    if (need_genome && !params.genome) {
        error "TEtrimmer and the evidence gate need the genome: pass --genome <genome.fa> " +
              "(or --skip_tetrimmer --evidence_gate off)"
    }

    def meta = [
        id            : params.id ?: file(params.input).getBaseName(),
        organism_type : (params.organism_type ?: 'animal').toLowerCase(),
    ]
    def raw = file(params.input, checkIfExists: true)
    // With --dedup, CD_HIT_EST checks the IDs (and renames user collisions).
    if (!params.dedup) checkLibraryIds(raw)
    ch_library = Channel.of([ meta, raw ])
    // Without --dedup the library is taken as is: no lineage, so no sequence is the user's.
    ch_lineage = Channel.of([ meta, [] ])

    ch_versions  = Channel.empty()
    ch_manifests = Channel.empty()

    // The validated, plain-FASTA genome: [ meta, genome ], or nothing without --genome.
    ch_genome = Channel.empty()
    if (params.genome) {
        INPUT_CHECK(Channel.of([ meta, file(params.genome, checkIfExists: true) ]))
        ch_genome   = INPUT_CHECK.out.genome
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    }

    if (params.dedup) {
        // with the genome, cluster representatives are chosen by full-length copies
        LIBRARY_DEDUP(ch_library, ch_genome.map { m, g -> [ m.id, g ] })
        ch_library  = LIBRARY_DEDUP.out.library
        ch_lineage  = LIBRARY_DEDUP.out.lineage
        ch_versions = ch_versions.mix(LIBRARY_DEDUP.out.versions)
    }

    // Pfam-A for TEtrimmer and for TE_EVIDENCE's domain search (fetched once, cached)
    def need_pfam = !params.skip_tetrimmer || (params.evidence_gate != 'off' && params.evidence_pfam)
    ch_pfam = Channel.value([])
    if (need_pfam) {
        if (params.pfam_db) {
            ch_pfam = Channel.value(file(params.pfam_db, type: 'dir', checkIfExists: true))
        } else {
            PFAM_FETCH()
            ch_pfam      = PFAM_FETCH.out.db.first()
            ch_versions  = ch_versions.mix(PFAM_FETCH.out.versions)
            ch_manifests = ch_manifests.mix(PFAM_FETCH.out.manifest)
        }
    }

    ch_evidence = Channel.empty()
    if (params.evidence_gate != 'off') {
        def prov    = params.containsKey('provenance') && params.provenance
                        ? file(params.provenance, checkIfExists: true) : []
        TE_EVIDENCE(
            ch_library.join(ch_lineage).join(ch_genome).map { m, lib, lin, g ->
                [ m, lib, raw, g, prov, lin ] },
            ch_pfam
        )
        ch_library  = TE_EVIDENCE.out.library
        ch_evidence = TE_EVIDENCE.out.evidence
        ch_versions = ch_versions.mix(TE_EVIDENCE.out.versions)
    }

    // With --skip_tetrimmer and no gate, the library stands in for the genome so the
    // CURATION input shape is unchanged (TEtrimmer is its only consumer there).
    ch_curation_input = params.genome
        ? ch_library.join(ch_genome)
        : ch_library.map { m, lib -> [ m, lib, lib ] }

    CURATION(ch_curation_input, ch_evidence, ch_lineage, ch_pfam)

    // The library the full pipeline would mask with, published as <id>.te_library.fa.
    FINAL_LIBRARY(CURATION.out.library)

    // Software versions and assets of the run, as in the full pipeline.
    ch_versions
        .mix(CURATION.out.versions, FINAL_LIBRARY.out.versions)
        .map { it.text.trim() }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info",
                     newLine: true, sort: true,
                     seed: "\"Workflow\":\n    nextflow: ${workflow.nextflow.version}\n" +
                           "    teannotate: run_reclassify.nf (${workflow.commitId ?: 'local checkout'})")
    ch_manifests
        .mix(CURATION.out.manifests)
        .map { f -> f.readLines().findAll { it }.collect { "${f.parent}\t${it}" }.join('\n') }
        .mix(Channel.fromList(suppliedAssets()))
        .collectFile(name: 'asset_manifest.tsv', storeDir: "${params.outdir}/pipeline_info",
                     newLine: true, sort: true, seed: 'asset\tfield\tvalue')
}
