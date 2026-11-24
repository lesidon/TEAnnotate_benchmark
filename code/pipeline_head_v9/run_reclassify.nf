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
// Output: <outdir>/<id>/05_classification/<id>.classified_library.fa (+ voting table);
// with --dedup also <outdir>/<id>/03_library/. Published by the same modules.config
// rules as the full pipeline.

nextflow.enable.dsl = 2

include { LIBRARY_DEDUP } from './subworkflows/local/library_dedup/main'
include { TE_EVIDENCE   } from './modules/local/te_evidence/main'
include { PFAM_FETCH    } from './modules/local/pfam_fetch/main'
include { CURATION      } from './subworkflows/local/curation/main'

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
    ch_library = Channel.of([ meta, raw ])

    if (params.dedup) {
        // with the genome, cluster representatives are chosen by full-length copies
        ch_genome_by_id = params.genome
            ? Channel.of([ meta.id, file(params.genome, checkIfExists: true) ])
            : Channel.empty()
        LIBRARY_DEDUP(ch_library, ch_genome_by_id)
        ch_library = LIBRARY_DEDUP.out.library
    }

    // Pfam-A for TEtrimmer and for TE_EVIDENCE's domain search (fetched once, cached)
    def need_pfam = !params.skip_tetrimmer || (params.evidence_gate != 'off' && params.evidence_pfam)
    ch_pfam = Channel.value([])
    if (need_pfam) {
        if (params.pfam_db) {
            ch_pfam = Channel.value(file(params.pfam_db, type: 'dir', checkIfExists: true))
        } else {
            PFAM_FETCH()
            ch_pfam = PFAM_FETCH.out.db.first()
        }
    }

    ch_evidence = Channel.empty()
    if (params.evidence_gate != 'off') {
        def prov    = params.containsKey('provenance') && params.provenance
                        ? file(params.provenance, checkIfExists: true) : []
        def trusted = params.user_te_library ? file(params.user_te_library, checkIfExists: true) : []
        TE_EVIDENCE(
            ch_library.map { m, lib -> [ m, lib, raw, file(params.genome, checkIfExists: true), prov ] },
            trusted,
            ch_pfam
        )
        ch_library  = TE_EVIDENCE.out.library
        ch_evidence = TE_EVIDENCE.out.evidence
    }

    // With --skip_tetrimmer and no gate, the library stands in for the genome so the
    // CURATION input shape is unchanged (TEtrimmer is its only consumer there).
    ch_curation_input = ch_library.map { m, lib ->
        [ m, lib, params.genome ? file(params.genome, checkIfExists: true) : lib ]
    }

    CURATION(ch_curation_input, ch_evidence, ch_pfam)
}
