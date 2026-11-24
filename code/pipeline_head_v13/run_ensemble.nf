#!/usr/bin/env nextflow
// run_ensemble.nf — run the re-classification ENSEMBLE on a ready-made consensus
// FASTA (e.g. a Dfam lineage partition) and emit voting_table.tsv.
//
// Reuses the pipeline's own classifier + ensemble modules and config, but skips
// TETRIMMER / EXTRACT_UNKNOWN so EVERY input sequence is classified and voted.
//
// Usage (run from the pipeline root so nextflow.config is picked up):
//   nextflow run run_ensemble.nf -profile singularity \
//       --input  dfam_partitions/dfam_fungi.fasta \
//       --organism_type fungus \
//       --id dfam_fungi \
//       --outdir ensemble_out
//
// Output: <outdir>/voting_table.tsv   (rows = sequence ids, cols = each model + Ensemble)
//
// Honours the usual switches: --skip_teclass2/-deepte/-create/-neuralte/-terrier,
// the ensemble controls (--ensemble_method two_step|weighted|hard|all, --ensemble_min_frac,
// --ensemble_min_votes, --ensemble_min_weight, --ensemble_weights,
// --ensemble_weights_from), and the *_repo / *_model overrides. With
// --ensemble_method all the voting table gains both Ensemble_hard and
// Ensemble_weighted columns. Sweeping --ensemble_min_frac here is the cheap way to
// pick a gate: it re-votes existing predictions without re-running the classifiers.

nextflow.enable.dsl = 2

include { TECLASS2_CLONE       } from './modules/local/teclass2_clone/main'
include { TECLASS2_FETCH_MODEL } from './modules/local/teclass2_fetch_model/main'
include { TECLASS2             } from './modules/local/teclass2/main'
include { DEEPTE_CLONE         } from './modules/local/deepte_clone/main'
include { DEEPTE_FETCH_MODEL   } from './modules/local/deepte_fetch_model/main'
include { DEEPTE               } from './modules/local/deepte/main'
include { CREATE_CLONE         } from './modules/local/create_clone/main'
include { CREATE_FETCH_MODEL   } from './modules/local/create_fetch_model/main'
include { CREATE               } from './modules/local/create/main'
include { NEURALTE_CLONE       } from './modules/local/neuralte_clone/main'
include { NEURALTE             } from './modules/local/neuralte/main'
include { TERRIER              } from './modules/local/terrier/main'
include { ENSEMBLE_VOTE        } from './modules/local/ensemble_vote/main'
include { checkLibraryIds      } from './modules/local/utils/library_ids'

// DeepTE organism-group selector (mirrors the CURATION subworkflow).
def deepteGroup(organism_type) {
    organism_type == 'plant'  ? 'Plants'
  : organism_type == 'fungus' ? 'Fungi'
  : organism_type == 'animal' ? 'Metazoans'
  :                             'Others'
}

workflow {

    if (!params.input) { error "Provide --input <consensus.fasta>" }

    def meta = [
        id            : params.id ?: file(params.input).getBaseName(),
        organism_type : (params.organism_type ?: 'animal').toLowerCase(),
    ]
    checkLibraryIds(file(params.input, checkIfExists: true))
    ch_input = Channel.value([ meta, file(params.input, checkIfExists: true) ])

    ch_preds = Channel.empty()

    // The classifier CLONE / FETCH_MODEL modules take an ordering-gate input (in the
    // full pipeline it defers them until after TEtrimmer). run_ensemble has no de
    // novo phase, so hand them a ready token to start immediately.
    ch_setup_gate = Channel.value(true)

    // ── TEclass2 ───────────────────────────────────────────────────────────
    if (!params.skip_teclass2) {
        if (params.teclass2_repo) {
            ch_tc_repo = Channel.value(file(params.teclass2_repo, type: 'dir', checkIfExists: true))
        } else {
            TECLASS2_CLONE(ch_setup_gate)
            ch_tc_repo = TECLASS2_CLONE.out.repo
        }
        if (params.teclass2_model) {
            ch_tc_model = Channel.value(file(params.teclass2_model, type: 'dir', checkIfExists: true))
        } else {
            TECLASS2_FETCH_MODEL(ch_setup_gate)
            ch_tc_model = TECLASS2_FETCH_MODEL.out.model
        }
        TECLASS2(ch_input, ch_tc_repo, ch_tc_model)
        ch_preds = ch_preds.mix(TECLASS2.out.predictions)
    }

    // ── DeepTE (organism-aware model) ──────────────────────────────────────
    if (!params.skip_deepte) {
        if (params.deepte_repo) {
            ch_dt_repo = Channel.value(file(params.deepte_repo, type: 'dir', checkIfExists: true))
        } else {
            DEEPTE_CLONE(ch_setup_gate)
            ch_dt_repo = DEEPTE_CLONE.out.repo
        }
        ch_dt_keyed = ch_input.map { m, lib -> [ deepteGroup(m.organism_type), m, lib ] }
        if (params.deepte_model) {
            ch_dt_model = Channel.value(file(params.deepte_model, type: 'dir', checkIfExists: true))
            ch_dt_in    = ch_input.combine(ch_dt_model)
        } else {
            DEEPTE_FETCH_MODEL(ch_dt_keyed.map { g, m, lib -> g }.unique())
            ch_dt_in = ch_dt_keyed
                .combine(DEEPTE_FETCH_MODEL.out.model, by: 0)
                .map { g, m, lib, model -> [ m, lib, model ] }
        }
        DEEPTE(ch_dt_in, ch_dt_repo)
        ch_preds = ch_preds.mix(DEEPTE.out.predictions)
    }

    // ── CREATE ─────────────────────────────────────────────────────────────
    if (!params.skip_create) {
        if (params.create_repo) {
            ch_cr_repo = Channel.value(file(params.create_repo, type: 'dir', checkIfExists: true))
        } else {
            CREATE_CLONE(ch_setup_gate)
            ch_cr_repo = CREATE_CLONE.out.repo
        }
        if (params.create_model) {
            ch_cr_model = Channel.value(file(params.create_model, type: 'dir', checkIfExists: true))
        } else {
            CREATE_FETCH_MODEL(ch_setup_gate)
            ch_cr_model = CREATE_FETCH_MODEL.out.model
        }
        CREATE(ch_input, ch_cr_repo, ch_cr_model)
        ch_preds = ch_preds.mix(CREATE.out.predictions)
    }

    // ── NeuralTE ───────────────────────────────────────────────────────────
    if (!params.skip_neuralte) {
        if (params.neuralte_repo) {
            ch_nte_repo = Channel.value(file(params.neuralte_repo, type: 'dir', checkIfExists: true))
        } else {
            NEURALTE_CLONE(ch_setup_gate)
            ch_nte_repo = NEURALTE_CLONE.out.repo
        }
        NEURALTE(ch_input, ch_nte_repo)
        ch_preds = ch_preds.mix(NEURALTE.out.predictions)
    }

    // ── Terrier (weights baked into the image) ─────────────────────────────
    if (!params.skip_terrier) {
        TERRIER(ch_input)
        ch_preds = ch_preds.mix(TERRIER.out.predictions)
    }

    // ── Vote: group each sample's per-model preds and majority-vote ────────
    ch_input
        .join(ch_preds.groupTuple(), remainder: true)
        .branch { m, lib, preds ->
            voted : preds
            novote: true
        }
        .set { ch_split }

    ch_split.novote.view { m, lib, preds ->
        "WARN | ${m.id}: no classifier produced predictions — no voting table" }

    ENSEMBLE_VOTE(
        ch_split.voted,
        params.ensemble_two_step_config ? file(params.ensemble_two_step_config, checkIfExists: true) : [],
        params.ensemble_weights_from    ? file(params.ensemble_weights_from,    checkIfExists: true) : []
    )

    // Land the table at <outdir>/voting_table.tsv exactly.
    ENSEMBLE_VOTE.out.table
        .map { m, tbl -> tbl }
        .collectFile(name: 'voting_table.tsv', storeDir: params.outdir)
}
