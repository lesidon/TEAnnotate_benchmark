include { TETRIMMER            } from '../../../modules/local/tetrimmer/main'
include { TETRIMMER_POSTPROCESS } from '../../../modules/local/tetrimmer_postprocess/main'
include { CD_HIT_EST as CD_HIT_EST_CURATED     } from '../../../modules/local/cd_hit_est/main'
include { DEDUP_REFINE as DEDUP_REFINE_CURATED } from '../../../modules/local/dedup_refine/main'
include { DOMAIN_RELABEL       } from '../../../modules/local/domain_relabel/main'
include { EXTRACT_UNKNOWN      } from '../../../modules/local/extract_unknown/main'
include { TECLASS2_CLONE       } from '../../../modules/local/teclass2_clone/main'
include { TECLASS2_FETCH_MODEL } from '../../../modules/local/teclass2_fetch_model/main'
include { TECLASS2             } from '../../../modules/local/teclass2/main'
include { DEEPTE_CLONE         } from '../../../modules/local/deepte_clone/main'
include { DEEPTE_FETCH_MODEL   } from '../../../modules/local/deepte_fetch_model/main'
include { DEEPTE               } from '../../../modules/local/deepte/main'
include { CREATE_CLONE         } from '../../../modules/local/create_clone/main'
include { CREATE_FETCH_MODEL   } from '../../../modules/local/create_fetch_model/main'
include { CREATE               } from '../../../modules/local/create/main'
include { NEURALTE_CLONE       } from '../../../modules/local/neuralte_clone/main'
include { NEURALTE             } from '../../../modules/local/neuralte/main'
include { TERRIER              } from '../../../modules/local/terrier/main'
include { ENSEMBLE_VOTE        } from '../../../modules/local/ensemble_vote/main'

// Map a sample's organism_type to DeepTE's pre-trained model group (the tarball /
// <Group>_model subdir name). Mirrors the -sp selection inside the DEEPTE module.
def deepteGroup(organism_type) {
    organism_type == 'plant'  ? 'Plants'
  : organism_type == 'fungus' ? 'Fungi'
  : organism_type == 'animal' ? 'Metazoans'
  :                             'Others'
}

// Kingdom-aware model SELECTION — which classifiers actually run for a sample.
// With ensemble_kingdom_aware on (default), only the subset per organism_type
// (ensemble_models_{plant,animal,fungus}; by default the models of the two-step
// settings) is launched; the other tools are skipped entirely rather than run and
// then dropped from the vote. Returns null => no restriction (run all).
def kingdomModels(organism_type) {
    if (!params.ensemble_kingdom_aware) return null
    def sel = organism_type == 'plant'  ? params.ensemble_models_plant
            : organism_type == 'fungus' ? params.ensemble_models_fungus
            : organism_type == 'animal' ? params.ensemble_models_animal
            :                             null
    if (!sel) return null
    return ((sel.toString().toLowerCase().split(/[,;\s]+/) as List).findAll { it }) as Set
}

// True if classifier `tool` (lowercase name) should run for this sample.
def runsFor(meta, tool) {
    def models = kingdomModels(meta.organism_type)
    models == null ? true : models.contains(tool)
}

// Stage 3 — Curation, Classification & Masking.
// TETRIMMER curates the deduplicated consensus library (carrying Earl Grey /
// RepeatClassifier labels, with an "Unknown" fraction). That Unknown fraction is
// re-classified by a voting ENSEMBLE of up to five models (TEclass2, DeepTE, CREATE,
// NeuralTE, Terrier), and ENSEMBLE_VOTE assigns the class (by default a two-step
// vote: Order, then a superfamily inside it). With ensemble_kingdom_aware (default),
// only the subset per organism_type actually RUNS — e.g. a fungal sample uses CREATE,
// NeuralTE and Terrier (see kingdomModels/runsFor) — so the run is not loaded with
// classifiers that would be dropped from the vote anyway. The step runs
// BY DEFAULT; skip the whole ensemble with --skip_reclassification, drop a member with
// --skip_teclass2 / --skip_deepte / --skip_create / --skip_neuralte / --skip_terrier,
// or set --ensemble_kingdom_aware false to run every enabled model everywhere. Each
// member is fetched + cached at run time and is resilient (errorStrategy 'ignore' ->
// abstains). REPEATMASKER (genome masking) follows.
workflow CURATION {

    take:
    ch_library_genome  // channel: [ val(meta), path(dedup_library.fa), path(genome.fa) ]
    ch_evidence        // channel: [ val(meta), path(evidence.tsv) ] from TE_EVIDENCE; may be empty
    ch_lineage         // channel: [ val(meta), path(lineage.tsv) | [] ] from LIBRARY_DEDUP; may be empty
    ch_pfam            // value channel: Pfam-A dir (mandatory for TEtrimmer); unused with --skip_tetrimmer

    main:
    ch_versions  = Channel.empty()
    ch_manifests = Channel.empty()   // source.tsv of every asset fetched here

    if (params.skip_tetrimmer) {
        // Ablation (--skip_tetrimmer): no TEtrimmer curation. The deduplicated de novo
        // library goes straight to re-classification, so the comparison with the
        // default run isolates what TEtrimmer adds.
        ch_curated = ch_library_genome.map { meta, lib, genome -> [ meta, lib ] }
    } else {
        // PFAM is mandatory for TEtrimmer (it orients TEs; curation aborts without it).
        // Fetched once by the calling workflow and shared with TE_EVIDENCE.
        TETRIMMER(ch_library_genome, ch_pfam)
        ch_versions = ch_versions.mix(TETRIMMER.out.versions.first())

        if (params.curation_class_aware) {
            // TEtrimmer's own final merge is length-based, like CD-HIT's: a consensus
            // of another class collapses into the longest cluster member. Take its
            // un-merged output instead, re-check labels of consensi TEtrimmer rebuilt
            // from another element, and merge with the class-aware Stage 2 steps.
            // With the evidence table, inputs TEtrimmer dropped (fewer than 10 genome
            // hits and fewer than 2 full-length copies) come back unchanged when they
            // carry a TE domain or another tool's support (--tetrimmer_reinject).
            // With the lineage, the user's records come back as given and are merged as
            // user records again, so the curated merge keeps them as representatives.
            TETRIMMER.out.output_dir
                .join(ch_library_genome.map { meta, lib, genome -> [ meta, lib ] })
                .join(ch_evidence, remainder: true)
                .filter { it[1] != null }
                .join(ch_lineage, remainder: true)
                .filter { it[1] != null }
                .map { meta, dir, lib, evidence, lineage -> [ meta, dir, lib, evidence ?: [], lineage ?: [] ] }
                .set { ch_tt_out }
            TETRIMMER_POSTPROCESS(ch_tt_out)
            CD_HIT_EST_CURATED(
                TETRIMMER_POSTPROCESS.out.library
                    .join(TETRIMMER_POSTPROCESS.out.user)
                    .map { meta, lib, user -> [ meta, [ lib ], user ] }
            )
            // The genome lets the refinement pick each cluster's representative by
            // full-length genomic copies rather than length (--dedup_rep_by_copies).
            DEDUP_REFINE_CURATED(
                CD_HIT_EST_CURATED.out.clustered
                    .join(ch_library_genome.map { meta, lib, genome -> [ meta, genome ] })
            )
            ch_versions = ch_versions
                .mix(TETRIMMER_POSTPROCESS.out.versions.first())
                .mix(CD_HIT_EST_CURATED.out.versions.first())
                .mix(DEDUP_REFINE_CURATED.out.versions.first())
            ch_curated         = DEDUP_REFINE_CURATED.out.library
            ch_curated_lineage = DEDUP_REFINE_CURATED.out.lineage
        } else {
            ch_curated         = TETRIMMER.out.library
            ch_curated_lineage = TETRIMMER.out.library.map { meta, lib -> [ meta, [] ] }
        }

        // A curated consensus is one element: when its TE domains point to another
        // class than its label (e.g. a HELIANO "Helitron" with a Gypsy RT), the domain
        // class is used. Not in the --skip_tetrimmer arm, where chimeras remain. The
        // user's records (from the lineage) keep their labels.
        if (params.curation_domain_relabel) {
            DOMAIN_RELABEL(ch_curated.join(ch_curated_lineage))
            ch_versions = ch_versions.mix(DOMAIN_RELABEL.out.versions.first())
            ch_curated  = DOMAIN_RELABEL.out.library
        }
    }
    ch_preds   = Channel.empty()
    def ran_any = false

    if (!params.skip_reclassification) {

        // Re-classify ONLY the Unknown fraction; classified consensi keep their
        // curated labels. Only the kingdom-selected classifiers run on this subset.
        EXTRACT_UNKNOWN(ch_curated)
        ch_versions = ch_versions.mix(EXTRACT_UNKNOWN.out.versions.first())
        ch_unknown = EXTRACT_UNKNOWN.out.unknown

        // Kingdom-aware EXECUTION: route each sample's Unknown library only to the
        // classifiers selected for its organism_type (see kingdomModels/runsFor).
        // A tool no sample selects gets an empty channel and is never launched, so a
        // fungus run (subset = Terrier) does not also spin up TEclass2/DeepTE/etc.
        ch_unk_teclass2 = ch_unknown.filter { meta, lib -> runsFor(meta, 'teclass2') }
        ch_unk_deepte   = ch_unknown.filter { meta, lib -> runsFor(meta, 'deepte')   }
        ch_unk_create   = ch_unknown.filter { meta, lib -> runsFor(meta, 'create')   }
        ch_unk_neuralte = ch_unknown.filter { meta, lib -> runsFor(meta, 'neuralte') }
        ch_unk_terrier  = ch_unknown.filter { meta, lib -> runsFor(meta, 'terrier')  }

        // Per-tool ordering gate: a single-emit trigger that fires only when >=1
        // sample selects the tool, so its no-data-input CLONE / FETCH_MODEL run AFTER
        // TEtrimmer (not at t=0) AND are skipped for tools no sample needs.
        ch_gate_teclass2 = ch_unk_teclass2.map { meta, lib -> true }.first()
        ch_gate_deepte   = ch_unk_deepte.map   { meta, lib -> true }.first()
        ch_gate_create   = ch_unk_create.map   { meta, lib -> true }.first()
        ch_gate_neuralte = ch_unk_neuralte.map { meta, lib -> true }.first()

        // ── TEclass2 (Transformer) ────────────────────────────────────────────
        if (!params.skip_teclass2) {
            if (params.teclass2_repo) {
                ch_tc_repo = Channel.value(file(params.teclass2_repo, type: 'dir', checkIfExists: true))
            } else {
                TECLASS2_CLONE(ch_gate_teclass2)
                ch_tc_repo  = TECLASS2_CLONE.out.repo
                ch_versions = ch_versions.mix(TECLASS2_CLONE.out.versions)
                ch_manifests = ch_manifests.mix(TECLASS2_CLONE.out.manifest)
            }
            if (params.teclass2_model) {
                ch_tc_model = Channel.value(file(params.teclass2_model, type: 'dir', checkIfExists: true))
            } else {
                TECLASS2_FETCH_MODEL(ch_gate_teclass2)
                ch_tc_model = TECLASS2_FETCH_MODEL.out.model
                ch_versions = ch_versions.mix(TECLASS2_FETCH_MODEL.out.versions)
                ch_manifests = ch_manifests.mix(TECLASS2_FETCH_MODEL.out.manifest)
            }
            TECLASS2(ch_unk_teclass2, ch_tc_repo, ch_tc_model)
            ch_versions = ch_versions.mix(TECLASS2.out.versions.first())
            ch_preds = ch_preds.mix(TECLASS2.out.predictions)
            ran_any = true
        }

        // ── DeepTE (CNN, organism-aware) ──────────────────────────────────────
        if (!params.skip_deepte) {
            if (params.deepte_repo) {
                ch_dt_repo = Channel.value(file(params.deepte_repo, type: 'dir', checkIfExists: true))
            } else {
                DEEPTE_CLONE(ch_gate_deepte)
                ch_dt_repo  = DEEPTE_CLONE.out.repo
                ch_versions = ch_versions.mix(DEEPTE_CLONE.out.versions)
                ch_manifests = ch_manifests.mix(DEEPTE_CLONE.out.manifest)
            }

            // Tag each sample with its DeepTE organism group, then attach the
            // matching model. Only the group(s) actually present are fetched, once
            // each (cached per group) — not the whole multi-GB suite.
            ch_dt_keyed = ch_unk_deepte.map { meta, lib -> [ deepteGroup(meta.organism_type), meta, lib ] }

            if (params.deepte_model) {
                // User-supplied dir holding all <Group>_model subdirs: broadcast it.
                ch_dt_model = Channel.value(file(params.deepte_model, type: 'dir', checkIfExists: true))
                ch_dt_in    = ch_unk_deepte.combine(ch_dt_model)
            } else {
                DEEPTE_FETCH_MODEL(ch_dt_keyed.map { group, meta, lib -> group }.unique())
                ch_versions = ch_versions.mix(DEEPTE_FETCH_MODEL.out.versions)
                ch_manifests = ch_manifests.mix(DEEPTE_FETCH_MODEL.out.manifest)
                // Join each sample to its group's fetched model by the group key.
                ch_dt_in = ch_dt_keyed
                    .combine(DEEPTE_FETCH_MODEL.out.model, by: 0)
                    .map { group, meta, lib, model -> [ meta, lib, model ] }
            }
            DEEPTE(ch_dt_in, ch_dt_repo)
            ch_versions = ch_versions.mix(DEEPTE.out.versions.first())
            ch_preds = ch_preds.mix(DEEPTE.out.predictions)
            ran_any = true
        }

        // ── CREATE (attention CNN+RNN; Wicker taxonomy) ───────────────────────
        if (!params.skip_create) {
            if (params.create_repo) {
                ch_cr_repo = Channel.value(file(params.create_repo, type: 'dir', checkIfExists: true))
            } else {
                CREATE_CLONE(ch_gate_create)
                ch_cr_repo  = CREATE_CLONE.out.repo
                ch_versions = ch_versions.mix(CREATE_CLONE.out.versions)
                ch_manifests = ch_manifests.mix(CREATE_CLONE.out.manifest)
            }
            if (params.create_model) {
                ch_cr_model = Channel.value(file(params.create_model, type: 'dir', checkIfExists: true))
            } else {
                CREATE_FETCH_MODEL(ch_gate_create)
                ch_cr_model = CREATE_FETCH_MODEL.out.model
                ch_versions = ch_versions.mix(CREATE_FETCH_MODEL.out.versions)
                ch_manifests = ch_manifests.mix(CREATE_FETCH_MODEL.out.manifest)
            }
            CREATE(ch_unk_create, ch_cr_repo, ch_cr_model)
            ch_versions = ch_versions.mix(CREATE.out.versions.first())
            ch_preds = ch_preds.mix(CREATE.out.predictions)
            ran_any = true
        }

        // ── NeuralTE (CNN; RepeatMasker labels via --is_wicker 0) ─────────────
        if (!params.skip_neuralte) {
            if (params.neuralte_repo) {
                ch_nte_repo = Channel.value(file(params.neuralte_repo, type: 'dir', checkIfExists: true))
            } else {
                NEURALTE_CLONE(ch_gate_neuralte)
                ch_nte_repo = NEURALTE_CLONE.out.repo
                ch_versions = ch_versions.mix(NEURALTE_CLONE.out.versions)
                ch_manifests = ch_manifests.mix(NEURALTE_CLONE.out.manifest)
            }
            NEURALTE(ch_unk_neuralte, ch_nte_repo)
            ch_versions = ch_versions.mix(NEURALTE.out.versions.first())
            ch_preds = ch_preds.mix(NEURALTE.out.predictions)
            ran_any = true
        }

        // ── Terrier (PyTorch/corgi; RepeatMasker labels; weights baked in image) ─
        if (!params.skip_terrier) {
            TERRIER(ch_unk_terrier)
            ch_versions = ch_versions.mix(TERRIER.out.versions.first())
            ch_preds = ch_preds.mix(TERRIER.out.predictions)
            ran_any = true
        }
    }

    if (ran_any) {
        // Gather each sample's per-model vote TSVs, then majority-vote. Samples
        // with no votes (all members skipped/failed) pass the curated library
        // through unchanged.
        ch_curated
            .join(ch_preds.groupTuple(), remainder: true)
            .branch { meta, lib, preds ->
                voted  : preds
                novote : true
            }
            .set { ch_split }

        ENSEMBLE_VOTE(
            ch_split.voted,
            params.ensemble_two_step_config ? file(params.ensemble_two_step_config, checkIfExists: true) : [],
            params.ensemble_weights_from    ? file(params.ensemble_weights_from,    checkIfExists: true) : []
        )
        ch_versions = ch_versions.mix(ENSEMBLE_VOTE.out.versions.first())

        ch_library = ENSEMBLE_VOTE.out.library
            .mix(ch_split.novote.map { meta, lib, preds -> [ meta, lib ] })
    }
    else {
        ch_curated.view { meta, lib ->
            "INFO | sample=${meta.id}: re-classification skipped (--skip_reclassification) " +
            "— keeping curated labels"
        }
        ch_library = ch_curated
    }

    emit:
    library   = ch_library    // channel: [ val(meta), path(*.{classified,curated}_library.fa) ]
    manifests = ch_manifests  // channel: path(source.tsv), one per fetched asset
    versions  = ch_versions
}
