// Ensemble for TE re-classification (Stage 3).
// Each enabled classifier (TEclass2, DeepTE, CREATE, NeuralTE, Terrier) contributes
// a normalized [base_id, RM_label] vote TSV. For every consensus still labelled
// Unknown, the class is chosen by params.ensemble_method:
//   two_step (default) — Order vote with per-kingdom weights and gate, then a
//                        superfamily vote INSIDE the winning Order with its own
//                        per-kingdom weights (bin/ensemble_vote.py TWO_STEP_DEFAULTS,
//                        overridable with params.ensemble_two_step_config).
//   weighted           — each model's vote scaled by its weight; the winning Order
//                        must clear the abstention gate (ensemble_min_frac/min_agree).
//   hard               — majority rule (enough agreeing votes, no tie).
//   all                — write all three ensemble columns; relabel with two_step.
// Only calls whose Order is in params.ensemble_relabel_orders (default DNA,RC) relabel
// a consensus; others leave it Unknown (see bin/ensemble_vote.py for why).
// When params.ensemble_kingdom_aware is set, only the model subset for this sample's
// organism_type votes (--include-models) — all models still appear in the table, but
// the subset alone decides. For weighted/hard the gate (ensemble_min_frac) is a
// fraction of what THAT subset could maximally cast, so it holds whether five models
// vote or one. Relabelled headers are RepeatModeler-
// style ">name#Class/Subfamily" for RepeatMasker; confident labels are preserved and
// no sequence is dropped. Also writes a voting-results table (rows = repeats, cols =
// models + ensemble winner(s)). See bin/ensemble_vote.py. The per-model vote TSVs are
// staged under their original "<id>.<model>.preds.tsv" names to label the columns.
process ENSEMBLE_VOTE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(library), path(predictions)
    path two_step_config   // optional JSON (params.ensemble_two_step_config), [] when unset
    path weights_from      // optional metrics TSV (params.ensemble_weights_from), [] when unset

    output:
    tuple val(meta), path("${meta.id}.classified_library.fa"), emit: library
    tuple val(meta), path("${meta.id}.voting_table.tsv")     , emit: table
    path 'versions.yml'                                      , emit: versions

    script:
    def weights_arg   = params.ensemble_weights      ? "--weights '${params.ensemble_weights}'" : ''
    // Config files are staged inputs: visible inside containers, resolved from the launch
    // directory, and a changed file re-runs the task on -resume.
    def wfrom_arg     = weights_from ? "--weights-from '${weights_from}' --weight-metric ${params.ensemble_weight_metric}" : ''
    def two_step_arg  = two_step_config ? "--two-step-config '${two_step_config}'" : ''
    // Kingdom-aware ensemble: restrict the ensemble decision to the model subset
    // for this sample's organism_type, which also selects the two-step settings.
    def kingdom = (meta.organism_type ?: 'animal').toLowerCase()
    def kmodels = ''
    if (params.ensemble_kingdom_aware) {
        kmodels = kingdom == 'plant'  ? (params.ensemble_models_plant  ?: '')
                : kingdom == 'fungus' ? (params.ensemble_models_fungus ?: '')
                :                       (params.ensemble_models_animal ?: '')
    }
    def include_arg = kmodels ? "--include-models '${kmodels}'" : ''
    """
    ensemble_vote.py \\
        --library ${library} \\
        --pred ${predictions} \\
        --out ${meta.id}.classified_library.fa \\
        --table ${meta.id}.voting_table.tsv \\
        --method ${params.ensemble_method} \\
        --kingdom ${kingdom} \\
        --relabel-orders '${params.ensemble_relabel_orders ?: 'all'}' \\
        --min-frac ${params.ensemble_min_frac} \\
        --min-votes ${params.ensemble_min_votes} \\
        --min-weight ${params.ensemble_min_weight} \\
        --min-agree ${params.ensemble_min_agree} \\
        ${weights_arg} ${wfrom_arg} ${two_step_arg} ${include_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${library} ${meta.id}.classified_library.fa
    printf 'repeat_id\\tEnsemble\\n' > ${meta.id}.voting_table.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
