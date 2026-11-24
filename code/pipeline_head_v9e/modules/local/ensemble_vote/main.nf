// Ensemble for TE re-classification (Stage 3): hard OR weighted voting.
// Each enabled classifier (TEclass2, DeepTE, CREATE, NeuralTE, Terrier) contributes
// a normalized [base_id, RM_label] vote TSV. For every consensus still labelled
// Unknown, the class is chosen by params.ensemble_method:
//   weighted (default) — each model's vote scaled by its weight (NeuralTE/Terrier
//                        lead); the winning Order must clear the abstention gate.
//   hard               — majority rule (enough agreeing votes, no tie).
//   all                — write BOTH ensemble columns; relabel the FASTA with weighted.
// Votes are cast at the Order level and refined to superfamily when the winners
// agree. When params.ensemble_kingdom_aware is set, only the best-performing model
// subset for this sample's organism_type votes (--include-models, from
// ensemble_diagnostics.py) — all models still appear in the table, but the subset
// alone decides. The gate (ensemble_min_frac) is a fraction of what THAT subset
// could maximally cast, so it holds whether five models vote or one (fungus =
// Terrier alone); an absolute floor could not, which is why it was previously
// disabled here — leaving no gate at all. Relabelled headers are RepeatModeler-
// style ">name#Class/Subfamily" for RepeatMasker; confident labels are preserved and
// no sequence is dropped. Also writes a voting-results table (rows = repeats, cols =
// models + ensemble winner(s)). See bin/ensemble_vote.py. The per-model vote TSVs are
// staged under their original "<id>.<model>.preds.tsv" names to label the columns.
process ENSEMBLE_VOTE {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(library), path(predictions)

    output:
    tuple val(meta), path("${meta.id}.classified_library.fa"), emit: library
    tuple val(meta), path("${meta.id}.voting_table.tsv")     , emit: table
    path 'versions.yml'                                      , emit: versions

    script:
    def weights_arg   = params.ensemble_weights      ? "--weights '${params.ensemble_weights}'" : ''
    def wfrom_arg     = params.ensemble_weights_from ? "--weights-from ${params.ensemble_weights_from} --weight-metric ${params.ensemble_weight_metric}" : ''
    // Kingdom-aware ensemble: restrict the ensemble decision to the best-performing
    // model subset for this sample's organism_type (from ensemble_diagnostics.py).
    def kmodels = ''
    if (params.ensemble_kingdom_aware) {
        def kingdom = (meta.organism_type ?: 'animal').toLowerCase()
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
        --min-frac ${params.ensemble_min_frac} \\
        --min-votes ${params.ensemble_min_votes} \\
        --min-weight ${params.ensemble_min_weight} \\
        --min-agree ${params.ensemble_min_agree} \\
        ${weights_arg} ${wfrom_arg} ${include_arg}

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
