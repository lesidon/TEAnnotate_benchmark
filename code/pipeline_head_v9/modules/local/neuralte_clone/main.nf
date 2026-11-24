// Clone the NeuralTE source (CSU-KangHu/NeuralTE) into a cached location. NeuralTE
// is script-based (no bioconda/registry image) and ships its pre-trained model in
// the repo (models/NeuralTE_model.h5), so the clone provides BOTH source and model.
// Fetched once and cached via storeDir (params.neuralte_cache); skipped when
// --neuralte_repo points at an existing clone. git is provided by a Conda env
// (conf/modules.config). errorStrategy 'ignore': a failure is non-fatal and the
// ensemble simply loses NeuralTE's vote.
process NEURALTE_CLONE {
    tag "${params.neuralte_revision ?: 'master'}"
    label 'process_low'
    storeDir "${params.neuralte_cache}/repo"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'NeuralTE'    , emit: repo
    path 'versions.yml', emit: versions

    script:
    def rev = params.neuralte_revision ?: 'master'
    """
    git clone "${params.neuralte_repo_url}" NeuralTE
    ( cd NeuralTE && git checkout "${rev}" )

    SHA=\$( ( cd NeuralTE && git rev-parse --short HEAD ) 2>/dev/null || true )
    [ -z "\${SHA}" ] && SHA=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neuralte_repo: \${SHA}
    END_VERSIONS
    """
}
