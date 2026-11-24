// Clone the DeepTE source (LiLabAtVT/DeepTE) into a cached location.
// DeepTE is script-based (no bioconda/registry image), so the repo is fetched at
// run time and cached via storeDir (params.deepte_cache), fetched once and reused
// (skipped when --deepte_repo points at an existing clone). git is provided by a
// Conda env (conf/modules.config) so it does not depend on the host PATH. The
// process is errorStrategy 'ignore' (modules.config): a failure is non-fatal and
// the cascade keeps the previous (e.g. TEclass2) library.
process DEEPTE_CLONE {
    tag "${params.deepte_revision ?: 'master'}"
    label 'process_low'
    storeDir "${params.deepte_cache}/repo"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'DeepTE'      , emit: repo
    path 'versions.yml', emit: versions

    script:
    def rev = params.deepte_revision ?: 'master'
    """
    git clone "${params.deepte_repo_url}" DeepTE
    ( cd DeepTE && git checkout "${rev}" )

    SHA=\$( ( cd DeepTE && git rev-parse --short HEAD ) 2>/dev/null || true )
    [ -z "\${SHA}" ] && SHA=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte_repo: \${SHA}
    END_VERSIONS
    """
}
