include { assetKey } from '../utils/asset_key'

// Clone the DeepTE source (LiLabAtVT/DeepTE) into a cached location.
// DeepTE is script-based (no bioconda/registry image), so the repo is fetched at
// run time and cached via storeDir (params.deepte_cache, in a directory named after
// the requested revision and URL: assetKey), fetched once per revision and reused
// (skipped when --deepte_repo points at an existing clone). source.tsv records the
// commit actually checked out. git is provided by a Conda env (conf/modules.config)
// so it does not depend on the host PATH. The process is errorStrategy 'ignore'
// (modules.config): a failure is non-fatal and the cascade keeps the previous
// (e.g. TEclass2) library.
process DEEPTE_CLONE {
    tag "${params.deepte_revision ?: 'master'}"
    label 'process_low'
    storeDir "${params.deepte_cache}/repo/${assetKey(params.deepte_revision ?: 'master', "${params.deepte_repo_url} ${params.deepte_revision ?: 'master'}")}"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'DeepTE'      , emit: repo
    path 'source.tsv'  , emit: manifest
    path 'versions.yml', emit: versions

    script:
    def rev = params.deepte_revision ?: 'master'
    """
    git clone "${params.deepte_repo_url}" DeepTE
    ( cd DeepTE && git -c advice.detachedHead=false checkout "${rev}" )

    SHA=\$(cd DeepTE && git rev-parse HEAD)
    printf 'repo\\t%s\\nrequested\\t%s\\ncommit\\t%s\\n' "${params.deepte_repo_url}" "${rev}" "\${SHA}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte_repo: \${SHA}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p DeepTE
    printf 'repo\\t%s\\nrequested\\t%s\\ncommit\\tstub\\n' "${params.deepte_repo_url}" "${params.deepte_revision}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte_repo: stub
    END_VERSIONS
    """
}
