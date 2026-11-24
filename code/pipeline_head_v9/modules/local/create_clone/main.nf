// Clone the CREATE source (yangqi-cs/CREATE) into a cached location. CREATE's
// pre-trained models are NOT in the repo (downloaded separately by
// CREATE_FETCH_MODEL). git is provided by a Conda env (conf/modules.config);
// errorStrategy 'ignore' (modules.config) makes a failure non-fatal — the ensemble
// degrades gracefully (CREATE simply abstains). storeDir caches the clone
// (params.create_cache).
process CREATE_CLONE {
    tag "${params.create_revision ?: 'master'}"
    label 'process_low'
    storeDir "${params.create_cache}/repo"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'CREATE'      , emit: repo
    path 'versions.yml', emit: versions

    script:
    def rev = params.create_revision ?: 'master'
    """
    git clone "${params.create_repo_url}" CREATE
    ( cd CREATE && git checkout "${rev}" )

    SHA=\$( ( cd CREATE && git rev-parse --short HEAD ) 2>/dev/null || true )
    [ -z "\${SHA}" ] && SHA=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_repo: \${SHA}
    END_VERSIONS
    """
}
