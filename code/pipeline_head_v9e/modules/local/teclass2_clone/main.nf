// Clone the TEclass2 source (IOB-Muenster/TEclass2) into a cached location.
// TEclass2 is script-based (no bioconda/registry image, no redistributable
// license), so the repo is fetched at run time instead of vendored. `storeDir`
// persists the clone under params.teclass2_cache so it is fetched once and reused
// across runs (and skipped entirely when --teclass2_repo points at an existing
// clone). git is provided by a Conda env (conf/modules.config) so it does not
// depend on the host PATH (HPC nodes often hide git behind `module load`). The
// process is errorStrategy 'ignore' (set in modules.config): a clone failure is
// non-fatal and CURATION falls back to the curated library.
process TECLASS2_CLONE {
    tag "${params.teclass2_revision ?: 'main'}"
    label 'process_low'
    storeDir "${params.teclass2_cache}/repo"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'TEclass2'    , emit: repo
    path 'versions.yml', emit: versions

    script:
    def rev = params.teclass2_revision ?: 'main'
    """
    git clone "${params.teclass2_repo_url}" TEclass2
    ( cd TEclass2 && git checkout "${rev}" )

    SHA=\$( ( cd TEclass2 && git rev-parse --short HEAD ) 2>/dev/null || true )
    [ -z "\${SHA}" ] && SHA=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teclass2_repo: \${SHA}
    END_VERSIONS
    """
}
