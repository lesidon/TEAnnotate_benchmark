// Download + extract the pre-trained TEclass2 model archive.
// The model (~hundreds of MB) is published as a .tgz by the Muenster group; it is
// fetched at run time and cached via `storeDir` (params.teclass2_cache) so it is
// downloaded once and reused (skipped when --teclass2_model points at an existing
// extracted dir). curl is provided by a Conda env (conf/modules.config) so it does
// not depend on the host PATH. The process is errorStrategy 'ignore' (set in
// modules.config): a download failure is non-fatal and CURATION falls back to the
// curated library.
process TECLASS2_FETCH_MODEL {
    label 'process_low'
    storeDir "${params.teclass2_cache}/model"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'teclass2_model', emit: model
    path 'versions.yml'  , emit: versions

    script:
    """
    mkdir -p teclass2_model
    curl -fL --retry 3 -o model.tgz "${params.teclass2_model_url}"

    # Guard against an HTML error page being saved as the archive.
    if ! tar -tzf model.tgz >/dev/null 2>&1; then
        echo "ERROR: downloaded file is not a valid gzip tar archive (check teclass2_model_url)" >&2
        exit 1
    fi
    tar -xzf model.tgz -C teclass2_model
    rm -f model.tgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teclass2_model_url: "${params.teclass2_model_url}"
    END_VERSIONS
    """
}
