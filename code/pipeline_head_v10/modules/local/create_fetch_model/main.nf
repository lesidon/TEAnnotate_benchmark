include { assetKey } from '../utils/asset_key'

// Download + extract CREATE's pre-trained models (the 9 hierarchical classifiers)
// into a cached directory. The models are hosted on Google Drive, so the download
// uses gdown (provided by a Conda env; conf/modules.config), which handles Drive's
// large-file confirmation token that plain curl cannot. Cached via storeDir
// (params.create_cache, in a directory named after the Drive ID: assetKey); skipped
// when --create_model points at an existing dir. source.tsv records the Drive ID and
// the archive's SHA-256. errorStrategy 'ignore' (modules.config): non-fatal — CREATE
// abstains on failure.
//
// The tarball extracts to a `model/` subdir in the demo; we flatten it so the
// emitted `create_model` dir directly contains the classifiers (CREATE's -d arg).
process CREATE_FETCH_MODEL {
    label 'process_low'
    storeDir "${params.create_cache}/model/${assetKey('model', params.create_model_gdrive_id)}"

    // Ordering gate (value unused): a trigger from CURATION so this run-time setup
    // starts only AFTER TEtrimmer. Without a data dependency Nextflow would launch
    // it at t=0, before de novo. run_ensemble.nf passes an immediate token.
    input:
    val ready

    output:
    path 'create_model', emit: model
    path 'source.tsv'  , emit: manifest
    path 'versions.yml', emit: versions

    script:
    """
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p create_model "\$XDG_CACHE_HOME"
    cd create_model

    gdown "https://drive.google.com/uc?id=${params.create_model_gdrive_id}" -O model.tar.gz

    if tar -tzf model.tar.gz >/dev/null 2>&1; then
        SUM=\$(sha256sum model.tar.gz | cut -d' ' -f1)
        tar -xzf model.tar.gz
        rm -f model.tar.gz
    else
        echo "ERROR: CREATE model did not download as a valid tar.gz (Drive quota?)" >&2
        exit 1
    fi

    # Flatten a single top-level 'model/' dir so the classifiers sit directly here.
    if [ -d model ]; then
        mv model/* . 2>/dev/null || true
        rmdir model 2>/dev/null || true
    fi
    cd ..

    printf 'google_drive_id\\t%s\\narchive_sha256\\t%s\\n' "${params.create_model_gdrive_id}" "\${SUM}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_model: "google-drive:${params.create_model_gdrive_id}"
        create_model_sha256: \${SUM}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p create_model
    printf 'google_drive_id\\t%s\\narchive_sha256\\tstub\\n' "${params.create_model_gdrive_id}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_model: "google-drive:${params.create_model_gdrive_id}"
    END_VERSIONS
    """
}
