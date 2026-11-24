include { assetKey } from '../utils/asset_key'

// Model archive URL for a DeepTE group (Plants / Metazoans / Fungi / Others).
def deepteModelUrl(group) {
    group == 'Plants'    ? params.deepte_model_url_plants
  : group == 'Metazoans' ? params.deepte_model_url_metazoans
  : group == 'Fungi'     ? params.deepte_model_url_fungi
  :                        params.deepte_model_url_others
}

// Download + extract ONE DeepTE pre-trained model group (Plants / Metazoans /
// Fungi / Others), selected by the sample's organism_type. DeepTE's models are
// already segmented by taxonomic lineage, so we fetch only the group actually
// needed instead of the full multi-GB suite. Cached per group and URL via storeDir
// (params.deepte_cache/model/<group>/<assetKey>), so a given group downloads once and
// is reused across samples and runs; skipped entirely when --deepte_model points at
// a pre-extracted dir. source.tsv records the URL and the archive's SHA-256. The
// download lands in the host cache, which Nextflow then bind-mounts into the DeepTE
// task container (the model is NOT baked into the image). curl is provided by a
// Conda env (conf/modules.config).
// errorStrategy 'ignore' (modules.config): non-fatal — the cascade degrades.
process DEEPTE_FETCH_MODEL {
    tag "${group}"
    label 'process_low'
    storeDir "${params.deepte_cache}/model/${group}/${assetKey('model', deepteModelUrl(group))}"

    input:
    val group   // Plants | Metazoans | Fungi | Others

    output:
    tuple val(group), path('deepte_models'), emit: model
    path 'source.tsv'                      , emit: manifest
    path 'versions.yml'                    , emit: versions

    script:
    def url = deepteModelUrl(group)
    """
    mkdir -p deepte_models
    cd deepte_models
    fname=\$(basename "${url}")
    curl -fL --retry 3 -o "\${fname}" "${url}"
    if tar -tzf "\${fname}" >/dev/null 2>&1; then
        SUM=\$(sha256sum "\${fname}" | cut -d' ' -f1)
        tar -xzf "\${fname}"
        rm -f "\${fname}"
    else
        echo "ERROR: ${url} did not download a valid tar.gz" >&2
        exit 1
    fi
    cd ..

    printf 'url\\t%s\\narchive_sha256\\t%s\\n' "${url}" "\${SUM}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte_model_${group}: "cyverse"
        deepte_model_${group}_sha256: \${SUM}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p deepte_models
    printf 'url\\t%s\\narchive_sha256\\tstub\\n' "${deepteModelUrl(group)}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte_model_${group}: "cyverse"
    END_VERSIONS
    """
}
