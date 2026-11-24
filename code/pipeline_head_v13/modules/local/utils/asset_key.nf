// Cache directory names for fetched assets (storeDir). storeDir reuses whatever it
// finds, so the name carries the requested identity: a readable label (release,
// commit or branch) plus a hash of the full identity (URL + revision). Changing a
// pinned revision, URL or Drive ID then fetches into a new directory instead of
// reusing the old asset. -stub runs get their own 'stub-' directories, so placeholder
// assets never reach a real run.
def assetKey(String label, String identity) {
    def l = (label ==~ /[0-9a-f]{40}/) ? label.take(12) : label.replaceAll(/[^A-Za-z0-9._-]/, '_')
    (workflow.stubRun ? 'stub-' : '') + l + '-' + identity.md5().take(8)
}

// Rows "supplied <param> <path>" of the asset manifest (pipeline_info/asset_manifest.tsv)
// for every asset given as a local path instead of being fetched: no source.tsv
// describes those.
def suppliedAssets() {
    ['pfam_db', 'hite_ltrdetector_dir', 'teclass2_repo', 'teclass2_model', 'deepte_repo',
     'deepte_model', 'create_repo', 'create_model', 'neuralte_repo']
        .findAll { params.containsKey(it) && params[it] }
        .collect { "supplied\t${it}\t${file(params[it])}" }
}
