include { assetKey } from '../utils/asset_key'

// Build the LtrDetector directory that kanghu/hite:3.3.3 lacks. The image ships
// bin/FiLTR-main/bin/LtrDetector/ without its bin/ subdirectory, so HiTE's default LTR
// step (FiLTR) finds nothing and HiTE still exits 0 with no LTR families. This fetches
// only that directory from the HiTE commit the image code matches (sparse, shallow:
// ~1.5 MB), makes the prebuilt binary executable and wraps it so it finds libgomp in the
// image's conda env. HITE then runs from a copy of /HiTE that uses this directory.
// Skipped when --hite_ltrdetector_dir points at a directory built the same way.
// Cached via storeDir in a directory named after the revision and URL (assetKey), so
// a different --hite_revision gets its own directory; source.tsv records the commit
// and the binary's SHA-256. git comes from a Conda env (conf/modules.config).
process HITE_LTRDETECTOR {
    tag "${params.hite_revision}"
    label 'process_single'
    storeDir "${params.hite_cache}/ltrdetector/${assetKey(params.hite_revision, "${params.hite_repo_url} ${params.hite_revision}")}"

    output:
    path 'LtrDetector'  , emit: dir
    path 'source.tsv'   , emit: manifest
    path 'versions.yml' , emit: versions

    script:
    def sub = 'bin/FiLTR-main/bin/LtrDetector'
    """
    git init -q hite_src
    cd hite_src
    git remote add origin "${params.hite_repo_url}"
    git sparse-checkout set "${sub}"
    git fetch -q --depth 1 --filter=blob:none origin "${params.hite_revision}"
    git -c advice.detachedHead=false checkout -q FETCH_HEAD
    SHA=\$(git rev-parse HEAD)
    cd ..
    mv "hite_src/${sub}" LtrDetector
    rm -rf hite_src

    # git stores the binary non-executable, and it needs libgomp from the image's env.
    D=LtrDetector/bin
    [[ -f \$D/LtrDetector ]] || { echo "ERROR: ${sub}/bin/LtrDetector not in ${params.hite_repo_url} at ${params.hite_revision}" >&2; exit 1; }
    mv \$D/LtrDetector \$D/LtrDetector.bin
    chmod +x \$D/LtrDetector.bin
    cat > \$D/LtrDetector <<-'EOF'
    #!/bin/sh
    # Added by nf-core/teannotate: the prebuilt LtrDetector needs libgomp.so.1, which is
    # in the HiTE image's conda env but not on the loader path.
    LD_LIBRARY_PATH=/opt/conda/envs/HiTE/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH} exec "\$(dirname "\$0")/LtrDetector.bin" "\$@"
    EOF
    chmod +x \$D/LtrDetector

    printf 'repo\\t%s\\nrequested\\t%s\\ncommit\\t%s\\nbinary_sha256\\t%s\\n' \\
        "${params.hite_repo_url}" "${params.hite_revision}" "\${SHA}" "\$(sha256sum \$D/LtrDetector.bin | cut -d' ' -f1)" > source.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hite_ltrdetector: \${SHA}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p LtrDetector/bin
    printf '#!/bin/sh\\necho " -fasta is a required argument"; exit 1\\n' > LtrDetector/bin/LtrDetector
    chmod +x LtrDetector/bin/LtrDetector
    printf 'repo\\t%s\\nrequested\\t%s\\ncommit\\tstub\\n' "${params.hite_repo_url}" "${params.hite_revision}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hite_ltrdetector: stub
    END_VERSIONS
    """
}
