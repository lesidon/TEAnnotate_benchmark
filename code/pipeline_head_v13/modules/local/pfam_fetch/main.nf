include { assetKey } from '../utils/asset_key'

// Cache label: the Pfam release named in the URL (releases/Pfam38.2/...), else 'pfam'.
def pfamLabel() {
    def m = params.pfam_hmm_url =~ /Pfam[0-9]+(\.[0-9]+)?/
    m.find() ? m.group() : 'pfam'
}

// Fetch + index the Pfam-A database for TEtrimmer. TEtrimmer needs PFAM to
// determine TE orientation, and analyze.create_dir() ABORTS (returns early, no
// consensus) when it is missing — "The database is therefore mandatory." Its own
// downloader (requests/urllib to EBI) fails on offline/locked-down HPC nodes, so
// we fetch the two files with curl (which works where TEtrimmer's downloader does
// not) and run hmmpress here, then hand TEtrimmer a ready, pre-indexed --pfam_dir
// (so it neither downloads nor writes into the staged, read-only copy).
//
// Cached via storeDir (params.pfam_cache), in a directory named after the release
// and URLs (assetKey); source.tsv records the URLs and the files' SHA-256. Override
// with --pfam_db to point at an already-downloaded local Pfam directory. NOT
// errorStrategy 'ignore': PFAM is mandatory, so a failure here should stop the run
// (use --pfam_db instead).
process PFAM_FETCH {
    label 'process_low'
    storeDir "${params.pfam_cache}/${assetKey(pfamLabel(), "${params.pfam_hmm_url} ${params.pfam_dat_url}")}"

    output:
    path 'pfam_db'     , emit: db
    path 'source.tsv'  , emit: manifest
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p pfam_db
    cd pfam_db

    curl -fL --retry 3 -o Pfam-A.hmm.gz     "${params.pfam_hmm_url}"
    curl -fL --retry 3 -o Pfam-A.hmm.dat.gz "${params.pfam_dat_url}"
    SUM_HMM=\$(sha256sum Pfam-A.hmm.gz | cut -d' ' -f1)
    SUM_DAT=\$(sha256sum Pfam-A.hmm.dat.gz | cut -d' ' -f1)
    gzip -df Pfam-A.hmm.gz
    gzip -df Pfam-A.hmm.dat.gz

    # Binary index files (.h3f/.h3i/.h3m/.h3p) for hmmscan / pfam_scan.pl. Having
    # them present means TEtrimmer skips hmmpress and never writes into this dir.
    hmmpress -f Pfam-A.hmm
    NFAM=\$(grep -c '^NAME' Pfam-A.hmm || true)
    cd ..

    printf 'hmm_url\\t%s\\nhmm_gz_sha256\\t%s\\ndat_url\\t%s\\ndat_gz_sha256\\t%s\\nfamilies\\t%s\\n' \\
        "${params.pfam_hmm_url}" "\${SUM_HMM}" "${params.pfam_dat_url}" "\${SUM_DAT}" "\${NFAM}" > source.tsv

    HVER=\$( (hmmpress -h 2>&1 || true) | grep -oE 'HMMER [0-9.]+' | head -1 | sed 's/HMMER //' || true)
    [ -z "\${HVER}" ] && HVER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \${HVER}
        pfam: "${params.pfam_hmm_url}"
        pfam_families: \${NFAM}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pfam_db
    touch pfam_db/Pfam-A.hmm pfam_db/Pfam-A.hmm.dat
    printf 'hmm_url\\t%s\\nhmm_gz_sha256\\tstub\\n' "${params.pfam_hmm_url}" > source.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pfam: "${params.pfam_hmm_url}"
    END_VERSIONS
    """
}
