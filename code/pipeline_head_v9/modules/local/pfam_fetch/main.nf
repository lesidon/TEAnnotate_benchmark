// Fetch + index the Pfam-A database for TEtrimmer. TEtrimmer needs PFAM to
// determine TE orientation, and analyze.create_dir() ABORTS (returns early, no
// consensus) when it is missing — "The database is therefore mandatory." Its own
// downloader (requests/urllib to EBI) fails on offline/locked-down HPC nodes, so
// we fetch the two files with curl (which works where TEtrimmer's downloader does
// not) and run hmmpress here, then hand TEtrimmer a ready, pre-indexed --pfam_dir
// (so it neither downloads nor writes into the staged, read-only copy).
//
// Cached via storeDir (params.pfam_cache). Override with --pfam_db to point at an
// already-downloaded local Pfam directory. NOT errorStrategy 'ignore': PFAM is
// mandatory, so a failure here should stop the run (use --pfam_db instead).
process PFAM_FETCH {
    label 'process_low'
    storeDir "${params.pfam_cache}"

    output:
    path 'pfam_db'     , emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p pfam_db
    cd pfam_db

    curl -fL --retry 3 -o Pfam-A.hmm.gz     "${params.pfam_hmm_url}"
    curl -fL --retry 3 -o Pfam-A.hmm.dat.gz "${params.pfam_dat_url}"
    gzip -df Pfam-A.hmm.gz
    gzip -df Pfam-A.hmm.dat.gz

    # Binary index files (.h3f/.h3i/.h3m/.h3p) for hmmscan / pfam_scan.pl. Having
    # them present means TEtrimmer skips hmmpress and never writes into this dir.
    hmmpress -f Pfam-A.hmm
    cd ..

    HVER=\$( (hmmpress -h 2>&1 || true) | grep -oE 'HMMER [0-9.]+' | head -1 | sed 's/HMMER //' || true)
    [ -z "\${HVER}" ] && HVER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \${HVER}
        pfam: "${params.pfam_hmm_url}"
    END_VERSIONS
    """
}
