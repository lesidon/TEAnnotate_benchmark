process EARLGREY {
    tag "$meta.id"
    label 'process_high'

    // No biocontainers image is available for earl-grey, so this process always
    // runs in a Conda environment built from assets/environments/earl-grey.yml.
    // The singularity profile keeps conda.enabled = true precisely so this
    // (container-less) process still resolves; all other processes prefer their
    // container and fall back to Conda only when the image is unavailable.
    conda "${projectDir}/assets/environments/earl-grey.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.earlgrey_library.fa"), emit: library
    path "earlgrey_out/"                                    , emit: output_dir
    path 'versions.yml'                                     , emit: versions

    script:
    // earlGreyLibConstruct runs ONLY de novo TE detection + BEAT consensus
    // curation (it calls RepeatModeler2 internally), skipping the genome-wide
    // RepeatMasker annotation, masking, defragmentation and plots — we want just
    // the curated consensus library.
    """
    # Some bioconda earlgrey builds do not symlink earlGreyLibConstruct onto PATH
    # even though the script ships inside the package — resolve it explicitly.
    EG_LIB=\$(command -v earlGreyLibConstruct 2>/dev/null || true)
    if [[ -z "\${EG_LIB}" ]]; then
        EG_LIB=\$(find "\${CONDA_PREFIX:-/opt}" -name earlGreyLibConstruct -type f 2>/dev/null | head -1)
    fi
    if [[ -z "\${EG_LIB}" ]]; then
        echo "ERROR [${meta.id}]: earlGreyLibConstruct not found in the EarlGrey install" >&2
        exit 1
    fi

    # RepeatMasker must have a real Dfam (famdb) library: RepeatModeler's RepeatClassifier
    # uses it to label the de novo families. bioconda repeatmasker >=4.2.4 ships NO Dfam and
    # has no post-link to fetch one (4.2.3 and earlier did), so an env solved against it has
    # an empty Libraries/. Fail loudly here rather than let EarlGrey die inside dfamCheck with
    # a bare "find: ... No such file or directory", or — worse — run with no library and
    # return everything as Unknown.
    RM_BIN=\$(command -v RepeatMasker || true)
    FAMDB=""
    [[ -n "\${RM_BIN}" ]] && FAMDB=\$(echo "\${RM_BIN}" | sed 's|/bin/RepeatMasker|/share/RepeatMasker/Libraries/famdb|')
    if [[ -z "\${FAMDB}" || ! -d "\${FAMDB}" ]]; then
        echo "ERROR [${meta.id}]: RepeatMasker has no Dfam/famdb library at '\${FAMDB:-<unresolved>}'." >&2
        echo "  bioconda repeatmasker >=4.2.4 ships no Dfam database and nothing downloads it." >&2
        echo "  Install the Dfam famdb partition(s) into that directory (and pin repeatmasker" >&2
        echo "  in assets/environments/earl-grey.yml) before running EarlGrey." >&2
        exit 1
    fi
    # dfamCheck() exits 2 when famdb holds EXACTLY 2 files (dfam*_full.0.h5 + rmlib.config),
    # printing a manual Dfam-3.9 setup guide. The sentinel makes the count != 2 so a properly
    # populated famdb is accepted and de novo proceeds.
    touch "\${FAMDB}/.earlgrey.config.complete" 2>/dev/null || true

    # -r <clade> (Dfam known-repeat pre-mask) is omitted: the minimal partition 0
    # has no clade families. De novo construction (RepeatModeler2 + BEAT) runs on
    # its own; supply a full Dfam to enable clade-aware masking.
    bash "\${EG_LIB}" \\
        -g "${fasta}" \\
        -s "${meta.id}" \\
        -o earlgrey_out \\
        -t ${task.cpus}

    # The curated consensus library lands in the summaryFiles directory. The exact
    # name varies slightly across versions (*-families.fa / *-families.fa.strained).
    LIBRARY=\$(find earlgrey_out -path "*summaryFiles*" \\( -name "*-families.fa" -o -name "*-families.fa.strained" \\) -print -quit 2>/dev/null || true)
    if [[ -z "\${LIBRARY}" ]]; then
        LIBRARY=\$(find earlgrey_out \\( -name "*-families.fa.strained" -o -name "*-families.fa" \\) -print -quit 2>/dev/null || true)
    fi

    if [[ -n "\${LIBRARY}" ]]; then
        cp "\${LIBRARY}" "${meta.id}.earlgrey_library.fa"
    else
        echo "WARN [${meta.id}]: earlGreyLibConstruct produced no TE library" >&2
        touch "${meta.id}.earlgrey_library.fa"
    fi

    # earlGrey --version is unreliable; read the version from conda metadata.
    EG_VER=\$(ls "\${CONDA_PREFIX:-}"/conda-meta/earlgrey-*.json 2>/dev/null | sed -E 's@.*/earlgrey-([0-9][^-]*)-.*@\\1@' | head -1 || true)
    [[ -z "\${EG_VER}" ]] && EG_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        earlgrey: \${EG_VER}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p earlgrey_out
    cat > "${meta.id}.earlgrey_library.fa" <<-EOF
    >${meta.id}_stub_EG#LTR/Copia
    TTGAATTGAATTGAATTGAATTGAATTGAATTGAATTGAATTGAATTGAATTGAATTGAA
    AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCG
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        earlgrey: "stub"
    END_VERSIONS
    """
}
