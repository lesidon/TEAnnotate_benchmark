// Stage 3 (classification) — TEclass2 (Transformer-based DL TE classifier).
// Classifies the curated consensus library's "Unknown" fraction into Wicker
// superfamilies, then emits normalized [id, RM_label] votes for the ensemble.
// GPU is selected via the config `device` key (cuda when --use_gpu, else cpu).
//
// TEclass2 (IOB-Muenster/TEclass2) has no bioconda package, no published image,
// and no redistributable license, so it cannot be vendored or containerised. The
// source and pre-trained model are fetched by TECLASS2_CLONE / TECLASS2_FETCH_MODEL
// (cached via storeDir) and passed in here as `repo` and `model`. Nextflow builds
// the Python/PyTorch env (assets/environments/teclass2.yml) on the driver. This
// module builds a run config from the repo's config.yml (architecture base) with
// the input FASTA, device, and auto-discovered checkpoint injected, copies the
// repo into the (writable) work dir, runs --classify, and normalizes to votes.
process TECLASS2 {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(library)
    path repo
    path model

    output:
    tuple val(meta), path("${meta.id}.teclass2.preds.tsv"), emit: predictions
    path "${meta.id}.teclass2.tsv"                        , emit: raw, optional: true
    path 'versions.yml'                                   , emit: versions

    script:
    def device = params.use_gpu ? 'cuda' : 'cpu'
    def ckpt   = params.teclass2_checkpoint ? "--checkpoint '${params.teclass2_checkpoint}'" : ''
    """
    WORK="\$PWD"

    # Inside the DL container the tool lives in a conda env at /opt/conda/envs/teclass2;
    # a read-only .sif cannot use `micromamba run`, so put that env's bin on PATH.
    # No-op under the Conda fallback (the dir does not exist there).
    if [ -d /opt/conda/envs/teclass2/bin ]; then
        export PATH=/opt/conda/envs/teclass2/bin:\$PATH
    fi
    # Ignore the host's ~/.local (bind-mounted with \$HOME) — a newer user-site numpy
    # etc. would shadow the image's pinned deps. Harmless under the Conda fallback.
    export PYTHONNOUSERSITE=1

    # \$HOME on HPC is often read-only from inside the env; redirect caches.
    export MPLCONFIGDIR="\$WORK/.mplconfig"
    export XDG_CACHE_HOME="\$WORK/.cache"
    export HF_HOME="\$WORK/.hf"
    export TRANSFORMERS_OFFLINE=1
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" "\$HF_HOME"

    # Build the run config from the repo's config.yml (architecture base): enable
    # classify, point the input fasta at our library, set the device, and resolve
    # the trained checkpoint inside the model dir.
    teclass2_make_config.py \\
        "${repo}/config.yml" \\
        "\$WORK/${library}" \\
        "${device}" \\
        "\$WORK/${model}" \\
        "\$WORK/run_config.yml" \\
        ${ckpt}

    # Run from a writable copy of the repo (TEclass2 expects to run from its root
    # and may write dataset caches alongside its sources).
    cp -rL "${repo}" teclass2_src
    cd teclass2_src
    # --classify needs the input FASTA via -f and the output via -o; -c supplies
    # the model/device/architecture (the input is NOT taken from the config).
    python TEclass2.py --classify \\
        -f "\$WORK/${library}" \\
        -o "\$WORK/${meta.id}.teclass2.tsv" \\
        -c "\$WORK/run_config.yml"
    cd "\$WORK"

    # Locate the prediction table (fall back to a glob if -o landed elsewhere).
    PRED=\$(ls "${meta.id}.teclass2.tsv" 2>/dev/null | head -1 || true)
    if [[ -z "\${PRED}" ]]; then
        PRED=\$(find . teclass2_src -maxdepth 3 -name '*.tsv' -newer "\$WORK/run_config.yml" 2>/dev/null | head -1 || true)
        [[ -n "\${PRED}" ]] && cp "\${PRED}" "${meta.id}.teclass2.tsv" && PRED="${meta.id}.teclass2.tsv"
    fi

    if [[ -n "\${PRED}" && -s "\${PRED}" ]]; then
        # Abstain (leave Unknown) when the top class probability is below the
        # threshold; only confident calls become votes.
        normalize_teclass2.py --pred "\${PRED}" --out ${meta.id}.teclass2.preds.tsv \\
            --min-score ${params.teclass2_min_score}
    else
        echo "WARN [${meta.id}]: TEclass2 produced no predictions; abstaining" >&2
        : > ${meta.id}.teclass2.preds.tsv
    fi

    TC_VER=\$( ( cd "${repo}" && git rev-parse --short HEAD ) 2>/dev/null || true )
    [[ -z "\${TC_VER}" ]] && TC_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teclass2: \${TC_VER}
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.teclass2.preds.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teclass2: "stub"
    END_VERSIONS
    """
}
