// Ensemble classifier — CREATE (attention-based CNN+RNN; Wicker taxonomy).
// Classifies the curated consensus library and emits normalized [id, RM_label]
// votes for the hard-voting ensemble. The source is cloned by CREATE_CLONE and the
// pre-trained models fetched by CREATE_FETCH_MODEL (passed in as repo / model).
// Conda env: assets/environments/create.yml (TensorFlow via pip).
//
// CREATE writes `relabeled_seqs.fasta` whose headers are
// ">{predicted_label}\t{original_description}"; normalize_create.py turns those into
// the per-sequence vote TSV.
process CREATE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(library)
    path repo
    path model

    output:
    tuple val(meta), path("${meta.id}.create.preds.tsv"), emit: predictions
    path "${meta.id}.create.relabeled.fasta"            , emit: raw, optional: true
    path 'versions.yml'                                 , emit: versions

    script:
    // NB: use an explicit null test, not `?:` — 0.0 is "falsy" in Groovy, so
    // `params.create_prob_thr ?: 0.70` would wrongly override a 0.0 (no-threshold) setting.
    def prob = params.create_prob_thr != null ? params.create_prob_thr : 0.70
    """
    WORK="\$PWD"

    # Inside the DL container the tool lives in a conda env at /opt/conda/envs/create;
    # a read-only .sif cannot use `micromamba run`, so put that env's bin on PATH.
    # No-op under the Conda fallback (the dir does not exist there).
    if [ -d /opt/conda/envs/create/bin ]; then
        export PATH=/opt/conda/envs/create/bin:\$PATH
    fi
    # Ignore host ~/.local (mounted with \$HOME) so it can't shadow pinned deps.
    export PYTHONNOUSERSITE=1

    export HOME="\$WORK"
    export MPLCONFIGDIR="\$WORK/.mplconfig"
    export XDG_CACHE_HOME="\$WORK/.cache"
    export TF_CPP_MIN_LOG_LEVEL=3
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" create_out

    # CREATE's -d expects the dir holding the classifiers; tolerate a nested model/.
    MDIR="\$WORK/${model}"
    [ -d "\$MDIR/model" ] && MDIR="\$MDIR/model"

    # CREATE reads each sequence's TRUE label from the text BEFORE '#' (a benchmark
    # convention) and crashes (find_node_path IndexError) when it isn't a taxonomy
    # node — which is exactly our '#Unknown'. Rewrite headers to ">TE#<orig_id>"
    # (valid root node "TE" before '#', real id after) so it runs; the id survives
    # into relabeled_seqs.fasta and normalize_create.py recovers it after '#'.
    create_prep_headers.py "\$WORK/${library}" "\$WORK/create_in.fa"

    # Run from a writable copy of the repo (scripts import sibling modules and may
    # write caches alongside their sources).
    cp -rL "${repo}" create_src
    cd create_src
    python create_test.py \\
        -i "\$WORK/create_in.fa" \\
        -d "\$MDIR" \\
        -o "\$WORK/create_out" \\
        -p ${prob}
    cd "\$WORK"

    PRED=\$(find create_out -name 'relabeled_seqs.fasta' -print -quit 2>/dev/null || true)
    if [[ -n "\${PRED}" && -s "\${PRED}" ]]; then
        cp "\${PRED}" ${meta.id}.create.relabeled.fasta
        normalize_create.py --pred "\${PRED}" --out ${meta.id}.create.preds.tsv
    else
        echo "WARN [${meta.id}]: CREATE produced no relabeled_seqs.fasta; abstaining" >&2
        : > ${meta.id}.create.preds.tsv
    fi

    CR_VER=\$( ( cd "${repo}" && git rev-parse --short HEAD ) 2>/dev/null || true )
    [[ -z "\${CR_VER}" ]] && CR_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create: \${CR_VER}
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.create.preds.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create: "stub"
    END_VERSIONS
    """
}
