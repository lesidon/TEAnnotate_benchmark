// Cascade classifier — DeepTE (CNN, organism-aware). Classifies the consensus
// library and emits normalized [id, RM_label] votes for the ensemble. The
// organism_type selects DeepTE's species model: plant -> P /
// Plants_model, fungus -> F / Fungi_model, animal -> M / Metazoans_model, else
// O / Others_model. Conda env (assets/environments/deepte.yml); the source is
// fetched by DEEPTE_CLONE (repo). The `model` dir is the single organism-matched
// group fetched by DEEPTE_FETCH_MODEL and attached per-sample by the subworkflow,
// so it holds exactly the <Group>_model subdir this sample needs.
process DEEPTE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(library), path(model)
    path repo

    output:
    tuple val(meta), path("${meta.id}.deepte.preds.tsv"), emit: predictions
    path "${meta.id}.deepte.txt"                        , emit: raw, optional: true
    path 'versions.yml'                                 , emit: versions

    script:
    def sp  = meta.organism_type == 'plant'  ? 'P'
            : meta.organism_type == 'fungus' ? 'F'
            : meta.organism_type == 'animal' ? 'M' : 'O'
    def mdl = meta.organism_type == 'plant'  ? 'Plants_model'
            : meta.organism_type == 'fungus' ? 'Fungi_model'
            : meta.organism_type == 'animal' ? 'Metazoans_model' : 'Others_model'
    """
    # Inside the DL container the tool lives in a conda env at /opt/conda/envs/deepte;
    # a read-only .sif cannot use `micromamba run`, so put that env's bin on PATH.
    # No-op under the Conda fallback (the dir does not exist there).
    if [ -d /opt/conda/envs/deepte/bin ]; then
        export PATH=/opt/conda/envs/deepte/bin:\$PATH
    fi
    # Ignore host ~/.local (mounted with \$HOME) so it can't shadow pinned deps.
    export PYTHONNOUSERSITE=1

    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" deepte_wd deepte_out

    python "${repo}/DeepTE.py" \\
        -d "\$PWD/deepte_wd" \\
        -o "\$PWD/deepte_out" \\
        -i ${library} \\
        -sp ${sp} \\
        -m_dir "${model}/${mdl}"

    PRED=\$(find deepte_out -name 'opt_DeepTE.txt' -print -quit 2>/dev/null || true)
    if [[ -n "\${PRED}" && -s "\${PRED}" ]]; then
        cp "\${PRED}" ${meta.id}.deepte.txt
        normalize_deepte.py --pred "\${PRED}" --out ${meta.id}.deepte.preds.tsv
    else
        echo "WARN [${meta.id}]: DeepTE produced no opt_DeepTE.txt; abstaining" >&2
        : > ${meta.id}.deepte.preds.tsv
    fi

    DT_VER=\$( ( cd "${repo}" && git rev-parse --short HEAD ) 2>/dev/null || true )
    [[ -z "\${DT_VER}" ]] && DT_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte: \${DT_VER}
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.deepte.preds.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepte: "stub"
    END_VERSIONS
    """
}
