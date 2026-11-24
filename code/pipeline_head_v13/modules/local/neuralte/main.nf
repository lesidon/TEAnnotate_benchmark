// Ensemble classifier — NeuralTE (CNN; RepeatMasker labels via --is_wicker 0).
// Classifies the curated library's Unknown fraction and emits normalized
// [base_id, RM_label] votes. Source + in-repo model come from NEURALTE_CLONE (repo).
// Conda env assets/environments/neuralte.yml (python + RepeatMasker + rmblast via
// conda; TF/Keras/etc. via pip) — or the per-tool container (--neuralte_container).
//
// NeuralTE writes `classified.info`, a CSV `#Seq_Name,Predict_Label`;
// normalize_neuralte.py maps the labels to RepeatMasker Order/Superfamily.
process NEURALTE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(library)
    path repo

    output:
    tuple val(meta), path("${meta.id}.neuralte.preds.tsv"), emit: predictions
    path "${meta.id}.neuralte.info"                       , emit: raw, optional: true
    path 'versions.yml'                                   , emit: versions

    script:
    """
    # Activate the container's baked env (read-only .sif cannot `micromamba run`);
    # no-op under the Conda fallback.
    if [ -d /opt/conda/envs/neuralte/bin ]; then
        export PATH=/opt/conda/envs/neuralte/bin:\$PATH
    fi
    # NeuralTE is python3.8; the host's ~/.local (bind-mounted with \$HOME) shadows the
    # image's pinned numpy 1.19.5 -> `np.object` crash. Ignore user-site.
    export PYTHONNOUSERSITE=1

    WORK="\$PWD"
    export HOME="\$WORK"
    export MPLCONFIGDIR="\$WORK/.mplconfig"
    export XDG_CACHE_HOME="\$WORK/.cache"
    export TF_CPP_MIN_LOG_LEVEL=3
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" neuralte_out

    # Run from a writable copy of the repo (imports from src/, writes domain files).
    cp -rL "${repo}" neuralte_src
    cd neuralte_src
    python src/Classifier.py \\
        --data "\$WORK/${library}" \\
        --model_path "\$PWD/models/NeuralTE_model.h5" \\
        --outdir "\$WORK/neuralte_out" \\
        --is_predict 1 \\
        --is_wicker 0 \\
        --thread ${task.cpus}
    cd "\$WORK"

    INFO=\$(find neuralte_out -name 'classified.info' -print -quit 2>/dev/null || true)
    if [[ -n "\${INFO}" && -s "\${INFO}" ]]; then
        cp "\${INFO}" ${meta.id}.neuralte.info
        normalize_neuralte.py --pred "\${INFO}" --out ${meta.id}.neuralte.preds.tsv
    else
        echo "WARN [${meta.id}]: NeuralTE produced no classified.info; abstaining" >&2
        : > ${meta.id}.neuralte.preds.tsv
    fi

    NTE_VER=\$( ( cd "${repo}" && git rev-parse --short HEAD ) 2>/dev/null || true )
    [[ -z "\${NTE_VER}" ]] && NTE_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neuralte: \${NTE_VER}
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.neuralte.preds.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        neuralte: "stub"
    END_VERSIONS
    """
}
