// Ensemble classifier — Terrier (PyTorch/corgi; RepeatMasker-schema labels).
// Classifies the curated library's Unknown fraction and emits normalized
// [base_id, RM_label] votes. Terrier is the `bio-terrier` pip package (no repo, no
// separate model input): its weights auto-download on first use and are baked into
// the per-tool image (/opt/cache/torchapps). Conda env assets/environments/
// terrier.yml, or the container (--terrier_container).
//
// Terrier writes --output-csv (one row per sequence, per-class probability columns);
// normalize_terrier.py takes the arg-max class and maps it to RepeatMasker labels.
process TERRIER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(library)

    output:
    tuple val(meta), path("${meta.id}.terrier.preds.tsv"), emit: predictions
    path "${meta.id}.terrier.csv"                        , emit: raw, optional: true
    path 'versions.yml'                                  , emit: versions

    script:
    def minp = params.terrier_min_prob ?: 0.0
    // Terrier's image bundles a Torch built against a newer CUDA than many cluster
    // drivers provide, so BOTH the checkpoint load (torch.load -> _cuda_init) and the
    // Lightning predict step (accelerator.setup_device -> get_device_capability) abort
    // with "driver too old", leaving a header-only CSV -> empty votes. Terrier is cheap
    // on a consensus library, so it runs on CPU by DEFAULT — independently of the global
    // --use_gpu (the other classifiers' images may tolerate the driver). Hide the GPU so
    // Lightning auto-selects the CPU accelerator. Opt back in only with --terrier_use_gpu.
    def cpu_guard = params.terrier_use_gpu ? '' : 'export CUDA_VISIBLE_DEVICES=""'
    """
    # Activate the container's baked env (read-only .sif cannot `micromamba run`);
    # no-op under the Conda fallback.
    if [ -d /opt/conda/envs/terrier/bin ]; then
        export PATH=/opt/conda/envs/terrier/bin:\$PATH
    fi
    export PYTHONNOUSERSITE=1
    ${cpu_guard}

    WORK="\$PWD"
    export HOME="\$WORK"
    export XDG_CACHE_HOME="\$WORK/.cache"
    export TORCH_HOME="\$WORK/.cache/torch"
    mkdir -p "\$XDG_CACHE_HOME" "\$TORCH_HOME"

    # torchapp writes its weight cache under XDG_CACHE_HOME. The image bakes the
    # weights at /opt/cache/torchapps (read-only), so copy them into the writable
    # cache — no runtime download. Under the Conda fallback the dir is absent and
    # torchapp downloads on first use.
    if [ -d /opt/cache/torchapps ]; then
        cp -r /opt/cache/torchapps "\$XDG_CACHE_HOME/" 2>/dev/null || true
    fi

    # Terrier trips an upstream hierarchicalsoftmax bug in greedy_predictions(): it
    # evaluates predictions[node.softmax_start_index + prediction_child_index] inside an
    # `if threshold and ...` guard, indexing one past the tensor at the last leaf
    # (IndexError "index N out of bounds for dimension 0 with size N"). Terrier's own
    # --threshold does not reach this call, and the image site-packages are read-only
    # (root-owned, so even --writable-tmpfs can't sed them). Instead monkey-patch the
    # function IN MEMORY to force threshold=0 (which short-circuits the buggy comparison
    # before the bad index), then run Terrier's CLI via runpy. Nothing is written to the
    # image. The per-class probability CSV is still produced; normalize_terrier.py
    # --min-prob filters afterward.
    TERRIER_BIN=\$(command -v terrier || true)
    python - "\${TERRIER_BIN}" "${library}" "${meta.id}.terrier.csv" <<'PYEOF' || echo "WARN [${meta.id}]: Terrier run failed; abstaining" >&2
import runpy
import sys
import hierarchicalsoftmax.inference as _inf

_orig_greedy = _inf.greedy_predictions

def _greedy_no_threshold(*args, **kwargs):
    # threshold=0 makes `if threshold and predictions[...]` short-circuit before the
    # out-of-bounds index, so the greedy walk simply descends to a leaf.
    kwargs["threshold"] = 0.0
    return _orig_greedy(*args, **kwargs)

_inf.greedy_predictions = _greedy_no_threshold

# Load the CUDA-saved checkpoint on CPU. Lightning passes torch.load a
# map_location that KEEPS the original device (a `lambda storage, loc: storage`),
# so tensors saved on cuda try storage.to('cuda:0') -> _cuda_init, which aborts
# when the driver is too old for this Torch build (and CUDA_VISIBLE_DEVICES="" does
# not help: _cuda_init probes the driver before checking device visibility). Force
# map_location='cpu' UNCONDITIONALLY so the checkpoint always deserialises on CPU,
# never touching CUDA. Loading on CPU is safe even under --use_gpu (Lightning moves
# the module to the GPU afterwards via the trainer).
import torch
_orig_load = torch.load

def _cpu_load(*args, **kwargs):
    kwargs["map_location"] = "cpu"      # override Lightning's device-preserving map
    return _orig_load(*args, **kwargs)

torch.load = _cpu_load

_bin, _lib, _out = sys.argv[1], sys.argv[2], sys.argv[3]
sys.argv = ["terrier", "--input", _lib, "--output-csv", _out]
runpy.run_path(_bin, run_name="__main__")
PYEOF

    if [[ -s "${meta.id}.terrier.csv" ]]; then
        normalize_terrier.py --pred "${meta.id}.terrier.csv" \\
            --out ${meta.id}.terrier.preds.tsv --min-prob ${minp}
    else
        echo "WARN [${meta.id}]: Terrier produced no CSV; abstaining" >&2
        : > ${meta.id}.terrier.preds.tsv
    fi

    TER_VER=\$(terrier --version 2>/dev/null | tr -d '\\n' || true)
    [[ -z "\${TER_VER}" ]] && TER_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        terrier: "\${TER_VER}"
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.terrier.preds.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        terrier: "stub"
    END_VERSIONS
    """
}
