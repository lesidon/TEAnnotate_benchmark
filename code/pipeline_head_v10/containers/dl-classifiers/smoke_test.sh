#!/usr/bin/env bash
# Smoke-test each per-tool DL classifier .sif: confirm its conda env imports and the
# tool actually runs on a tiny Unknown library, dumping each tool's RAW output so its
# on-disk format is visible. The NeuralTE / Terrier dumps are what we need to lock
# their ensemble vote normalizers (normalize_neuralte.py / normalize_terrier.py).
#
# Run from a WRITABLE directory (outputs land under ./_smoke, bind-mounted to /work
# inside each image). SIF_DIR defaults to the current dir.
#
#   ./containers/dl-classifiers/smoke_test.sh [SIF_DIR]
#
# Each image ships the tool in a conda env at /opt/conda/envs/<tool>; a read-only
# .sif cannot use `micromamba run`, so we activate by putting that env's bin on PATH
# (exactly what the pipeline modules do).
set -uo pipefail        # deliberately NOT -e: run every tool even if one fails

SIF_DIR="${1:-$PWD}"
WORK="$PWD/_smoke"
mkdir -p "$WORK"

# Tiny 2-sequence Unknown library (random 600 bp; a smoke test, not a benchmark).
python3 - > "$WORK/test.fa" <<'PY'
import random
random.seed(1)
for i in (1, 2):
    print(f">rnd_1_family_{i}#Unknown")
    print(''.join(random.choice('ACGT') for _ in range(600)))
PY
echo "test library: $WORK/test.fa"

# CREATE reads each sequence's TRUE label from the part BEFORE '#' (a benchmark
# convention) and crashes if it isn't a taxonomy node. Give the valid root node "TE"
# and move the real family id AFTER '#', which survives into the relabeled headers so
# predictions can be mapped back. This mirrors what the CREATE module will do.
python3 - "$WORK/test.fa" "$WORK/create_in.fa" <<'PY'
import sys
inp, out = sys.argv[1], sys.argv[2]
with open(inp) as fh, open(out, 'w') as o:
    for line in fh:
        if line.startswith('>'):
            base = line[1:].split()[0].split('#')[0]
            o.write('>TE#%s\n' % base)
        else:
            o.write(line)
PY

# Bind $WORK to /work and expose it as $WORK inside the container, so the inner
# (single-quoted) scripts never interpolate host paths.
#   PYTHONNOUSERSITE=1  stops the host's ~/.local/lib/pythonX.Y/site-packages from
#                       leaking in (Apptainer mounts $HOME) and shadowing the image's
#                       pinned numpy/etc. — this is what broke NeuralTE.
#   HOME=/work          a writable HOME whose (absent) .local can't shadow anything.
# The tool caches (XDG_CACHE_HOME/TORCH_HOME) are redirected to writable /work below,
# since /opt/cache in the .sif is read-only (this is what broke Terrier).
run() {
    local name="$1" sif="$2"; shift 2
    echo; echo "############################ ${name} ############################"
    if [ ! -f "$sif" ]; then echo "SKIP: $sif not found"; return; fi
    apptainer exec --bind "$WORK":/work \
        --env WORK=/work --env HOME=/work --env PYTHONNOUSERSITE=1 \
        --env XDG_CACHE_HOME=/work/cache --env TORCH_HOME=/work/cache/torch \
        "$sif" "$@"
}
mkdir -p "$WORK/cache/torch"

# ── TEclass2 (Transformer) ───────────────────────────────────────────────────
run teclass2 "$SIF_DIR/teannotate-teclass2.sif" bash -lc '
  export PATH=/opt/conda/envs/teclass2/bin:$PATH
  python -c "import torch, transformers, numpy; print(\"import OK  torch\", torch.__version__)"
  python /opt/teclass2/TEclass2.py --help 2>&1 | head -20'

# ── CREATE (attention CNN+RNN; Wicker) ───────────────────────────────────────
run create "$SIF_DIR/teannotate-create.sif" bash -lc '
  export PATH=/opt/conda/envs/create/bin:$PATH
  python -c "import tensorflow, pandas, numpy; print(\"import OK  tf\", tensorflow.__version__)"
  MDIR=/opt/models/create; [ -d "$MDIR/model" ] && MDIR="$MDIR/model"
  cd /opt/create
  python create_test.py -i "$WORK/create_in.fa" -d "$MDIR" -o "$WORK/create_out" -p 0.5 || true
  echo "---- relabeled_seqs.fasta headers (label <TAB> orig) ----"
  find "$WORK/create_out" -name relabeled_seqs.fasta -exec grep "^>" {} \; | head'

# ── DeepTE (models are fetched by the pipeline at runtime; only --help here) ──
run deepte "$SIF_DIR/teannotate-deepte.sif" bash -lc '
  export PATH=/opt/conda/envs/deepte/bin:$PATH
  python -c "import tensorflow, keras, numpy; print(\"import OK  tf\", tensorflow.__version__)"
  python /opt/deepte/DeepTE.py --help 2>&1 | head -20'

# ── NeuralTE (RepeatMasker labels via --is_wicker 0) ─────────────────────────
run neuralte "$SIF_DIR/teannotate-neuralte.sif" bash -lc '
  export PATH=/opt/conda/envs/neuralte/bin:$PATH
  python -c "import tensorflow, keras, Bio, sklearn; print(\"import OK  tf\", tensorflow.__version__)"
  mkdir -p "$WORK/nte_out"; cd /opt/neuralte
  python src/Classifier.py --data "$WORK/test.fa" \
      --model_path /opt/neuralte/models/NeuralTE_model.h5 \
      --outdir "$WORK/nte_out" --is_predict 1 --is_wicker 0 --thread 4 || true
  echo "---- NeuralTE output files ----"; find "$WORK/nte_out" -type f
  echo "---- first lines of each ----"
  find "$WORK/nte_out" -type f \( -name "*.info" -o -name "*.tsv" -o -name "*.txt" -o -name "*.fasta" \) \
      -exec sh -c "echo; echo == \$1 ==; head -5 \$1" _ {} \;'

# ── Terrier (RepeatMasker-schema CSV) ────────────────────────────────────────
run terrier "$SIF_DIR/teannotate-terrier.sif" bash -lc '
  export PATH=/opt/conda/envs/terrier/bin:$PATH
  terrier --input "$WORK/test.fa" --output-csv "$WORK/terrier_out.csv" || true
  echo "---- terrier_out.csv ----"; head -5 "$WORK/terrier_out.csv"'

echo; echo "Done. Outputs under: $WORK"
echo "Please paste the NeuralTE output block and the terrier_out.csv header so the"
echo "two vote normalizers can be finalized."
