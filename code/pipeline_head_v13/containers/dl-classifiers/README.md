# DL TE-classifier container (TEclass2 + CREATE + DeepTE + NeuralTE + Terrier)

One image, five conda envs — built so the network-restricted HPC compute node
needs **zero internet** for the Stage 3 re-classification ensemble. Replaces the
per-tool conda builds (which repeatedly failed on Zhores: solver kills, mamba
`Multi-download failed`, missing pip deps).

## Why one image with five envs

The five classifiers cannot share a Python environment:

| env (`-n`) | Python | Framework | source | model |
|---|---|---|---|---|
| `teclass2` | 3.10 | PyTorch `<1.13`, transformers 4.17 | `/opt/teclass2` | `/opt/models/teclass2` |
| `create`   | 3.6.13 | TensorFlow 2.4 | `/opt/create` | `/opt/models/create` |
| `deepte`   | 3.7 | TensorFlow 2.2 / Keras 2.3.1 (protobuf<3.20) | `/opt/deepte` | `/opt/models/deepte` |
| `neuralte` | 3.8 | TensorFlow / Keras 2.6 | `/opt/neuralte` | `/opt/neuralte/models/NeuralTE_model.h5` (in-repo) |
| `terrier`  | 3.11 | PyTorch (`bio-terrier`) | pip package | weights pre-cached at build (`/opt/cache`) |

Each tool is invoked with `micromamba run -n <env> …`. NeuralTE outputs RepeatMasker
labels with `--is_wicker 0`; Terrier outputs RepeatMasker-schema labels natively.

## Two packagings

You can build **one combined image** (all five envs) or **five per-tool images**.
They expose the **same** env names (`micromamba run -n <env>`) and `TEANNOTATE_*`
paths, so the pipeline works with either. Build everything **from the repository
root** (the `%files`/COPY paths resolve from there); it needs internet and clones the
repos + downloads the models.

### Recommended on a flaky network — per-tool images (Apptainer)

A monolithic build is fragile: a single stalled conda download aborts 20+ minutes of
work. Five small independent builds isolate that — re-run only the one that stalled.

```bash
# from the repo root, on the login node (has internet):
bash containers/dl-classifiers/build_all.sh
# -> teannotate-teclass2.sif, teannotate-create.sif, teannotate-deepte.sif,
#    teannotate-neuralte.sif, teannotate-terrier.sif
```

Per-tool defs: [`teclass2.def`](teclass2.def), [`create.def`](create.def),
[`deepte.def`](deepte.def), [`neuralte.def`](neuralte.def), [`terrier.def`](terrier.def).
Build one by hand with
`apptainer build --fakeroot teannotate-deepte.sif containers/dl-classifiers/deepte.def`.
Each `.def` sets conda retry/backoff to ride out brief stalls.

### Publishing (Singularity Hub is retired)

Singularity Hub shut down in 2021. Publish so others can pull:

```bash
# Sylabs Cloud Library (closest successor; free): apptainer remote login first
apptainer push teannotate-deepte.sif library://<user>/teannotate/deepte:1.0
#   -> pull/use as  library://<user>/teannotate/deepte:1.0

# or GitHub Container Registry via ORAS
apptainer push teannotate-deepte.sif oras://ghcr.io/<user>/teannotate-deepte:1.0
```

Or simply keep the local `.sif` files and point the pipeline at their absolute paths.

---

### Combined single image

Both recipes below build the all-in-one image. Use them if you prefer one artifact.
The build clones four repos and downloads all pretrained models (~15–30 min).

### Option A — Apptainer/Singularity, directly on the HPC (no Docker)

Use [`teannotate-dl.def`](teannotate-dl.def). Put the build scratch on a **writable**
filesystem (not the read-only `/trinity` home):

```bash
export APPTAINER_TMPDIR=/beegfs/home/$USER/tmp/apptainer_tmp
export APPTAINER_CACHEDIR=/beegfs/home/$USER/tmp/apptainer_cache
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
apptainer build --fakeroot teannotate-dl.sif containers/dl-classifiers/teannotate-dl.def
```

If `--fakeroot` is not permitted: try plain `apptainer build` (newer Apptainer falls
back to `proot`), or `apptainer build --remote …` (cloud build, needs a free Sylabs
token), or ask an admin to enable fakeroot (`subuid`/`subgid`). `singularity build`
takes the same arguments. The result is a local `teannotate-dl.sif`.

### Option B — Docker, then push to a registry

```bash
docker build -f containers/dl-classifiers/Dockerfile -t <dockerhub-user>/teannotate-dl:1.0 .
docker push <dockerhub-user>/teannotate-dl:1.0
```

Pin tool revisions in either recipe with `--build-arg TECLASS2_REV=… CREATE_REV=…
DEEPTE_REV=… NEURALTE_REV=…` (Docker) or by editing the `git clone` lines (`.def`).

## Use it in the pipeline

Once the image exists, tell me its reference and I'll wire it into the five
classifier modules (a `container` directive with the conda envs kept as the offline
fallback, plus the `*_CLONE` / `*_FETCH_MODEL` steps skipped when the image provides
the baked-in `/opt/...` paths), add the NeuralTE + Terrier modules, and set
`ensemble_min_votes = 3`. Reference either a local `.sif`
(`--dl_container /abs/path/teannotate-dl.sif`) or a registry image
(`--dl_container docker://<dockerhub-user>/teannotate-dl:1.0`).

## Smaller image (optional)

To fetch models at runtime instead of baking them in (image ~a few GB smaller),
delete **section 3** of the `Dockerfile` and the `TEANNOTATE_*_MODEL` `ENV` lines;
the pipeline's `*_FETCH_MODEL` steps then download them (those fetches currently
work on the node and are storeDir-cached).
