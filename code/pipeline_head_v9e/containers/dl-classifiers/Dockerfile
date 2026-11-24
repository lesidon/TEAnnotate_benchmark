# DL TE-classifier image for nf-core/teannotate:
#   TEclass2 + CREATE + DeepTE + NeuralTE + Terrier.
#
# These five tools have mutually-incompatible stacks (Python 3.6/3.7/3.8/3.10/3.11,
# PyTorch / TF 2.2 / 2.4 / 2.6), so they live in FIVE separate conda envs in one
# image and are each invoked via `micromamba run -n <env> ...`. Source repos and
# pretrained models are baked in, so the (network-restricted) HPC compute node needs
# ZERO internet at runtime — no conda solve, no pip download, no git clone, no
# Google-Drive/cyVerse model fetch, no Terrier weight download.
#
# BUILD FROM THE REPO ROOT (the env ymls must be in the build context):
#   docker build -f containers/dl-classifiers/Dockerfile -t <dockerhub-user>/teannotate-dl:1.0 .
#   docker push <dockerhub-user>/teannotate-dl:1.0
#
# The pipeline pulls it via Singularity (docker://<dockerhub-user>/teannotate-dl:1.0).

FROM mambaorg/micromamba:1.5.8

# Build-time tools to clone repos and download models (not needed at runtime).
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        git curl ca-certificates tar gzip procps && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Fixed cache so model weights downloaded at build time are found at runtime too.
ENV XDG_CACHE_HOME=/opt/cache \
    TORCH_HOME=/opt/cache/torch \
    HF_HOME=/opt/cache/hf \
    MAMBA_ROOT_PREFIX=/opt/conda

# ── 1. Conda envs (one per tool) ──────────────────────────────────────────────
COPY --chown=$MAMBA_USER:$MAMBA_USER assets/environments/teclass2.yml /tmp/teclass2.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER assets/environments/create.yml   /tmp/create.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER assets/environments/deepte.yml   /tmp/deepte.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER assets/environments/neuralte.yml /tmp/neuralte.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER assets/environments/terrier.yml  /tmp/terrier.yml
RUN micromamba create -y -n teclass2 -f /tmp/teclass2.yml && \
    micromamba create -y -n create   -f /tmp/create.yml   && \
    micromamba create -y -n deepte   -f /tmp/deepte.yml   && \
    micromamba create -y -n neuralte -f /tmp/neuralte.yml && \
    micromamba create -y -n terrier  -f /tmp/terrier.yml  && \
    micromamba clean -ay

# Install each env's pip deps EXPLICITLY (do not rely on micromamba processing the
# yml `pip:` subsection — that silent gap is what left numpy/pandas missing on the
# cluster). protobuf<3.20 is required for DeepTE's TF 2.2 (its *_pb2 code breaks with
# protobuf 4.x). gdown is for the CREATE model below.
RUN micromamba run -n teclass2 pip install --no-cache-dir \
        "torch<1.13" torchvision tokenizers "transformers==4.17" accelerate && \
    micromamba run -n create   pip install --no-cache-dir \
        "tensorflow-gpu==2.4.0" gdown && \
    micromamba run -n deepte   pip install --no-cache-dir \
        "numpy==1.18.5" "tensorflow==2.2.0" "keras==2.3.1" \
        "scikit-learn==0.24.2" "biopython==1.79" "protobuf==3.19.6" && \
    micromamba run -n neuralte pip install --no-cache-dir \
        "numpy==1.19.5" "tensorflow==2.6.0" "keras==2.6.0" "protobuf==3.19.6" \
        "h5py==3.1.0" "scikit-learn==1.2.0" "biopython==1.78" \
        "seaborn==0.11.2" "matplotlib==3.5.3" openpyxl && \
    micromamba run -n terrier  pip install --no-cache-dir bio-terrier "click<8.2"

# Fail the build NOW if anything is missing (so we never ship a half-populated env).
RUN micromamba run -n teclass2 python -c "import torch, transformers, numpy" && \
    micromamba run -n create   python -c "import pandas, numpy, tensorflow" && \
    micromamba run -n deepte   python -c "import numpy, tensorflow, keras" && \
    micromamba run -n neuralte python -c "import tensorflow, keras, Bio, sklearn" && \
    micromamba run -n terrier  terrier --help >/dev/null

# ── 2. Tool source ────────────────────────────────────────────────────────────
USER root
ARG TECLASS2_REV=main
ARG CREATE_REV=master
ARG DEEPTE_REV=master
ARG NEURALTE_REV=master
RUN git clone https://github.com/IOB-Muenster/TEclass2.git /opt/teclass2 && \
    git -C /opt/teclass2 checkout "$TECLASS2_REV" && \
    git clone https://github.com/yangqi-cs/CREATE.git /opt/create && \
    git -C /opt/create checkout "$CREATE_REV" && \
    git clone https://github.com/LiLabAtVT/DeepTE.git /opt/deepte && \
    git -C /opt/deepte checkout "$DEEPTE_REV" && \
    git clone https://github.com/CSU-KangHu/NeuralTE.git /opt/neuralte && \
    git -C /opt/neuralte checkout "$NEURALTE_REV"
# (Terrier has no repo to clone — it is the pip package; its weights are fetched below.)

# ── 3. Pretrained models ──────────────────────────────────────────────────────
# NeuralTE's model ships in its repo (/opt/neuralte/models/NeuralTE_model.h5).
# DeepTE's models are NOT baked: they are segmented by organism group and the
# pipeline's DEEPTE_FETCH_MODEL pulls only the needed group at runtime into the host
# cache (bind-mounted in) — keeps the image small and the build fast. TEclass2/CREATE
# downloads can likewise be dropped in favour of their *_FETCH_MODEL steps.
RUN mkdir -p /opt/models/teclass2 /opt/models/create /opt/cache && \
    # TEclass2 model
    curl -fL "https://bioinformatics.uni-muenster.de/tools/teclass2/models/teclass2-models.tgz" \
        -o /tmp/tc2.tgz && tar -xzf /tmp/tc2.tgz -C /opt/models/teclass2 && rm /tmp/tc2.tgz && \
    # CREATE model (Google Drive)
    cd /opt/models/create && \
    micromamba run -n create gdown "https://drive.google.com/uc?id=1OWzmD5sMM7kAJjKmq0LS3HibqzNVFThe" \
        -O model.tar.gz && tar -xzf model.tar.gz && rm model.tar.gz

# Terrier weights auto-download on first run — trigger that now so they are cached in
# the image (under $XDG_CACHE_HOME/$TORCH_HOME). The seed run is best-effort: the
# weight download happens at model-load, before inference, so `|| true` is safe.
RUN micromamba run -n create python -c "import random; random.seed(0); \
print('>seed'); print(''.join(random.choice('ACGT') for _ in range(600)))" > /tmp/seed.fa && \
    (micromamba run -n terrier terrier --input /tmp/seed.fa --output-csv /tmp/seed.csv || true) && \
    rm -f /tmp/seed.fa /tmp/seed.csv && \
    chmod -R a+rX /opt/teclass2 /opt/create /opt/deepte /opt/neuralte /opt/models /opt/cache

# Fixed in-image locations the pipeline modules read (TEANNOTATE_* env vars).
ENV TEANNOTATE_TECLASS2=/opt/teclass2 \
    TEANNOTATE_CREATE=/opt/create \
    TEANNOTATE_DEEPTE=/opt/deepte \
    TEANNOTATE_NEURALTE=/opt/neuralte \
    TEANNOTATE_TECLASS2_MODEL=/opt/models/teclass2 \
    TEANNOTATE_CREATE_MODEL=/opt/models/create \
    TEANNOTATE_NEURALTE_MODEL=/opt/neuralte/models/NeuralTE_model.h5
USER $MAMBA_USER

LABEL org.opencontainers.image.title="teannotate-dl" \
      org.opencontainers.image.description="TEclass2 + CREATE + DeepTE + NeuralTE + Terrier TE classifiers (5 conda envs) for nf-core/teannotate"
