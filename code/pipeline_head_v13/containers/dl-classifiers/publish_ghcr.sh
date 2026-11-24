#!/usr/bin/env bash
# Publish the per-tool DL classifier .sif images to the GitHub Container Registry
# (GHCR) as ORAS artifacts, so the pipeline — and anyone else — can pull them with
# Apptainer/Singularity via  oras://ghcr.io/<owner>/teannotate-<tool>:<tag>.
#
# GHCR is free and public, ties to the same GitHub account as the repo, and needs
# no separate service (unlike the retired Singularity Hub). Apptainer speaks OCI/ORAS
# natively, so a .sif pushes and pulls with no conversion.
#
# Prerequisites
#   * Apptainer (auth is `apptainer remote login` — there is no `registry` subcommand).
#   * A GitHub Personal Access Token (classic) with `write:packages` + `read:packages`
#     (or a fine-grained token with Packages: Read and write). Create at
#     https://github.com/settings/tokens .
#
# Usage
#   export GHCR_USER=lesidon           # your GitHub username or org
#   export GHCR_TOKEN=ghp_xxxxxxxx     # PAT with write:packages
#   ./containers/dl-classifiers/publish_ghcr.sh [TAG] [SIF_DIR]
#       TAG      image tag       (default: 1.0)
#       SIF_DIR  dir with .sifs  (default: current dir)
set -euo pipefail

TAG="${1:-1.0}"
SIF_DIR="${2:-$PWD}"
: "${GHCR_USER:?set GHCR_USER to your GitHub username/org}"
: "${GHCR_TOKEN:?set GHCR_TOKEN to a PAT with write:packages}"

TOOLS=(teclass2 create deepte neuralte terrier)

# One-time login. Credentials are stored per registry host, so logging in to
# oras://ghcr.io authorises the oras:// pushes below. If your Apptainer is too old
# for --password-stdin, use:  apptainer remote login -u "$GHCR_USER" -p "$GHCR_TOKEN" oras://ghcr.io
echo "$GHCR_TOKEN" | apptainer remote login --username "$GHCR_USER" --password-stdin oras://ghcr.io

for tool in "${TOOLS[@]}"; do
    sif="${SIF_DIR}/teannotate-${tool}.sif"
    ref="oras://ghcr.io/${GHCR_USER}/teannotate-${tool}:${TAG}"
    if [ ! -f "$sif" ]; then echo "SKIP: $sif not found"; continue; fi
    echo ">> pushing ${sif}  ->  ${ref}"
    apptainer push "$sif" "$ref"
done

cat <<EOF

── Next (one-time) ──────────────────────────────────────────────────────────────
Make each package PUBLIC in GitHub (packages start private):
  https://github.com/${GHCR_USER}?tab=packages
  -> each teannotate-* package -> Package settings -> Change visibility -> Public

── Pull test ────────────────────────────────────────────────────────────────────
  apptainer pull oras://ghcr.io/${GHCR_USER}/teannotate-teclass2:${TAG}

── Wire into the pipeline (per tool) ─────────────────────────────────────────────
  --teclass2_container oras://ghcr.io/${GHCR_USER}/teannotate-teclass2:${TAG}
  --create_container   oras://ghcr.io/${GHCR_USER}/teannotate-create:${TAG}
  --deepte_container   oras://ghcr.io/${GHCR_USER}/teannotate-deepte:${TAG}
  --neuralte_container oras://ghcr.io/${GHCR_USER}/teannotate-neuralte:${TAG}
  --terrier_container  oras://ghcr.io/${GHCR_USER}/teannotate-terrier:${TAG}
Point --singularity_cache_dir at a shared dir so each image pulls once and is reused.
EOF
