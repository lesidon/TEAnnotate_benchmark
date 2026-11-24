# nf-core/teannotate

**End-to-end transposable element (TE) and tandem-repeat annotation for animal, plant, and fungal genomes.**

`nf-core/teannotate` is a [Nextflow](https://www.nextflow.io/) DSL2 pipeline that takes one or more genome assemblies and produces a curated, classified TE consensus library, a soft-masked genome, repeat annotations (GFF3 + `.tbl`), and publication-ready visualizations — all reproducibly, on a laptop, a workstation, or an HPC cluster.

---

## Contents

- [Pipeline overview](#pipeline-overview)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Samplesheet](#samplesheet)
- [Profiles (container / conda / scheduler)](#profiles)
- [Key parameters](#key-parameters)
- [Outputs](#outputs)
- [HPC / offline notes](#hpc--offline-notes)
- [Troubleshooting](#troubleshooting)

---

## Pipeline overview

The workflow is organism-aware: it routes each genome through the discovery engines best suited to its clade, then converges on a common curation → classification → masking → visualization path.

| Stage | What happens | Key tools |
|------:|--------------|-----------|
| **0** | Input validation, genome statistics, assembly-quality tiering | `seqkit` |
| **1** | Organism-aware de novo TE discovery. **Animal (default):** Earl Grey + HiTE. **Plant:** EDTA (+ HELIANO) + HiTE. **Fungus:** Earl Grey + HiTE (+ optional Starfish/Starship, which normalizes any GTF/GFF gene annotation to GFF3 via AGAT) | Earl Grey, HiTE, EDTA, HELIANO, Starfish, AGAT |
| **2** | Library deduplication and merge (per-genome consensi + optional user library), then an **evidence gate**. A consensus is kept if it carries a TE protein domain, if another discovery tool found the same element, or if it has ≥ 2 full-length genomic copies; non-TE repeats (satellites, rRNA, …) go to a separate file | CD-HIT-EST, TEsorter, BLAST+ |
| **3** | Curation (boundary trimming, MSA cleanup), with dropped families that carry evidence brought back. Then a DL **re-classification ensemble** over the `Unknown` fraction: a two-step vote picks the Order, then a superfamily inside that Order, each step with its own per-kingdom model weights (switchable to single-step weighted or hard majority voting); by default it may only assign DNA-transposon and Helitron classes | TEtrimmer, TEclass2, CREATE, NeuralTE, Terrier (DeepTE optional) |
| **3-final** | Single-pass genome soft-masking; one RepeatMasker run reports **simple/tandem/satellite AND interspersed** repeats together | RepeatMasker (`dfam/tetools`) |
| **4** | Nature-styled visualization + QC: Kimura divergence landscape, TE-density karyoplot, and a repeat detection/masking QC figure — each ships the exact R script that drew it | R / ggplot2 |

**Organism routing** is chosen per sample from the samplesheet `organism_type` column (`animal` / `plant` / `fungus`); an unspecified type takes the animal route.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) `>= 24.04.0`
- **Conda or Mamba**, with every profile: Earl Grey and four small helper steps (CD-HIT-EST, the Pfam download, and two Python steps of the ensemble) have no container image and always run in a Conda environment the pipeline builds.
- A container engine for everything else: **Singularity / Apptainer** (recommended on HPC; it can pull the DL classifiers' `oras://` images), **Docker**, or **Podman**. Under Docker/Podman the DL classifiers also run from Conda, because those engines cannot pull `oras://` images.
- `-profile conda` / `mamba` alone runs only with `--skip_hite`: HiTE exists only as a container image, and the pipeline stops at launch if it would be needed.

You do **not** need to install the bioinformatics tools yourself: the pipeline pulls the images and builds the Conda environments it needs.

---

## Quick start

```bash
# 1. Create a samplesheet (see below)
# 2. Launch — Singularity backend, local executor
nextflow run . \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results
```

Test the installation with the bundled dataset (*Zymoseptoria tritici* chromosome 21, 409 kb, fungus route; see `tests/test_data/README.md`). A stub run checks the wiring of every step in minutes: every step writes placeholder outputs, and nothing (Pfam, classifier sources or models) is downloaded, but Nextflow still pulls the container images and builds the Conda environments the steps declare. The full test runs every tool on the test genome: about 15 minutes on 4 CPUs once images, models and Conda environments are in place; the first run also downloads several GB of images and classifier models and builds the environments. One chromosome holds few copies of each TE family, so expect only a handful of consensi; the test checks that every step runs, not annotation quality:

```bash
nextflow run . -profile test,singularity -stub --outdir results_stub   # wiring only
nextflow run . -profile test,singularity --outdir results_test         # full pipeline
```

Conda instead of containers (HiTE needs its image, so it is skipped):

```bash
nextflow run . -profile conda --skip_hite --input samplesheet.csv --outdir results
```

On a SLURM cluster (combine a container profile with a scheduler profile):

```bash
nextflow run . -profile singularity,slurm \
    --input samplesheet.csv --outdir results \
    --singularity_cache_dir /scratch/shared/singularity
```

---

## Samplesheet

A CSV with a header row. One row per genome.

```csv
sample,genome,organism_type,gff
speciesA,/data/speciesA.fasta,animal,
speciesB,/data/speciesB.fa.gz,plant,
fungusC,/data/fungusC.fna,fungus,/data/fungusC.braker.gtf
```

| Column | Required | Description |
|--------|:--------:|-------------|
| `sample` | ✔ | Unique sample id (`A-Z a-z 0-9 _ - .`). |
| `genome` | ✔ | Path to the assembly FASTA (`.fa`, `.fasta`, `.fna`, optionally `.gz`): absolute, or relative to the launch directory or to the samplesheet's own directory. Sequence IDs must be unique, and line endings Unix-style. |
| `organism_type` | – | `animal`, `plant`, or `fungus`. Defaults to `animal` when omitted. |
| `gff` | – | Gene annotation (`.gff3`, `.gff`, `.gtf`, optionally `.gz`). Used **only** by the fungal **Starfish** route (`--run_starfish`); ignored on other routes. |

> **Gene-annotation normalization** applies only to the fungal Starfish route: a file already in `.gff3` is used as-is, while `.gtf`/`.gff` are converted to spec-compliant GFF3 (AGAT) before Starfish runs. Genomes with no `gff` skip Starfish.

---

## Profiles

Combine one **backend** profile with, optionally, one **scheduler** profile (e.g. `-profile singularity,slurm`).

**Backend:** `singularity` · `docker` · `podman` · `conda` · `mamba`

**Scheduler:** `slurm` · `pbs` · `lsf` · `sge`

**Other:** `test` (bundled 409 kb test genome) · `debug`

> Under every container profile, Conda is also enabled: steps without an image always use it, and a step whose image is unavailable (offline, unreachable registry, or an `oras://` image under Docker/Podman) falls back to it. You never need to null containers per-process.

---

## Key parameters

Run `nextflow run . --help` for the full, schema-validated list. The most common:

### Input / output
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | – | Samplesheet CSV (required). |
| `--outdir` | `./results` | Output directory. |

### Discovery & curation
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_hite` | `false` | Skip the HiTE discovery engine. |
| `--run_starfish` | `false` | Enable fungal Starship (Starfish) annotation (needs a `gff`; heavy Conda solve). |
| `--user_te_library` | `null` | Optional curated TE consensus FASTA merged into the library. Its sequences are never removed by the evidence gate, and with `--dedup_class_aware` (default) a CD-HIT cluster that contains one keeps the user's sequence as its representative. An ID that is already used by a de novo consensus is renamed `<id>_user`. Which sequences are the user's, and which piece of a split chimera came from which consensus, is recorded in `03_library/<id>.lineage.tsv` / `<id>.refined_lineage.tsv` (`id`, `origin` denovo or user, `parent`, `source_id`, the ID in the input library); the evidence gate reads its trust from there. The user's records are kept as given through curation as well: a consensus TEtrimmer rebuilds from one is dropped in favour of the record, the post-TEtrimmer merge keeps them as representatives, and their labels are only changed when they are `Unknown` (`04_curation/` holds the curated lineage). |
| `--dedup_class_aware` | `true` | After CD-HIT-EST, repair clusters that mix TE classes. CD-HIT keeps each cluster's longest sequence with its own label, so, for example, a LINE collapsed into a chimeric LTR-int consensus is lost and later masked as LTR. The repair relabels Unknown representatives from their classified members and keeps foreign-class members as their own families. Actions are listed in `03_library/<id>.dedup_refine.tsv`. |
| `--dedup_split_chimeras` | `true` | With `--dedup_class_aware`: align the library against itself (blastn). Segments of another TE class (≥500 bp at ≥80% identity) are cut out of the consensus they sit in, when they cover at most half of it. An example is an LTR-int that spans a LINE inserted next to it. The remaining pieces keep the consensus's label. |
| `--dedup_rep_by_copies` | `false` | Experimental. Make the member with the most full-length genomic copies (≥ 75% identity over ≥ 90% of it) each CD-HIT cluster's representative, within the representative's class, instead of the longest one. Off by default because it is harmful in LTR-rich genomes. CD-HIT clusters an LTR terminus with its complete element, and the terminus, having the most copies, replaces it. In the benchmark this cost Drosophila and Zymoseptoria LTR F1 and perfect families. Swaps are listed in `03_library/<id>.dedup_refine.tsv` (action `representative`). |
| `--curation_domain_relabel` | `true` | After TEtrimmer, a curated consensus takes the class of its TE protein domains (TEsorter) when they point to another class than its label. An example is a HELIANO "Helitron" that carries a Gypsy reverse transcriptase. Changes are listed in `04_curation/<id>.domain_relabel.tsv`. Not applied with `--skip_tetrimmer`, where uncurated chimeric consensi carry domains of only part of their sequence. |
| `--curation_class_aware` | `true` | Replaces TEtrimmer's own length-based merge with the class-aware merge above. First, it relabels TEtrimmer consensi that were rebuilt from an element of another class: TEtrimmer keeps the input's label even then. Relabels are listed in `04_curation/<id>.tetrimmer_postprocess.tsv`; merge actions in `04_curation/<id>.curated.dedup_refine.tsv`. |
| `--evidence_gate` | `lite` | Evidence gate after deduplication (`TE_EVIDENCE`). Each consensus is scored on three kinds of evidence: a TE protein domain (TEsorter/REXdb, or a TE-specific Pfam family such as reverse transcriptase, integrase, a DDE transposase or Helitron Rep/Hel on any six-frame ORF fragment, `--evidence_pfam`; Pfam's broad models still catch divergent TEs of unexplored lineages that REXdb's clade profiles miss); the same element in another discovery tool's library (≥ 50% of it at ≥ 80% identity), or TEsorter's homology-only second pass (80-80-80 to a consensus that does carry a domain; column `homology`, counted with this support evidence, never as a domain, and never overriding a label); and genomic copies: ≥ 2 full-length ones (≥ 75% identity over ≥ 90% of it) or ≥ 5 covering half of it (`--evidence_min_half_copies`). Old families decay by mutation and deletion and keep only partial copies; RIP, the fungal defence, makes this extreme. In unexplored genomes such families also lack a known protein domain and a second tool to vouch for them: in the sponge *Oscarella lobularis*, 76 such multi-copy families would otherwise be dropped. `lite` sends non-TE repeats (Satellite, rRNA, tRNA, simple repeats, …) to `03_library/<id>.non_te_repeats.fa` and drops consensi with none of the three. `strict` also drops consensi whose only evidence is full-length copies. `off` skips the step. A union of discovery tools carries every tool's false positives, and most of them sit in consensi with no such evidence. The per-consensus table is `03_library/<id>.evidence.tsv`. Sequences of `--user_te_library` are never gated. |
| `--evidence_relabel_unknown` | `true` | An `Unknown` consensus with a TE protein domain takes the domain's class. Multi-class ("mixture") calls are not used. |
| `--tetrimmer_reinject` | `true` | TEtrimmer writes back a skipped input only if it has ≥ 2 full-length copies, so it discards real low-copy and degraded families. With the evidence table, any input left without a consensus is brought back unchanged if it carries a TE domain or another tool's support. |
| `--pfam_db` | `null` | Pre-downloaded, `hmmpress`-indexed Pfam-A dir (TEtrimmer requires Pfam; otherwise fetched + cached automatically). |
| `--hite_ltrdetector_dir` | `null` | Prebuilt copy of HiTE's `LtrDetector` directory, for hosts that cannot reach GitHub. `kanghu/hite:3.3.3` lacks the binary; by default the pipeline fetches the directory from HiTE commit `--hite_revision` (cached in `--hite_cache`). See the troubleshooting entry below. |
| `--tetrimmer_teaid_dir` | `null` | Host dir with a writable copy of the image's `TE-Aid-master`, bind-mounted over it. Needed only where Singularity has no overlay support — see the troubleshooting entry below. |
| `--skip_tetrimmer` | `false` | Ablation: skip TEtrimmer and re-classify the deduplicated de novo library directly. To run only the library stages on an existing raw library, use `run_reclassify.nf`, which avoids re-running de novo discovery: `nextflow run run_reclassify.nf -profile singularity --input <sample>.raw_library.fa --dedup --genome <genome.fa> --id <sample> --organism_type plant [--provenance <sample>.provenance.tsv] [--skip_tetrimmer] --outdir <dir>`. It publishes the final library as `05_classification/<id>.te_library.fa`, like the full pipeline. Without `--dedup` the library's IDs (before `#`) must be unique; the run stops at launch otherwise. |

### Re-classification ensemble
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_reclassification` | `false` | Skip the whole DL ensemble. |
| `--ensemble_method` | `two_step` | `two_step`: an Order vote, then a superfamily vote **inside the winning Order** (only models whose call lies in that Order take part; DIRS/Ngaro/VIPER count as DIRS, Penelope as PLE), each step with its own per-kingdom model weights and gate. The settings were tuned for support-weighted F1 on the cluster-split train set of leakage-filtered curated consensi (PanTEon + MycoMobilome, held out from RepBase/Dfam); see `TWO_STEP_DEFAULTS` in `bin/ensemble_vote.py`. `weighted` (single Order vote with per-model weights) and `hard` (plain majority) are the older rules; `all` writes every rule as its own voting-table column and relabels with the two-step call. |
| `--ensemble_relabel_orders` | `DNA,RC` | Orders the ensemble may assign; a call of another Order leaves the consensus `Unknown`. The `Unknown` fraction reaching the ensemble has already been through domain-based labelling, so it is mostly consensi without a detectable coding domain — plausible for non-autonomous DNA transposons (MITEs, other TIR elements) and Helitrons, whereas `LTR`/`LINE` calls on such consensi fell mostly outside reference elements of that class on the development genomes. `all` lifts the restriction. |
| `--ensemble_two_step_config` | `null` | Path to a JSON file overriding the two-step settings, with any subset of kingdoms, steps and keys; for example a file containing `{"fungus": {"order": {"weights": {"neuralte": 1, "terrier": 2}}}}`. |
| `--ensemble_min_frac` | `0.5` | Weighted / hard only. Abstention gate: the winning class must carry this fraction of the vote weight the kingdom's **voting subset** could maximally cast; otherwise the consensus stays `Unknown`. Scales with the subset, so one value fits both a 5-model set (`0.5` ⇒ 3-of-5 votes) and a single-model subset. Abstentions count against the winner; a crashed model does not. `0` = no gate (plain argmax). |
| `--ensemble_min_votes` | `0` | Hard voting: optional **absolute** floor on agreeing votes, on top of `--ensemble_min_frac`. `0` = off. |
| `--ensemble_min_agree` | `2` | Weighted / hard only. A relabel needs at least this many models of the voting subset on the winning Order. Unlike `--ensemble_min_votes`, it is never relaxed for a small subset, so a single-model subset relabels nothing. `0` = off. |
| `--ensemble_min_weight` | `0` | Weighted voting: optional **absolute** floor on the winning Order's summed weight, on top of `--ensemble_min_frac`. `0` = off. Prefer `--ensemble_min_frac`: an absolute floor is subset-blind, so a value tuned for five models can be unreachable for a one-model kingdom subset. A floor the subset cannot reach is ignored with a warning. |
| `--ensemble_weights` | `null` | Weighted voting: override weights, e.g. `'neuralte=1.5,terrier=1.5,create=1.3,teclass2=1.2,deepte=1.1'`. Unset uses defaults prioritising NeuralTE/Terrier. |
| `--ensemble_weights_from` | `null` | Weighted voting: derive weights from an order-level metrics TSV (`evaluate_classification.py *.all_models.summary_metrics.tsv`), normalised to mean 1.0. |
| `--ensemble_weight_metric` | `macro_f1` | Metric column read from `--ensemble_weights_from`. |
| `--ensemble_kingdom_aware` | `true` | Run and vote with only the model subset per `organism_type`; classifiers outside it are not launched. `false` = all enabled models run everywhere (the two-step vote still uses only models with a two-step weight). |
| `--ensemble_models_animal` | `teclass2,neuralte,terrier` | Model subset that runs and votes for animal samples (the two-step animal set). |
| `--ensemble_models_plant` | `teclass2,create,neuralte,terrier` | Model subset that runs and votes for plant samples (the two-step plant set). |
| `--ensemble_models_fungus` | `create,neuralte,terrier` | Model subset that runs and votes for fungal samples (the two-step fungus set). |
| `--skip_teclass2` / `--skip_deepte` / `--skip_create` / `--skip_neuralte` / `--skip_terrier` | `false` | Drop an individual ensemble member. |
| `--use_gpu` | `false` | Run GPU-capable classifiers on a GPU (adds Singularity `--nv`). |

### Masking & visualization
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_masking` | `false` | Skip RepeatMasker soft-masking (and, consequently, visualization). |
| `--skip_visualization` | `false` | Skip Stage 4 figures/QC. |
| `--plot_formats` | `png,svg` | Output image formats (`png`, `svg`, `pdf`). |
| `--karyoplot_window_kb` | `500` | Karyoplot bin size (kb). |
| `--kimura_faceted` | `false` | Facet the Kimura landscape per class instead of one stacked area. |

### Environment / caching
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--use_mamba` | `false` | Solve Conda envs with mamba. |
| `--singularity_cache_dir` | `null` | Shared `.sif` cache (or set `$NXF_SINGULARITY_CACHEDIR`). |
| `--conda_cachedir` | `null` | Shared Conda-env cache. |
| `--pfam_cache` / `--hite_cache` / `--teclass2_cache` / `--create_cache` / `--deepte_cache` / `--neuralte_cache` | `<launchDir>/<name>` | Caches for assets fetched at run time (see the HPC notes). |
| `--auto_bind_ps` | `true` | Auto-detect the host `ps` and bind it into Singularity containers (see below). |

---

## Outputs

Results are written under `--outdir/<sample>/`, in numbered stage folders for intuitive navigation:

```
results/<sample>/
├── 01_genome_stats/       # seqkit stats + assembly-quality JSON
├── 02_denovo_te/          # per-engine raw libraries + merged raw_library.fa
│                          #   (fungus: starfish/ holds the normalized *.gff3)
├── 03_library/            # CD-HIT-EST deduplicated library; evidence.tsv, gated_library.fa,
│                          #   non_te_repeats.fa and gate_removed.fa from the evidence gate;
│                          #   lineage.tsv: origin and parent of every sequence
├── 04_curation/           # TEtrimmer curated library (Need_check consensi removed)
├── 05_classification/     # per-tool predictions, ensemble voting_table.tsv, classified_library.fa,
│                          #   and te_library.fa: the final library, the one masking uses
├── 06_masking/            # softmasked.fa, repeats.{tbl,out,align,bed}, simple/te/combined GFF3
└── 07_visualization/      # Kimura landscape, karyoplot, QC figure + the R script for each
results/pipeline_info/     # software_versions.yml (every tool + Nextflow + pipeline commit),
                           #   asset_manifest.tsv (Pfam, classifier sources/models, LtrDetector:
                           #   cache directory, requested + resolved identity, checksum)
```

Highlights:

- **`05_classification/<sample>.te_library.fa`** — the final TE library on every route (re-classified, curated, or deduplicated when steps are skipped), with RepeatMasker-style `name#Class/Superfamily` headers.

- **`06_masking/`** — the `.tbl` reports **both** simple/satellite **and** interspersed repeats (single-pass RepeatMasker; with an empty TE library it still masks simple and low-complexity repeats). Three GFF3s are produced: simple-only, TE-only, and combined, with the repeat in `Name` and its class in `repeat_class`. The `.align` feeds the Kimura landscape.
- **`05_classification/…voting_table.tsv`** — a model-voting matrix (rows = repeat ids, columns = each classifier + the ensemble decision).
- **`07_visualization/`** — every figure ships alongside the exact R script (`plot_kimura.R`, `plot_karyotype.R`, `repeat_qc.R`) so you can restyle and re-run it.

---

## HPC / offline notes

- **Shared caches:** point `--singularity_cache_dir` / `--conda_cachedir` at shared scratch so images/envs build once and are reused across runs.
- **Fetched assets are pinned and cached per identity:** Pfam-A (release 38.2), the TEclass2 / CREATE / DeepTE / NeuralTE sources (the commits the ensemble was tuned with, `--<tool>_revision`), their model archives and HiTE's LtrDetector. Each lands in a subdirectory of its cache named after the release, commit or URL, so changing one of these parameters downloads the new asset instead of reusing the old one. Each subdirectory holds a `source.tsv` with the resolved commit or the archive's SHA-256. Caches filled by earlier versions of the pipeline (without these subdirectories) are not reused; the assets download once more. A cached asset is reused as found (its files are not re-verified); delete its subdirectory to fetch it again. `pipeline_info/asset_manifest.tsv` lists every asset a run used, including those given as local paths.
- **Compute nodes without internet:** Conda environments are built where the Nextflow driver runs, and some packages download data while they install (the Earl Grey environment fetches an R data package from bioconductor.org). Build the environments once from a machine with internet access into the shared `--conda_cachedir` (for example a `-profile test` run on a login or transfer node); later runs on compute nodes reuse them. Environments use conda-forge and bioconda only, never Anaconda's `defaults` channel.
- **Offline:** the container ⇄ Conda resolver probes the registry once; force Conda with `-profile conda`, or set `--assume_offline true`. Supply local `.sif` images via `--<tool>_container /abs/path/img.sif`.
- **`ps` metric collector:** Nextflow calls `ps` to collect task metrics, and on some clusters a container that lacks a working `ps` aborts the task. `--auto_bind_ps` (default **on**) auto-detects the host `ps` and its `procps` library and binds them into every Singularity/Apptainer container. Under a remote scheduler where the driver and compute nodes differ, either keep classifier images (which ship `procps-ng`) or pass `--host_ps <path>`; disable entirely with `--auto_bind_ps false`.

---

## Troubleshooting

- **A classifier step fails but the run continues** — by design: each DL classifier is `errorStrategy 'ignore'`; if a member fails or abstains, the ensemble decides from the remaining votes and curation falls back to the curated (non-reclassified) library.
- **`calcDivergenceFromAlign.pl` errors / empty Kimura matrix** — the RepeatMasker steps use `dfam/tetools`, and the pipeline sets `PERL5LIB` so the utility finds its modules. If the divergence step still yields nothing, the landscape falls back to binning the RepeatMasker `.out` divergence column, so a plot is still produced whenever repeats exist.
- **TEtrimmer / Pfam** — TEtrimmer requires Pfam-A and aborts without it; the pipeline fetches and indexes it once (cached), or pass `--pfam_db` to a pre-indexed directory.
- **TEtrimmer `OSError: [Errno 30] Read-only file system: '.../TE-Aid-master'`** — TEtrimmer `chmod`s its bundled TE-Aid at startup. The `--writable-tmpfs` the pipeline already passes fixes this *only where Singularity has overlay support*; an unprivileged Singularity with no setuid path uses **underlay** mode and drops the flag (`WARNING: Disabling --writable-tmpfs as it can't be used in conjunction with underlay`), leaving the image read-only. Bind a writable copy over that path instead — populate it once from the image, then pass `--tetrimmer_teaid_dir`:
  ```bash
  IMG=<NXF_SINGULARITY_CACHEDIR>/quay.io-biocontainers-tetrimmer-1.7.2--hdfd78af_0.img
  mkdir -p "$HOME/containers/teaid_rw"
  singularity exec -B "$HOME/containers":/mnt "$IMG" \
      cp -r /usr/local/share/tetrimmer/TE-Aid-master/. /mnt/teaid_rw/
  chmod -R u+rwX "$HOME/containers/teaid_rw"
  nextflow run . ... --tetrimmer_teaid_dir "$HOME/containers/teaid_rw"
  ```
  The same missing-overlay limitation is why `--writable-tmpfs` cannot rescue container-based EarlGrey on such clusters.
- **HiTE: "LtrDetector does not run in this image"** — `kanghu/hite:3.3.3` ships without `bin/FiLTR-main/bin/LtrDetector/bin/`, the `LtrDetector` binary that HiTE's default LTR step (FiLTR) runs. HiTE does not check that step's exit status, so it would finish "successfully" with an LTR-less library. The pipeline therefore fetches that directory from the HiTE commit the image code matches (`HITE_LTRDETECTOR`, about 1.5 MB, cached), runs HiTE from a linked copy of `/HiTE` that uses it, and checks that `LtrDetector` runs before discovery starts. Where the fetch cannot reach GitHub, build the directory once elsewhere and pass it with `--hite_ltrdetector_dir`:
  ```bash
  git clone https://github.com/CSU-KangHu/HiTE.git hite_src && git -C hite_src checkout 03c2e77
  D=$PWD/hite_src/bin/FiLTR-main/bin/LtrDetector
  # git stores the binary non-executable, and it needs libgomp from the image's conda env
  mv $D/bin/LtrDetector $D/bin/LtrDetector.bin && chmod +x $D/bin/LtrDetector.bin
  printf '#!/bin/sh\nLD_LIBRARY_PATH=/opt/conda/envs/HiTE/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} exec "$(dirname "$0")/LtrDetector.bin" "$@"\n' > $D/bin/LtrDetector
  chmod +x $D/bin/LtrDetector
  nextflow run . ... --hite_ltrdetector_dir "$D"
  ```
  Do not switch FiLTR off instead (`--use_FiLTR 0`): HiTE's classic LTR_harvest/LTR_finder → LTR_retriever path fails in this image too (`The RMblast engine is not installed in RepeatMasker!`).
- **Masking produced an empty annotation / `createLib(): Error invoking … makeblastdb … general.working/is.lib`** — RepeatMasker builds its general libraries inside its install directory, which in the container lives in the `--writable-tmpfs` overlay. On hosts where Apptainer falls back to fuse-overlayfs (`df -h /` inside the container shows `fuse`), makeblastdb cannot write there. The masking step now runs RepeatMasker with `-libdir` pointing at a task-local copy (symlinks) of its `Libraries`, so this build happens in the work directory on every host. The step also exits with an error, and never passes the genome through unmasked, if RepeatMasker fails; Nextflow retries it up to twice.
- **EDTA: "TIR results not found" / rice TEs in the library** — The TIR-Learner 3.0.1 bundled with `edta:2.3.0` crashes under the image's pandas 3 in two places:
  - `KeyError: 0` in `check_TIR_TSD.py`;
  - `KeyError: '<seqid>_split_NofM'` for chromosome IDs containing dots, such as NCBI accessions.

  EDTA then skips HelitronScanner. With `--force 1` it silently substitutes its bundled **rice** TIR and Helitron library (sequences named `Os00xx`). The EDTA step now runs a copy of TIR-Learner patched at run time, and fails with an error if any EDTA module still falls back to the rice library.
- **EDTA: whole multi-copy (mostly LTR) families missing from the library, differently from run to run** — EDTA's redundancy filter `cleanup_nested.pl` runs its comparisons in threads. When two near-identical copies of a family are processed at the same time, each removes the other, and the family disappears from `*.EDTA.TElib.fa`. In *Zymoseptoria tritici*, a 13 kb Gypsy and a 6 kb Copia internal region, about 1 Mb of genomic LTR, vanished this way. The EDTA step now runs EDTA from a copy in the work directory (its `database/` is linked, not copied). In that copy, `cleanup_nested.pl` applies each query's edits under a lock and skips a query that another thread has already changed or removed. This gives the single-threaded result while keeping the threads.
- **TEtrimmer drops low-copy families ("check_low_copy is False") / TE-Aid "BLAST Database error … RepeatPeps.lib"** — TE-Aid downloads `RepeatPeps.lib` from a URL that now returns 404, so every TE-Aid call fails. TEtrimmer's low-copy check depends on TE-Aid, so it discards nearly every family with fewer than 10 genomic hits. The TEtrimmer step now runs TEtrimmer from a copy of its own source in the work directory. The copy's TE-Aid gets a `db/` with the `RepeatPeps.lib` shipped with the image's RepeatMasker. Because the copy is writable, this works whether or not the container image is writable.

---

## Citation

This pipeline follows the [nf-core](https://nf-co.re) template and conventions. Please also cite the underlying tools (Earl Grey, HiTE, EDTA, HELIANO, Starfish, CD-HIT, TEtrimmer, TEclass2, DeepTE, CREATE, NeuralTE, Terrier, RepeatMasker, AGAT, seqkit) — versions are recorded per run in `pipeline_info/`.
