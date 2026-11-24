# TEAnnotate benchmark

This source tree accompanies the TEAnnotate manuscript. It covers the 15-genome library benchmark, the maize large-library run, the HELIANO and ensemble comparisons, the classifier-only test sets, and the manuscript figures. The [TEAnnotate repository](https://github.com/lesidon/TEAnnotate) contains the current workflow for new analyses.

## Layout

| Directory | Contents |
| --- | --- |
| `code/` | Benchmark workflow snapshots and EDTA scoring code |
| `benchmark/` | Discovery and reference preparation, base- and family-level scoring, discovery ablations, evidence-gate tests, classifier evaluation and scheduler scripts |
| `benchmark/maize_mask/` and `benchmark/tetrimmer_chunks/` | Maize splitting, curation, merging, masking and parallel scoring |
| `figures/scripts/` | Manuscript figure and table scripts with the small aggregate tables used for rendering |
| `supp/reference_libraries.tsv` | Assemblies are listed in the manuscript. This manifest identifies the reference consensus libraries, releases and scoring checksums. |

The v10 and v13 snapshots contain the exact workflow source collected from the benchmark checkouts. The v13 TEtrimmer module includes the patch for finishing from merged chunk results. The v9 and v9e snapshots support the original discovery runs. The snapshots are preserved so the reported results can be traced to the code that produced them.

## Review the results and regenerate figures

The companion reviewer backup contains the source archives, consensus FASTAs, exact genome assemblies, score tables, classifier datasets and maize masks. Extract the archives to a common directory. See its README and SHA-256 manifest. The five RepBase 29.10 reference libraries require licensed access. Follow the reference manifest and extraction script to obtain the scored versions.

Install the figure dependencies in a Python environment, then render the manuscript figures from the included summary tables.

```bash
python -m pip install -r requirements-figures.txt
python figures/scripts/fig1_pipeline.py
python figures/scripts/fig2_accuracy.py
python figures/scripts/fig3_classifiers.py
python figures/scripts/fig4_compute.py
```

To rebuild the summary tables from the archive, pass its extracted `benchmark/` directory to `python figures/scripts/extract_data.py <benchmark-root>`. The supplementary table and figure generator also reads raw benchmark files. See the source header of `figures/scripts/supplement.py` for its input layout. The included summary tables are the ones used for the manuscript.

## Repeat the benchmark runs

Use Linux, Nextflow and a compatible Java runtime. The workflow snapshots define the tool versions and environments. The discovery, masking and classifier scripts preserve the original Zhores paths and Slurm settings for provenance. Review the data roots, container and Conda paths, job accounts and resource settings before submitting jobs on another cluster. Start with `benchmark/heldout/run_heldout_v10.sh`, `benchmark/heldout/run_heldout_v13_maize.sh` and the maize controller and scoring scripts. Per-order base and family scoring are under `benchmark/edta_manual_bench/`. The classifier pool, split and evaluation scripts are under `benchmark/ensemble_v2/scripts/` and `benchmark/ensemble_v3_dfam37/`.

The eventual release should link the versioned benchmark code to the deposited reviewer data. See the archived versions' license files for supplied third-party notices.
