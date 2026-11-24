#!/usr/bin/env python3
"""Aggregate gold (lib-test) + proxy metrics into one tidy master table.

Reads:
    results/gold/<species>.<method>.metrics.tsv   (from parse_libtest.py)
    results/proxy/<species>.coverage.tsv          (from proxy_metrics.py)
    results/proxy/<species>.concordance.tsv
    results/proxy/<species>.lai.tsv

Writes (into --outdir):
    benchmark_metrics.tsv            long: kingdom species method scoring_type category metric_name value
    benchmark_gold_F1.tsv            wide: gold F1 by species x method x category
    benchmark_proxy_coverage.tsv     wide: proxy coverage_pct by species x method
    concordance_all.tsv              all pairwise Jaccard rows, kingdom-annotated
"""
import argparse
import glob
import os

import pandas as pd

GOLD_METRICS = ["sensitivity", "specificity", "accuracy", "precision", "FDR", "F1"]


def species_kingdom_map(genomes_tsv):
    m = {}
    with open(genomes_tsv, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) >= 2:
                m[f[1]] = f[0]
    return m


def load_gold(gold_dir):
    rows = []
    for path in sorted(glob.glob(os.path.join(gold_dir, "*.metrics.tsv"))):
        df = pd.read_csv(path, sep="\t", dtype=str)
        if df.empty:
            continue
        long = df.melt(
            id_vars=["kingdom", "species", "method", "scoring_type", "category"],
            value_vars=[c for c in GOLD_METRICS if c in df.columns],
            var_name="metric_name", value_name="value")
        rows.append(long)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def load_proxy(proxy_dir, k_of):
    long_rows, concord = [], []

    for path in sorted(glob.glob(os.path.join(proxy_dir, "*.coverage.tsv"))):
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            for metric in ("coverage_pct", "unique_bp", "masked_bp"):
                long_rows.append(dict(kingdom=r["kingdom"], species=r["species"],
                                      method=r["method"], scoring_type="proxy",
                                      category="genome", metric_name=metric,
                                      value=r[metric]))

    for path in sorted(glob.glob(os.path.join(proxy_dir, "*.concordance.tsv"))):
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        sp = df["species"].iloc[0]
        df["kingdom"] = k_of.get(sp, "NA")
        concord.append(df)
        # per-method mean pairwise Jaccard (exclude self) -> a summary metric.
        off = df[df["method_a"] != df["method_b"]]
        jac = {}
        for _, r in off.iterrows():
            jac.setdefault(r["method_a"], []).append(r["jaccard"])
            jac.setdefault(r["method_b"], []).append(r["jaccard"])
        for method, vals in jac.items():
            long_rows.append(dict(kingdom=k_of.get(sp, "NA"), species=sp, method=method,
                                  scoring_type="proxy", category="genome",
                                  metric_name="mean_jaccard",
                                  value=round(sum(vals) / len(vals), 4)))

    for path in sorted(glob.glob(os.path.join(proxy_dir, "*.lai.tsv"))):
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            long_rows.append(dict(kingdom=k_of.get(r["species"], "NA"), species=r["species"],
                                  method="(genome)", scoring_type="proxy", category="genome",
                                  metric_name="LAI", value=r["lai"]))

    long = pd.DataFrame(long_rows)
    concord = pd.concat(concord, ignore_index=True) if concord else pd.DataFrame()
    return long, concord


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--results", default=os.path.join(here, "results"),
                    help="results dir containing gold/ and proxy/")
    ap.add_argument("--genomes", default=os.path.join(here, "config", "genomes.tsv"))
    ap.add_argument("--outdir", default=None, help="default: the results dir")
    args = ap.parse_args(argv)
    outdir = args.outdir or args.results
    os.makedirs(outdir, exist_ok=True)

    k_of = species_kingdom_map(args.genomes)
    gold = load_gold(os.path.join(args.results, "gold"))
    proxy, concord = load_proxy(os.path.join(args.results, "proxy"), k_of)

    master = pd.concat([d for d in (gold, proxy) if not d.empty], ignore_index=True)
    if master.empty:
        raise SystemExit("[aggregate] no metrics found under " + args.results)
    cols = ["kingdom", "species", "method", "scoring_type", "category", "metric_name", "value"]
    # Sort by the descriptive keys only — NOT 'value', which mixes gold strings and
    # proxy floats (comparing str vs float raises TypeError in Python 3).
    master = master[cols].sort_values(cols[:-1]).reset_index(drop=True)
    master_path = os.path.join(outdir, "benchmark_metrics.tsv")
    master.to_csv(master_path, sep="\t", index=False)

    # Convenience wides.
    if not gold.empty:
        f1 = gold[gold["metric_name"] == "F1"]
        (f1.pivot_table(index=["species", "method"], columns="category",
                        values="value", aggfunc="first")
           .to_csv(os.path.join(outdir, "benchmark_gold_F1.tsv"), sep="\t"))
    if not proxy.empty:
        covw = proxy[proxy["metric_name"] == "coverage_pct"]
        (covw.pivot_table(index="species", columns="method", values="value", aggfunc="first")
             .to_csv(os.path.join(outdir, "benchmark_proxy_coverage.tsv"), sep="\t"))
    if not concord.empty:
        concord.to_csv(os.path.join(outdir, "concordance_all.tsv"), sep="\t", index=False)

    print(f"[aggregate] {len(master)} rows -> {master_path}")
    print(f"[aggregate]   gold rows:  {len(gold)}")
    print(f"[aggregate]   proxy rows: {len(proxy)}")


if __name__ == "__main__":
    main()
