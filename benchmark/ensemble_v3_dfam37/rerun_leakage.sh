#!/bin/bash
#SBATCH --job-name=leak_dfam37
#SBATCH --partition=mem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360G
#SBATCH --time=3-00:00:00
#SBATCH --output=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v3_dfam37/logs/rerun_%j.out
# Classifier benchmark with the leakage filter extended to Dfam 3.7 (curated + uncurated families),
# the training data of TEclass2 that the first filter (RepBase 29.10 + curated Dfam 3.8) missed.
#  1. Dfam.embl.gz (md5-checked) -> FASTA
#  2. leakage_filter.py of each kingdom pool of ensemble_v2 against it (same 0.7 identity / 0.5 coverage)
#  3. new workdir = ensemble_v2 truth + split minus the leaked consensi (both halves; same clusters)
#  4. tune_test.py + two_step.py (support-weighted objective, .wf1) on the filtered halves
#  5. paired cluster bootstrap: re-tuned settings, and the deployed settings (ensemble_v2/tune)
# ensemble_v2 itself is only read.
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
X=$B/ensemble_v2
W=$B/ensemble_v3_dfam37
V=$B/pipeline_head_v9/bin
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
MMSEQS=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/mmseqs
T=${SLURM_CPUS_PER_TASK:-32}
cd $W
echo "[$(date)] start on $(hostname), $T CPUs"

# 1. download check + conversion
(cd data/dfam37 && md5sum -c Dfam.embl.gz.md5sum)
if [ ! -s data/dfam37/dfam37.fasta ]; then
    $PY $X/scripts/prepare_training_fasta.py --source dfam37=data/dfam37/Dfam.embl.gz --outdir data/dfam37
fi
echo "[$(date)] Dfam 3.7 FASTA: $(grep -c '^>' data/dfam37/dfam37.fasta) families" \
     "($(grep -c '^>dfam37__DF' data/dfam37/dfam37.fasta) DF curated, $(grep -c '^>dfam37__DR' data/dfam37/dfam37.fasta) DR uncurated)"

# 2. leakage filter per kingdom (resumable: a finished kingdom is skipped)
mkdir -p leakage/tmp
for k in plant animal fungi; do
    [ -s leakage/$k.leakage_summary.tsv ] && { echo "$k: done before"; continue; }
    echo "[$(date)] leakage filter: $k"
    $PY $X/scripts/leakage_filter.py --query $X/pools/$k.fasta --training dfam37=data/dfam37/dfam37.fasta \
        --identity 0.7 --coverage 0.5 --threads $T --mmseqs $MMSEQS \
        --tmpdir $W/leakage/tmp --outdir leakage --prefix $k
    rm -f leakage/$k.heldout.fasta        # a copy of the pool minus the leaked; the id lists suffice
    cat leakage/$k.leakage_summary.tsv
done

# 3. filtered workdir (same layout as ensemble_v2)
mkdir -p pools split preds tune test results
ln -sfn $X/scripts scripts
for k in plant animal fungi; do ln -sfn $X/preds/$k preds/$k; done
$PY - <<'PY'
import csv, collections
X = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/ensemble_v2"
out = open("results/filter_summary.tsv", "w")
out.write("kingdom\tsource\tsplit\tbefore\tremoved\tafter\tremoved_by_DR_uncurated\n")
for k in ("plant", "animal", "fungi"):
    rep = {r["family_id"]: r for r in csv.DictReader(open("leakage/%s.leakage_report.tsv" % k), delimiter="\t")}
    removed = {f for f, r in rep.items() if r["status"] == "removed"}
    truth = list(csv.DictReader(open("%s/pools/%s.truth.tsv" % (X, k)), delimiter="\t"))
    split = list(csv.DictReader(open("%s/split/%s.split.tsv" % (X, k)), delimiter="\t"))
    miss = sum(1 for r in truth if r["family_id"] not in rep)
    print("%s: %d truth ids, %d not in the leakage report, %d removed" % (k, len(truth), miss, len(removed)))
    if miss > 0.01 * len(truth):
        raise SystemExit("ERROR: truth ids do not match the pool FASTA ids")
    for name, rows in (("pools/%s.truth.tsv" % k, truth), ("split/%s.split.tsv" % k, split)):
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rows:
                if r["family_id"] not in removed:
                    w.writerow(r)
    c = collections.Counter((r["source"], r["split"]) for r in split)
    d = collections.Counter((r["source"], r["split"]) for r in split if r["family_id"] in removed)
    dr = collections.Counter((r["source"], r["split"]) for r in split if r["family_id"] in removed
                             and rep[r["family_id"]]["best_target"].startswith("dfam37__DR"))
    for key in sorted(c):
        out.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (k, key[0], key[1], c[key], d[key], c[key] - d[key], dr[key]))
out.close()
print(open("results/filter_summary.tsv").read())
PY

# 4. tuning on the filtered training halves, one evaluation on the filtered test halves
for k in plant animal fungi; do
    echo "[$(date)] tune_test + two_step: $k"
    $PY scripts/tune_test.py --kingdom $k --workdir $W --voter $V --label-map scripts/label_map.json \
        --objective weighted --tag .wf1 --threads $T > logs/tune_test_$k.log 2>&1
    $PY scripts/two_step.py --kingdom $k --workdir $W --voter $V --label-map scripts/label_map.json \
        --objective weighted --tag .wf1 --threads $T > logs/two_step_$k.log 2>&1
done

# 5. paired cluster bootstrap: re-tuned settings, then the deployed (ensemble_v2) settings
$PY $W/paired_boot.py $W $V 2000 results/paired_boot_retuned.tsv $W/tune
$PY $W/paired_boot.py $W $V 2000 results/paired_boot_deployed.tsv $X/tune
for k in plant animal fungi; do
    echo "== $k: deployed vs re-tuned settings"; cat $X/tune/$k.wf1.two_step.best.json | tr -d ' \n'; echo
    cat tune/$k.wf1.two_step.best.json | tr -d ' \n'; echo
done
rm -rf leakage/tmp
echo "[$(date)] done"
