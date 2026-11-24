#!/bin/bash
# Code-path check of tune_test.py on the July PanTEon fungi voting table (random split).
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
X=$B/ensemble_v2; P=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/ensemble_prediction
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
S=$X/selftest; rm -rf $S; mkdir -p $S/pools $S/split $S/preds/fungi
cp $(find $P/predictions/panteon_fungi_final -name '*voting_table.tsv' | head -1) $S/preds/fungi/voting_table.tsv
$PY - <<PY
import csv, random, sys
sys.path.insert(0, "$X/scripts")
from ens_labels import order_sf
random.seed(1)
rows = list(csv.DictReader(open("$P/panteon_partitions/panteon_fungi.tsv"), delimiter="\t"))
with open("$S/pools/fungi.truth.tsv", "w") as t, open("$S/split/fungi.split.tsv", "w") as s:
    t.write("family_id\tspecies\tkingdom\tsource\tlabel_raw\torder\tsuperfamily\theldout_species\n")
    s.write("family_id\tcluster\tsplit\torder\tsource\theldout_species\n")
    for r in rows:
        raw = r["label"].split("#", 1)[-1]; o, sf = order_sf(raw)
        if o in ("Unknown", "Other"): continue
        t.write("\t".join([r["family_id"], r["species"], "fungi", "PanTEon", raw, o, sf, "0"]) + "\n")
        s.write("\t".join([r["family_id"], r["family_id"], random.choice(["train", "test"]), o, "PanTEon", "0"]) + "\n")
PY
cd $S
time $PY $X/scripts/tune_test.py --kingdom fungi --workdir $S --voter $B/pipeline_head_v9/bin \
    --label-map $X/scripts/label_map.json --threads 8 --n-boot 50 2>&1 | tail -25
head -5 $S/tune/fungi.top_configs.tsv

