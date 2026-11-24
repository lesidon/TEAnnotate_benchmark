#!/bin/bash
set -euo pipefail
R=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/maize_mask/parallel_scoring_20261120
cd "$R"
test ! -e job_receipts.tsv
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
# Refuse to schedule a supposedly ready method with missing chunks or masks.
"$PY" - <<'PY'
import parallel_scoring as p
for method in p.METHODS[:5]:
    mod = p.scope(p.original(), p.ROOT, method)
    mod.validate_masks()
print('Five ready methods validated')
PY
record() { printf '%s\t%s\n' "$1" "$2" | tee -a job_receipts.tsv; }
INDEX=$(sbatch --hold --parsable -J mz_score_index -p htc -c 1 --mem=16G -t 1-00:00:00 --exclude='ct[01-10],gn[19,26]' -o "$R/index_%j.out" score.sbatch index)
record index "$INDEX"
FAMILY=$(sbatch --hold --parsable --array=0-4 -J mz_family_ready -p htc -c 16 --mem=64G -t 1-00:00:00 --exclude='ct[01-10],gn[19,26]' -o "$R/family_%A_%a.out" score.sbatch family)
record family_ready "$FAMILY"
BASE=$(sbatch --hold --parsable --array=0-4 -J mz_base_ready -p mem -c 4 --mem=120G -t 3-00:00:00 --exclude='ct[01-04]' --dependency=afterok:$INDEX -o "$R/base_%A_%a.out" score.sbatch base)
record base_ready "$BASE"
LATE_FAMILY=$(sbatch --hold --parsable --array=5 -J mz_family_tetrimmer -p htc -c 16 --mem=64G -t 1-00:00:00 --exclude='ct[01-10],gn[19,26]' --dependency=afterok:8559231 -o "$R/family_%A_%a.out" score.sbatch family)
record family_tetrimmer "$LATE_FAMILY"
LATE_BASE=$(sbatch --hold --parsable --array=5 -J mz_base_tetrimmer -p mem -c 4 --mem=120G -t 3-00:00:00 --exclude='ct[01-04]' --dependency=afterok:8559231:$INDEX -o "$R/base_%A_%a.out" score.sbatch base)
record base_tetrimmer "$LATE_BASE"
MERGE=$(sbatch --hold --parsable -J mz_scores_merge -p mem -c 1 --mem=8G -t 03:00:00 --exclude=ct01 --dependency=afterok:$FAMILY:$BASE:$LATE_FAMILY:$LATE_BASE -o "$R/merge_%j.out" score.sbatch merge)
record merge "$MERGE"
# Replace the obsolete serial branches only once every replacement job exists.
scontrol hold 8557098 8557099
scontrol update JobId=8557100 Dependency=afterok:$MERGE
scancel 8557098 8557099
for id in "$INDEX" "$FAMILY" "$BASE" "$LATE_FAMILY" "$LATE_BASE" "$MERGE"; do
    scontrol release "$id"
done
date
squeue -r -j "$INDEX,$FAMILY,$BASE,$LATE_FAMILY,$LATE_BASE,$MERGE,8557100,8559231,8562907" -o '%.18i %.24j %.12T %.100E %.30R' | tee jobs_after.txt
