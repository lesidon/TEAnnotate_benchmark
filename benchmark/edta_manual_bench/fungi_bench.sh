#!/bin/bash
# After the fungal pipeline runs: snapshot the 5 libraries per species into specs.tsv,
# RepeatMasker them (rm_one.sh, EDTA-manual options), then lib-test + CSV.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
cd $E
FIRST=$(( $(grep -c . specs.tsv) ))      # 0-based index of the first new row
DONE=()
for SP in septoria_linicola zymoseptoria_tritici; do
    O=$B/work/pipeline/$SP/out/$SP/02_denovo_te
    declare -A L=(
        [reference]=$B/curated/$SP/$SP.curated.lib.fa
        [teannotate]=$B/libraries/$SP/pipeline/pipeline_library.fa
        [hite]=$O/hite/$SP.hite_library.fa
        [edta]=$B/libraries/$SP/edta/edta_library.fa
        [earlgrey]=$O/earlgrey/$SP.earlgrey_library.fa )
    ok=1
    for m in reference teannotate hite edta earlgrey; do
        [[ -s ${L[$m]} ]] || { echo "SKIP $SP: missing $m library ${L[$m]}"; ok=0; }
    done
    # the pipeline library must come from THIS run, not an older one
    [[ ${L[teannotate]} -nt $B/run_fungus.sh ]] || { echo "SKIP $SP: pipeline library is stale"; ok=0; }
    (( ok )) || continue
    grep -q "^$SP	" specs.tsv && { echo "$SP already in specs.tsv"; continue; }
    for m in reference teannotate hite edta earlgrey; do
        mkdir -p $SP/$m; cp "${L[$m]}" $SP/$m/library.fa
        printf '%s\t%s\t%s\n' "$SP" "$m" "${L[$m]}" >> specs.tsv
        echo "$SP/$m: $(grep -c '>' $SP/$m/library.fa) seqs"
    done
    DONE+=($SP)
done
LAST=$(( $(grep -c . specs.tsv) - 1 ))
(( ${#DONE[@]} )) || { echo "nothing to benchmark"; exit 1; }
A=$(sbatch --parsable -p htc --array=$FIRST-$LAST -c 16 --mem=32G -t 12:00:00 --job-name=edta_rm_fungi \
      -o "$E/logs/rm_%a.out" -e "$E/logs/rm_%a.err" "$E/rm_one.sh")
L=$(sbatch --parsable -p htc -c 4 --mem=16G -t 04:00:00 --job-name=edta_libtest_fungi --dependency=afterany:$A \
      -o "$E/logs/libtest_fungi.out" -e "$E/logs/libtest_fungi.err" "$E/libtest_species.sh" "${DONE[@]}")
echo "RM array $A (rows $FIRST-$LAST), lib-test $L for ${DONE[*]}"
