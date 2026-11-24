#!/bin/bash
# refined_bench_v10.sh - final development rerun (pipeline_head_v10) into the benchmark, as
# refined_bench_v9.sh: specs rows -> RepeatMasker -> lib-test (default + harmonised) ->
# family level -> summary of v10 next to v9 and the standalone tools.
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
SPECIES="zymoseptoria_tritici septoria_linicola drosophila_melanogaster arabidopsis_thaliana"
NEW="teannotate_ref10 teannotate_wo_ref10"
COMPARE="teannotate_ref9 teannotate_wo_ref9 teannotate_ref7 teannotate_wo_ref7 teannotate_ref4 teannotate_wo_ref4 hite earlgrey edta2 heliano"
cd $E
ROWS=(); DONE=()
for SP in $SPECIES; do
    ok=0
    for arm in teannotate_ref10:pipeline_refined_v10 teannotate_wo_ref10:pipeline_wo_tetrimmer_refined_v10; do
        M=${arm%%:*}; L=$B/libraries/$SP/${arm#*:}/pipeline_library.fa
        [[ -s $L ]] || { echo "SKIP $SP/$M: no library"; continue; }
        ok=1
        mkdir -p $SP/$M; cp "$L" $SP/$M/library.fa
        n=$(grep -nP "^$SP\t$M\t" specs.tsv | cut -d: -f1)
        if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$SP" "$M" "$L" >> specs.tsv; n=$(grep -c . specs.tsv); fi
        ROWS+=($(( n - 1 )))
        echo "$SP/$M: $(grep -c '>' $SP/$M/library.fa) seqs (specs row $(( n - 1 )))"
    done
    (( ok )) && DONE+=($SP)
done
(( ${#DONE[@]} )) || { echo "nothing to benchmark"; exit 1; }
echo "[$(date)] RepeatMasker rows: ${ROWS[*]}"
for r in "${ROWS[@]}"; do SLURM_ARRAY_TASK_ID=$r bash $E/rm_one.sh > $E/logs/rm_${r}_v10.out 2>&1 & done
wait
for r in "${ROWS[@]}"; do tail -n 1 $E/logs/rm_${r}_v10.out; done
echo "[$(date)] lib-test (default, new rows) + harmonised"
for SP in "${DONE[@]}"; do
    KING=$(grep -vE '^[[:space:]]*(#|$)' $B/config/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
    G=$B/genomes/$KING/$SP/$SP.fna
    for M in $NEW; do
        [[ -s $E/$SP/$M/$SP.fna.out ]] || continue
        ( cd $E/$SP/$M && $PY $B/tools/EDTA/lib-test.py --genome $G --reference $E/$SP/reference/$SP.fna.out \
              --test $E/$SP/$M/$SP.fna.out --extended_report ) > $E/$SP/$M/libtest.log 2>&1 &
    done
    $PY $E/libtest_harmonised.py $SP $NEW $COMPARE > $E/logs/harm_v10_$SP.log 2>&1 &
done
wait
echo "[$(date)] family level"
SLURM_CPUS_PER_TASK=16 bash $E/family_level.sh "${DONE[@]}" > $E/logs/famlevel_v10.log 2>&1
echo "[$(date)] summary"
$PY - "$NEW teannotate_ref9 teannotate_wo_ref9 hite earlgrey edta2 heliano" "${DONE[@]}" <<'PY'
import csv, os, re, sys
E = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/edta_manual_bench"
LABEL = {"teannotate_ref9": "teannotate with tetrimmer (v9)", "teannotate_wo_ref9": "teannotate w/o tetrimmer (v9)",
         "edta2": "edta"}
methods = sys.argv[1].split()
for sp in sys.argv[2:]:
    nice = sp[0].upper() + sp[1:]
    print(f"\n{nice}")
    harm = list(csv.DictReader(open(f"{E}/results/{nice}_bench_harmonised.csv")))
    fam = list(csv.DictReader(open(f"{E}/results/{nice}_family_level.csv")))
    for m in methods:
        h = {x["TE type"]: x for x in harm if x["Method"] == LABEL.get(m, m)}
        rep = f"{E}/{sp}/{m}/{sp}.fna.out.repeatmasker.lib.report"
        d = {}
        if os.path.exists(rep):
            for line in open(rep):
                f = line.rstrip("\n").split("\t")
                mm = re.search(r"\.(\w+)\.lib\.report$", f[0])
                if mm and len(f) >= 11:
                    d[mm.group(1)] = f
        if not h:
            continue
        def g(c):
            x = h.get(c)
            return f"{float(x['F1']):.3f} (s{float(x['sens']):.2f}/p{float(x['prec']):.2f})" if x else "  -  "
        fa = [x for x in fam if x["Method"] == m and x["class"] == "all"]
        fs = f"fam {fa[0]['perfect']}/{fa[0]['good']}/{fa[0]['present_only']} ({fa[0]['n_library']})" if fa else ""
        dany = f"{float(d['any'][6]):.3f}" if 'any' in d else "-"
        print(f"  {m:20s} harm any {g('any')}  any_dc {g('any_dontcare')}  ltr {g('ltr')}  nonltr {g('nonltr')}  "
              f"tir {g('tir')}  helitron {g('helitron')} | default any F1 {dany} | {fs}")
PY
echo "[$(date)] all done"

