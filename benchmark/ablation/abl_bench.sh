#!/bin/bash
# abl_bench.sh [arm] - score the discovery-tool ablation like the dev bench: every
# combination's final library (combos/<sp>/<combo>/library_<arm>.fa, method abl_<arm>_<combo>)
# and every fresh single-tool raw library (libs/<sp>/<tool>.fa, method fresh_<tool>: the
# standalone baselines, EDTA included) -> rm_one.sh -> harmonised lib-test (+ any_dontcare)
# -> family level. Results: results/<Sp>_bench_harmonised_abl.csv (a copy of the scorer;
# the shared CSV is untouched) and ablation/results/ablation_<arm>.tsv with compute cost:
# discovery CPU-h (sum over the combination's tools of cpus x realtime), discovery wall-h
# (the longest tool: they run in parallel) and library-stage CPU-h (run_reclassify trace).
ARM="${1:-wo}"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
E=$B/edta_manual_bench; A=$B/ablation
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
source $B/env.sh
# genomes and extra comparison methods can be set from the environment (default: the
# first four development genomes, next to their v10 pipeline libraries)
SPECIES="${SPECIES:-septoria_linicola zymoseptoria_tritici drosophila_melanogaster arabidopsis_thaliana}"
COMPARE="${COMPARE-teannotate_ref10 teannotate_wo_ref10}"
sed 's/_bench_harmonised.csv/_bench_harmonised_abl.csv/' $E/libtest_harmonised.py > $A/libtest_harmonised_abl.py
mkdir -p $A/results $E/logs
cd $E
ROWS=()
add() {  # <sp> <method> <library>
    local SP=$1 M=$2 L=$3 n
    [[ -s $L ]] || return 0
    mkdir -p $SP/$M
    if cmp -s $L $SP/$M/library.fa && [[ -s $SP/$M/$SP.fna.out ]]; then return 0; fi
    cp $L $SP/$M/library.fa
    n=$(grep -nF -- "$SP	$M	" specs.tsv | head -1 | cut -d: -f1)
    if [[ -z $n ]]; then printf '%s\t%s\t%s\n' "$SP" "$M" "$L" >> specs.tsv; n=$(grep -c . specs.tsv); fi
    ROWS+=($(( n - 1 )))
}
for SP in $SPECIES; do
    for t in earlgrey edta heliano hite; do add $SP fresh_$t $A/libs/$SP/$t.fa; done
    for D in $A/combos/$SP/*/; do CB=$(basename $D); add $SP abl_${ARM}_$CB $D/library_$ARM.fa; done
done
echo "[$(date)] RepeatMasker rows (${#ROWS[@]}): ${ROWS[*]}"
for r in "${ROWS[@]}"; do
    SLURM_ARRAY_TASK_ID=$r bash $E/rm_one.sh > $E/logs/rm_${r}_abl.out 2>&1 &
    while (( $(jobs -rp | wc -l) >= 4 )); do sleep 20; done
done
wait
echo "[$(date)] harmonised lib-test + family level"
for SP in $SPECIES; do
    M=$(ls -d $E/$SP/fresh_* $E/$SP/abl_${ARM}_* 2>/dev/null | xargs -n1 basename | tr '\n' ' ')
    $PY $A/libtest_harmonised_abl.py $SP $COMPARE $M > $A/results/harm_${SP}_$ARM.log 2>&1 &
done
wait
SLURM_CPUS_PER_TASK=16 bash $E/family_level.sh $SPECIES > $A/results/famlevel_$ARM.log 2>&1
echo "[$(date)] summary"
$PY - $ARM $SPECIES <<'PY'
import csv, glob, os, sys
arm, species = sys.argv[1], sys.argv[2:]
B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E, A = f"{B}/edta_manual_bench", f"{B}/ablation"

def cpu_h(trace, names=None):
    """CPU-hours and wall-hours (longest task) of a raw Nextflow trace."""
    cpu = wall = 0.0
    if not os.path.exists(trace):
        return None, None
    for r in csv.DictReader(open(trace), delimiter='\t'):
        if r.get('status') not in ('COMPLETED', 'CACHED') or r.get('realtime') in (None, '', '-'):
            continue
        h = float(r['realtime']) / 3.6e6
        cpu += h * float(r.get('cpus') or 1)
        wall = max(wall, h)
    return cpu, wall

rows = []
for sp in species:
    nice = sp[0].upper() + sp[1:]
    p = f"{E}/results/{nice}_bench_harmonised_abl.csv"
    if not os.path.exists(p):
        continue
    harm = list(csv.DictReader(open(p)))
    fam = list(csv.DictReader(open(f"{E}/results/{nice}_family_level.csv")))
    tool = {t: cpu_h(f"{A}/libs/{sp}/{t}.trace.tsv") for t in ('earlgrey', 'edta', 'heliano', 'hite')}
    for m in sorted({x['Method'] for x in harm}):
        h = {x['TE type']: float(x['F1']) for x in harm if x['Method'] == m}
        fa = [x for x in fam if x['Method'] == m and x['class'] == 'all']
        r = dict(species=sp, method=m, any=h.get('any'), any_dc=h.get('any_dontcare'), ltr=h.get('ltr'),
                 nonltr=h.get('nonltr'), tir=h.get('tir'), helitron=h.get('helitron'),
                 perfect=int(fa[0]['perfect']) if fa else None, good=int(fa[0]['good']) if fa else None,
                 n_lib=int(fa[0]['n_library']) if fa else None)
        tools = []
        if m.startswith(f'abl_{arm}_'):
            tools = m[len(f'abl_{arm}_'):].split('+')
            r['lib_cpu_h'], _ = cpu_h(f"{A}/combos/{sp}/{'+'.join(tools)}/trace_{arm}.tsv")
        elif m.startswith('fresh_'):
            tools = [m[6:]]
        if tools and all(tool[t][0] is not None for t in tools):
            r['disc_cpu_h'] = sum(tool[t][0] for t in tools)
            r['disc_wall_h'] = max(tool[t][1] for t in tools)
        rows.append(r)

# Selection rule (fixed before looking at the results): among the combinations of a genome,
# 'ok' = harmonised any_dontcare F1 within 0.01 (run-to-run noise) of the best combination and
# perfect families >= 95% of the best; 'pick' = the ok combination with the fewest discovery
# CPU-hours (ties: fewer library-stage CPU-hours).
for sp in species:
    cs = [r for r in rows if r['species'] == sp and r['method'].startswith(f'abl_{arm}_')
          and r['any_dc'] is not None and r['perfect'] is not None]
    if not cs:
        continue
    bf, bp = max(r['any_dc'] for r in cs), max(r['perfect'] for r in cs)
    for r in cs:
        r['ok'] = r['any_dc'] >= bf - 0.01 and r['perfect'] >= 0.95 * bp
    ok = [r for r in cs if r['ok'] and r.get('disc_cpu_h') is not None]
    if ok:
        min(ok, key=lambda r: (r['disc_cpu_h'], r.get('lib_cpu_h') or 0))['pick'] = True

cols = ['species', 'method', 'any', 'any_dc', 'ltr', 'nonltr', 'tir', 'helitron', 'perfect', 'good',
        'n_lib', 'disc_cpu_h', 'disc_wall_h', 'lib_cpu_h', 'ok', 'pick']
out = f"{A}/results/ablation_{arm}{os.environ.get('ABL_TAG', '')}.tsv"
with open(out, 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for r in rows:
        fh.write('\t'.join('' if r.get(c) is None else (f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r[c]))
                           for c in cols) + '\n')
f = lambda v, d=3: '-' if v is None else f"{v:.{d}f}"
for sp in species:
    print(f"\n{sp}")
    for r in sorted((r for r in rows if r['species'] == sp), key=lambda r: -(r['any_dc'] or 0)):
        flag = 'PICK' if r.get('pick') else ('ok' if r.get('ok') else '')
        print(f"  {r['method']:38s} any_dc {f(r['any_dc'])} tir {f(r['tir'])} hel {f(r['helitron'])}"
              f" | perf {r['perfect']} | disc {f(r.get('disc_cpu_h'),1)} CPU-h, {f(r.get('disc_wall_h'),1)} h wall"
              f" | lib {f(r.get('lib_cpu_h'),1)} CPU-h  {flag}")
print(f"\n-> {out}")
PY
echo "[$(date)] all done"
