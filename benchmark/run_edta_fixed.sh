#!/bin/bash
# run_edta_fixed.sh <species> - EDTA 2.3.0 exactly as the pipeline's EDTA module (88bb57e):
# TIR-Learner patched at run time, fails on rice fallback. Output:
#   work/edta_fixed/<sp>/<sp>.edta_library.fa  (+ raw_library.fa for plant-route genomes)
set -eo pipefail
SP="$1"
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
source $B/env.sh
IMG=$WORK_DIR/singularity_cache/quay.io-biocontainers-edta-2.3.0--hdfd78af_0.img
KING=$(grep -vE '^[[:space:]]*(#|$)' $CONFIG_DIR/genomes.tsv | awk -F'\t' -v s="$SP" '$2==s{print $1}')
G=$(readlink -f $GENOME_DIR/$KING/$SP/$SP.fna)
D=$WORK_DIR/edta_fixed/$SP; mkdir -p $D; cd $D
ln -sf $G $SP.fna
CPUS=${SLURM_CPUS_PER_TASK:-16}
cat > edta_task.sh <<TASK
set -eo pipefail
cd $D
export HOME="\$PWD" MPLCONFIGDIR="\$PWD/.mplconfig" XDG_CACHE_HOME="\$PWD/.cache"
mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"
TL_WRAPPER=\$(command -v TIR-Learner || true)
TL_SRC=\$(grep -oE '/[^ ]*TIR-Learner[^ /]*\.py' "\${TL_WRAPPER:-/dev/null}" 2>/dev/null | head -1 || true)
rm -rf tirlearner; mkdir -p tirlearner/bin
cp -r "\$(dirname "\${TL_SRC}")" tirlearner/TIR-Learner
sed -i 's/^\(\s*\)family = x\[0\]\$/\1family = x.iloc[0]/' tirlearner/TIR-Learner/bin/check_TIR_TSD.py
sed -i 's/id_pattern = r"^(\\\\w+)_split_/id_pattern = r"^(\\\\S+?)_split_/' tirlearner/TIR-Learner/bin/run_GRF.py
printf '#!/bin/bash\npython3 %s "\$@"\n' "\$PWD/tirlearner/TIR-Learner/\$(basename "\${TL_SRC}")" > tirlearner/bin/TIR-Learner
chmod +x tirlearner/bin/TIR-Learner
export PATH="\$PWD/tirlearner/bin:\$PATH"
grep -c 'x.iloc\[0\]' tirlearner/TIR-Learner/bin/check_TIR_TSD.py; grep '^ *id_pattern' tirlearner/TIR-Learner/bin/run_GRF.py
EDTA.pl --genome $SP.fna --species others --step all --anno 0 --force 1 --threads $CPUS --overwrite 0
for t in LTR:LTR.raw.fa LINE:LINE.raw.fa SINE:SINE.raw.fa TIR:TIR.intact.raw.fa Helitron:Helitron.intact.raw.fa; do
    RAW=\$(find . -maxdepth 2 -path '*.EDTA.raw/*' -name "*.mod.\${t#*:}" -print -quit)
    RICE=/usr/local/share/EDTA/database/rice7.0.0.liban.\${t%%:*}
    if [[ -n "\$RAW" ]] && cmp -s "\$RAW" "\$RICE"; then echo "RICE FALLBACK: \${t%%:*}"; exit 3; fi
done
cp \$(find . -maxdepth 1 -name '*.mod.EDTA.TElib.fa' -print -quit) $SP.edta_library.fa
TASK
echo "[$(date)] $SP EDTA on $(hostname), $CPUS cpus"
singularity exec --no-home --bind $B --bind $(readlink -f $B) $IMG bash $D/edta_task.sh
echo "[$(date)] EDTA library: $(grep -c '>' $SP.edta_library.fa) families; raw TIR $(grep -c '>' $SP.fna.mod.EDTA.raw/$SP.fna.mod.TIR.intact.raw.fa), Helitron $(grep -c '>' $SP.fna.mod.EDTA.raw/$SP.fna.mod.Helitron.intact.raw.fa)"
# plant route: rebuild the raw de novo library with the fixed EDTA, in the pipeline's order
O=$WORK_DIR/pipeline/$SP/out/$SP/02_denovo_te
if [[ -s $O/$SP.raw_library.fa && -d $O/edta ]]; then
    ORDER=$(for f in $O/*/$SP.*_library.fa; do
        h=$(grep -m1 '>' $f); n=$(grep -nFm1 "$h" $O/$SP.raw_library.fa | cut -d: -f1); echo "$n $f"; done | sort -n | awk '{print $2}')
    : > raw_library.fa
    for f in $ORDER; do
        if [[ $f == */edta/* ]]; then cat $SP.edta_library.fa >> raw_library.fa; else cat $f >> raw_library.fa; fi
        echo "  raw part: $f"
    done
    echo "[$(date)] raw library: $(grep -c '>' raw_library.fa) seqs (old: $(grep -c '>' $O/$SP.raw_library.fa))"
fi
