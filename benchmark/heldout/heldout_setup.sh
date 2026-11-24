#!/bin/bash
# heldout_setup.sh - genomes.tsv rows, image cache union, Blumeria/Ectocarpus curated libs,
# EarlGrey standalone launcher; then submit the stub wiring test.
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
W=$B/work; H=$B/heldout
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
cd $B

# 1. manifest rows (appended; addressed by key, never by position here)
cp -n config/genomes.tsv config/genomes.tsv.pre_heldout
add() { grep -qP "^[^\t]+\t$2\t" config/genomes.tsv && { echo "row exists: $2"; return; }; printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> config/genomes.tsv; echo "added: $2"; }
add animal       caenorhabditis_elegans    "Caenorhabditis elegans (WBcel235)"          GCF_000002985.6 animal gold repbase_6239
add animal       caenorhabditis_briggsae   "Caenorhabditis briggsae (CB4)"              GCF_000004555.2 animal gold repbase_6238
add fungus       blumeria_hordei           "Blumeria hordei DH14 (BGH_DH14_v4)"         GCA_900239735.1 fungus gold qian2023
add plant        chlamydomonas_reinhardtii "Chlamydomonas reinhardtii (v5.5)"           GCF_000002595.2 plant  gold craig
add stramenopile ectocarpus_sp7            "Ectocarpus sp. 7 Ec32 (instaGRAAL)"         GCA_004764655.1 animal gold dinatale2025
add animal       gallus_gallus             "Gallus gallus (GRCg6a)"                     GCF_000002315.5 animal gold repbase_9031
add animal       taeniopygia_guttata       "Taeniopygia guttata (taeGut3.2.4)"          GCF_000151805.1 animal gold repbase_59729
add animal       mus_musculus              "Mus musculus (GRCm38)"                      GCA_000001635.2 animal gold repbase_10090
tail -n 9 config/genomes.tsv

# 2. image cache for the held-out runs: union of the per-run caches (private classifier
#    images) + the shared cache, as symlinks, so no task has to pull
mkdir -p $W/singularity_heldout
for d in $W/pipeline/drosophila_melanogaster/work/singularity $W/pipeline/arabidopsis_thaliana/work/singularity $W/singularity_cache; do
    for f in $d/*.img; do [[ -e $W/singularity_heldout/$(basename $f) ]] || ln -s $(readlink -f $f) $W/singularity_heldout/; done
done
ls $W/singularity_heldout

# 3. curated libraries with publisher classification (labels normalised only)
mkdir -p curated/blumeria_hordei curated/ectocarpus_sp7
$PY $H/normalize_lib.py $H/src/B_hordei_TE_manually_curated_consensus_library.fasta curated/blumeria_hordei/blumeria_hordei.curated.lib.fa
echo "source: github.com/qjiangzhao/TransposableElement-identification-alpha B_hordei_TE_manually_curated_consensus_library.fasta (Qian et al. 2023 Mob DNA 14:17; EDTA-seeded manual curation, TEtrimmer group); labels normalised by heldout/normalize_lib.py" > curated/blumeria_hordei/SOURCE.txt
$PY $H/normalize_lib.py $H/src/ectocarpus_MOESM1_ESM.fa curated/ectocarpus_sp7/ectocarpus_sp7.curated.lib.fa
echo "source: Dinatale et al. 2025 Genome Biol 26:320, Additional file 1 (13059_2025_3742_MOESM1_ESM.fa, Springer static content); curated on Ec32 V5, applied to GCA_004764655.1 (same strain; V5 not public); labels normalised by heldout/normalize_lib.py" > curated/ectocarpus_sp7/SOURCE.txt

# 4. EarlGrey standalone launcher (plant-route genomes have no EarlGrey in the pipeline)
cat > $H/run_earlgrey_standalone.sh <<'EOS'
#!/bin/bash
# run_earlgrey_standalone.sh <species> - EarlGrey baseline (04_run_tools.sh, mirrors the
# pipeline's EarlGrey module) with the pipeline's cached EarlGrey conda env.
set -eo pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
. /trinity/shared/opt/anaconda3/etc/profile.d/conda.sh
export EARLGREY_ENV=$B/work/conda_cache/env-e8cd3caef48f718c-8aca2d793016020c698f928e1801d622
export THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "[$(date)] EarlGrey standalone $1 on $(hostname), $THREADS threads"
bash $B/04_run_tools.sh "$1" earlgrey
L=$B/libraries/$1/earlgrey/earlgrey_library.fa
[[ -s $L ]] || { echo "no EarlGrey library for $1" >&2; exit 1; }
echo "[$(date)] EarlGrey library: $(grep -c '>' $L) seqs"
EOS
chmod +x $H/*.sh
echo "setup done"

