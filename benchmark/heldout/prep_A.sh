#!/bin/bash
# heldout/prep_A.sh - download held-out genomes (NCBI) and curated libraries (GitHub / PMC
# supplements) into the benchmark tree. Re-entrant: finished items are skipped.
set -o pipefail
B=/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark
PY=/beegfs/home/leonid.sidorov/envs/teannot_bench/bin/python
H=$B/heldout; mkdir -p $H/dl $H/src
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"
log() { echo "[$(date '+%F %T')] $*"; }

dl_genome() {  # dl_genome <kingdom> <species_key> <accession>
    local k=$1 sp=$2 acc=$3
    local d=$B/genomes/$k/$sp z=$H/dl/$sp.zip t=$H/dl/$sp
    if [[ -s $d/$sp.fna && -s $d/provenance.json ]]; then log "skip genome $sp"; return 0; fi
    mkdir -p $d
    log "genome $sp ($acc) ..."
    for try in 1 2 3; do
        curl -sSL --retry 3 --max-time 3600 -o $z \
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" && break
        sleep 20
    done
    rm -rf $t; mkdir -p $t
    unzip -q -o $z -d $t || { log "ERROR unzip $sp"; return 1; }
    if [[ -f $t/md5sum.txt ]]; then
        ( cd $t && md5sum -c --quiet md5sum.txt ) || { log "ERROR md5 $sp"; return 1; }
    fi
    local fnas; fnas=$(find $t/ncbi_dataset/data/$acc -name '*.fna' | sort)
    [[ -n $fnas ]] || { log "ERROR no fna for $sp"; return 1; }
    cat $fnas > $d/$sp.fna
    local md5; md5=$(md5sum $d/$sp.fna | awk '{print $1}')
    local rep; rep=$(curl -s -m 60 "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/dataset_report")
    echo "$rep" | $PY -c "
import json,sys,datetime
r=json.load(sys.stdin)['reports'][0]
print(json.dumps({'species_key':'$sp','accession':'$acc','assembly_name':r['assembly_info'].get('assembly_name'),
 'organism':r['organism'].get('organism_name'),'total_length':r['assembly_stats'].get('total_sequence_length'),
 'fna_md5':'$md5','downloaded':datetime.date.today().isoformat(),'source':'NCBI Datasets API v2 (accession $acc)'},indent=1))
" > $d/provenance.json
    log "genome $sp: $(grep -c '^>' $d/$sp.fna) seqs, $(grep -v '^>' $d/$sp.fna | tr -d '\n' | wc -c) bp"
    rm -rf $t $z
}

dl_file() {  # dl_file <url> <dest>
    local url=$1 dest=$2
    if [[ -s $dest ]]; then log "skip $(basename $dest)"; return 0; fi
    curl -sSL --retry 3 --max-time 1800 -A "$UA" -o $dest "$url" || { log "ERROR download $url"; rm -f $dest; return 1; }
    log "got $(basename $dest): $(wc -c < $dest) bytes, type: $(file -b $dest | cut -c1-60)"
}

# --- genomes ---------------------------------------------------------------
dl_genome animal       caenorhabditis_elegans     GCF_000002985.6
dl_genome animal       caenorhabditis_briggsae    GCF_000004555.2
dl_genome fungus       blumeria_hordei            GCA_900239735.1
dl_genome plant        chlamydomonas_reinhardtii  GCF_000002595.2
dl_genome stramenopile ectocarpus_sp7             GCA_004764655.1
dl_genome animal       gallus_gallus              GCF_000002315.5
dl_genome animal       taeniopygia_guttata        GCF_000151805.1
dl_genome animal       mus_musculus               GCA_000001635.2

# --- curated library sources -------------------------------------------------
dl_file "https://raw.githubusercontent.com/qjiangzhao/TransposableElement-identification-alpha/main/B_hordei_TE_manually_curated_consensus_library.fasta" \
        $H/src/B_hordei_TE_manually_curated_consensus_library.fasta
dl_file "https://pmc.ncbi.nlm.nih.gov/articles/instance/12477816/bin/13059_2025_3742_MOESM1_ESM.fa" $H/src/ectocarpus_MOESM1_ESM.fa
dl_file "https://pmc.ncbi.nlm.nih.gov/articles/instance/12477816/bin/13059_2025_3742_MOESM3_ESM.xlsx" $H/src/ectocarpus_MOESM3_ESM.xlsx
dl_file "https://pmc.ncbi.nlm.nih.gov/articles/instance/8226300/bin/koab026_supplementary_data.zip" $H/src/craig2021_koab026_supplementary_data.zip
dl_file "https://pmc.ncbi.nlm.nih.gov/articles/instance/9940879/bin/koac347_supplementary_data.zip" $H/src/craig2023_koac347_supplementary_data.zip

log "### inspection"
for f in $H/src/*.zip; do echo "== $f"; unzip -l $f 2>&1 | head -60; done
for f in $H/src/*.fa $H/src/*.fasta; do
    echo "== $f: $(grep -c '^>' $f) seqs"; grep -m5 '^>' $f
    grep '^>' $f | sed 's/^>[^#]*#\{0,1\}//' | awk '{print $1}' | sort | uniq -c | sort -rn | head -30
done
log "prep_A done"

