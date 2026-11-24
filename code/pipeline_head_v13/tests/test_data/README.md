# Test data

`zymoseptoria_tritici_chr21.fa.gz` is the test genome of `-profile test`
(`assets/test_samplesheet.csv`): chromosome 21 of *Zymoseptoria tritici* IPO323,
RefSeq `NC_018198.1` from assembly `GCF_000219625.1` (MG2), 409,213 bp. It is one of
the species' accessory chromosomes, which are short and rich in transposable elements:
the curated Dfam library of the species (Baril & Croll 2023) masks about 153 kb (37%)
of it, so every discovery, curation and classification step has real work to do.

| File | MD5 |
|---|---|
| `zymoseptoria_tritici_chr21.fa.gz` | `1ede9f4a6b7d77bbbadd1c2df2adde43` |
| uncompressed | `33111a648e76960171f571534fe2499c` |

## Regenerating it

```bash
curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/219/625/GCF_000219625.1_MYCGR_v2.0/GCF_000219625.1_MYCGR_v2.0_genomic.fna.gz" \
  | gzip -cd > zt.fna
{
  awk '/^>/ { p = (substr($1, 2) == "NC_018198.1"); if (p) print; next }' zt.fna
  awk '/^>/ { p = (substr($1, 2) == "NC_018198.1"); next } p' zt.fna | tr -d '\n' | fold -w 80
  echo
} > zymoseptoria_tritici_chr21.fa
gzip -n -9 zymoseptoria_tritici_chr21.fa
md5sum zymoseptoria_tritici_chr21.fa.gz
```

## Running the test

```bash
nextflow run . -profile test,singularity -stub --outdir results_stub   # wiring only, minutes
nextflow run . -profile test,singularity --outdir results_test         # every tool
```

The genome path in the samplesheet is relative to the samplesheet, so both commands work
from any directory (point `nextflow run` at the checkout or the repository).
