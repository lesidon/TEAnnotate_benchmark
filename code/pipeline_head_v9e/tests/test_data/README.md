# Test data

The test profile (`-profile test`) expects a small FASTA file at:

    tests/test_data/Neurospora_crassa_chr7_1Mb.fa

## Obtaining the test genome

Download and subset the N. crassa OR74A assembly (NC12, GCF_000182925.2):

```bash
# Requires: samtools or seqkit, wget/curl
curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz" \
  | zcat \
  | seqkit grep -p "NC_026510.1" \  # chromosome 7
  | seqkit subseq -r 1:1000000 \
  > tests/test_data/Neurospora_crassa_chr7_1Mb.fa
```

The resulting file is ~1 Mb and will be classified as `draft`
(N50 = 1 Mb < 10 Mb threshold), which is sufficient to exercise
all input_validation logic without triggering satellite/HOR Stage B.
