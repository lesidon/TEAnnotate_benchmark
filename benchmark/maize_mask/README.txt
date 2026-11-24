maize_mask/ - chunked full-genome RepeatMasker of maize for the held-out benchmark (benchmark only,
not part of the pipeline). Started 16 Nov 2026 so that masking runs while Zhores is still up.

Why: one whole-genome RepeatMasker run per library (edta_manual_bench/rm_one.sh) would take days per
library on the 2.2 Gb genome (about 45,000 CPU-h for all libraries). The genome is split into
~50 Mb chunks and many workers mask <library, chunk> items in parallel.

  split_genome.py   genome -> chunks/cNNN.fa, pieces.tsv, chunks.tsv. Long chromosomes are cut at
                    the middle of the nearest >=100-N gap within 2 Mb of an equal split (so no hit can
                    span a cut at a gap); short sequences are packed whole.
  methods.tsv       method -> library, in priority order (reference first: every score needs it).
  queue.py          file-based work queue (claims/done/attempts in queue/, mkdir mutex); a library is
                    snapshotted into libs/ when first used (and only when its file is >10 min old).
  mask_worker.sh    claims items whose estimated time fits its remaining wall time; RepeatMasker with
                    exactly the options of rm_one.sh (-e ncbi -q -no_is -norna -nolow -div 40
                    -cutoff 225, -libdir symlinks, same Dfam TE Tools 1.89.2 image), -pa = CPUs/4.
                    The worker that finishes a method's last chunk runs merge_out.py.
  merge_out.py      per-chunk .out -> merged/<method>/zea_mays.fna.out in genome coordinates
                    (query name, begin/end shifted, (left) recomputed, IDs offset per chunk).
  mask_feed.sh      feeder job: keeps up to 5 htc (48 CPUs, 1 day, gn nodes) and 2 mem (24 CPUs,
                    3 days, not ct02-04) workers queued while claimable items exist; leaves submit
                    slots for other jobs; exits when every item is done or failed 3 times.

Status:   python queue.py status        Logs: logs/worker_<job>.out, logs/feed_<job>.out
Scoring:  register the methods in edta_manual_bench (copy library.fa first, then the merged .out
          into edta_manual_bench/zea_mays/<method>/zea_mays.fna.out) and run the harmonised
          lib-test + family level as for the other held-out genomes.
