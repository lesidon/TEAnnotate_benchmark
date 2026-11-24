"""vote_v10_common.py - genomes, arms and paths of the v10-vote replay on frozen v9 runs."""
B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E = f"{B}/edta_manual_bench"
V = f"{B}/vote_v10"
DEV = f"{B}/work/refined_dedup_v9"
HO = f"{B}/work/heldout_v9"

# sp: (code, organism_type, v9 run root, library dir with / w/o TEtrimmer)
GENOMES = {
    "arabidopsis_thaliana":      ("ath", "plant",  DEV, "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "drosophila_melanogaster":   ("dme", "animal", DEV, "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "zymoseptoria_tritici":      ("ztr", "fungus", DEV, "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "septoria_linicola":         ("sli", "fungus", DEV, "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "caenorhabditis_briggsae":   ("cbr", "animal", HO,  "pipeline_v9", "pipeline_wo_tetrimmer_v9"),
    "caenorhabditis_elegans":    ("cel", "animal", HO,  "pipeline_v9", "pipeline_wo_tetrimmer_v9"),
    "blumeria_hordei":           ("bho", "fungus", HO,  "pipeline_v9", "pipeline_wo_tetrimmer_v9"),
    # v9e = v9 + EDTA sequence-ID fix (16d03c3); v9 failed in EDTA on this genome
    "chlamydomonas_reinhardtii": ("chl", "plant",  f"{B}/work/heldout_v9e", "pipeline_v9", "pipeline_wo_tetrimmer_v9"),
    "ectocarpus_sp7":            ("ect", "animal", HO,  "pipeline_v9", "pipeline_wo_tetrimmer_v9"),
}
# arm: (run outdir, library dir index in GENOMES, bench method dir under E/<sp>/)
ARMS = {"w": ("out_with", 3, "teannotate_ref9"), "o": ("out_wo", 4, "teannotate_wo_ref9")}
# default kingdom model sets of the two-step vote (nextflow.config, commit 27320e2)
SUBSET = {"plant": "teclass2,create,neuralte,terrier", "animal": "teclass2,neuralte,terrier",
          "fungus": "create,neuralte,terrier"}


def read_fasta(path):
    recs, name, seq = [], None, []
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if name is not None:
                recs.append((name, seq))
            name, seq = l[1:], []
        else:
            seq.append(l)
    if name is not None:
        recs.append((name, seq))
    return recs


def paths(sp, arm):
    code, org, root, lw, lo = GENOMES[sp]
    outd, li, bdir = ARMS[arm]
    lib = f"{B}/libraries/{sp}/{GENOMES[sp][li]}/pipeline_library.fa"
    vt = f"{root}/{sp}/{outd}/{sp}/05_classification/{sp}.voting_table.tsv"
    return code + arm, org, lib, vt, f"{E}/{sp}/{bdir}/{sp}.fna.out"

