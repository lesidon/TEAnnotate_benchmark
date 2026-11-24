#!/usr/bin/env python3
"""How large is the Unknown fraction the ensemble works on, per v9 library, and how much of
the genome do those consensi mask (v9 RepeatMasker .out of the bench)?"""
import os
B = "/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark"
E = f"{B}/edta_manual_bench"
RUNS = {  # sp: (v9 run root, library dir prefix)
    "arabidopsis_thaliana": (f"{B}/work/refined_dedup_v9", "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "drosophila_melanogaster": (f"{B}/work/refined_dedup_v9", "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "zymoseptoria_tritici": (f"{B}/work/refined_dedup_v9", "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "septoria_linicola": (f"{B}/work/refined_dedup_v9", "pipeline_refined_v9", "pipeline_wo_tetrimmer_refined_v9"),
    "caenorhabditis_briggsae": (f"{B}/work/heldout_v9", "pipeline_v9", "pipeline_wo_tetrimmer_v9"),
}
for sp, (root, lw, lo) in RUNS.items():
    for arm, libdir, bdir in (("out_with", lw, "teannotate_ref9"), ("out_wo", lo, "teannotate_wo_ref9")):
        lib = f"{B}/libraries/{sp}/{libdir}/pipeline_library.fa"
        vt = f"{root}/{sp}/{arm}/{sp}/05_classification/{sp}.voting_table.tsv"
        labels = {}
        for l in open(lib):
            if l.startswith(">"):
                n, _, lab = l[1:].split()[0].partition("#")
                labels[n] = lab
        voted = {l.split("\t")[0]: l.rstrip("\n").split("\t")[-1] for i, l in enumerate(open(vt)) if i}
        unk = {n for n, lab in labels.items() if lab.lower().startswith("unknown")}
        U = unk | set(voted)
        relab = {n for n, c in voted.items() if c != "Unknown"}
        out = f"{E}/{sp}/{bdir}/{sp}.fna.out"
        tot = bp_u = bp_rel = 0
        if os.path.exists(out):
            for i, l in enumerate(open(out)):
                f = l.split()
                if i < 3 or len(f) < 11:
                    continue
                L = int(f[6]) - int(f[5]) + 1
                tot += L
                if f[9] in U: bp_u += L
                if f[9] in relab: bp_rel += L
        print(f"{sp:24s} {arm:8s} library {len(labels):5d}; Unknown-at-vote {len(U):4d} "
              f"(v9 relabelled {len(relab):3d}, still Unknown {len(unk):4d}); masked by them "
              f"{bp_u/1e6:6.2f} Mb of {tot/1e6:7.2f} Mb ({100*bp_u/max(1,tot):4.1f}%), "
              f"v9-relabelled {bp_rel/1e6:5.2f} Mb")

