// De novo TE library construction with HiTE (CSU-KangHu/HiTE, Nat Commun 2024).
// HiTE detects full-length TEs (LTR / TIR / Helitron / non-LTR) with dynamic
// boundary adjustment and emits a classified consensus library that is merged
// with the Earl Grey library before Stage 2 dedup. Run library-only (--annotate 0).
//
// Container-only: the self-contained image kanghu/hite:3.3.3 bundles the source at
// /HiTE plus all tool/database deps. HiTE has no bioconda package, and its conda
// environment.yml provides only dependencies (not the main.py source), so there is
// no usable Conda fallback — the image is required. Disable with --skip_hite.
process HITE {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}.hite_library.fa"), emit: library
    path 'versions.yml'                                , emit: versions

    script:
    // --plant: HiTE defaults to 1 ("Is it a plant genome") and passes it to its TIR
    // (judge_TIR_transposons.py) and non-LTR (judge_Non_LTR_transposons.py) detectors.
    // Left at the default, every animal and fungal genome ran with plant-tuned TIR and
    // non-LTR detection — the run log showed "Is plant genome = [ 1 ]" for Drosophila.
    // Derive it from the same organism_type that selects the de novo route.
    def is_plant = meta.organism_type == 'plant' ? 1 : 0
    //
    // LTR detection: HiTE's default LTR step (FiLTR) runs
    // bin/FiLTR-main/bin/LtrDetector/bin/LtrDetector, and kanghu/hite:3.3.3 ships
    // WITHOUT that whole bin/ directory (the image code matches HiTE commits
    // 14cb9e0..03c2e77, 2025-10-17/19, which do contain it). FiLTR then gets zero
    // candidates and, because main.py never checks that step's status, HiTE exits 0
    // with no LTR families at all. Supply the missing directory with
    // --hite_ltrdetector_dir (bound over the image path in conf/modules.config; see
    // the README for the three-step recipe). Do NOT switch FiLTR off instead: the
    // classic LTR_harvest/LTR_finder -> LTR_retriever path also fails in this image
    // ("The RMblast engine is not installed in RepeatMasker!").
    """
    # \$HOME on this HPC is read-only from inside the container; redirect HOME and
    # tool caches to the (writable) task work dir so RepeatMasker/matplotlib/etc.
    # can write their dotfiles.
    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    OUT="${meta.id}_hite"
    mkdir -p "\$OUT"

    python /HiTE/main.py \\
        --genome ${genome} \\
        --out_dir "\$PWD/\$OUT" \\
        --thread ${task.cpus} \\
        --plant ${is_plant} \\
        --annotate 0

    # HiTE runs its structural LTR step without checking the exit status, so a broken
    # LTR toolchain still "succeeds" with zero LTR families. Make that visible.
    if [[ ! -s "\$OUT/confident_ltr_cut.fa" ]]; then
        echo "WARN [${meta.id}]: HiTE found no LTR families (\$OUT/confident_ltr_cut.fa empty or missing)." >&2
        echo "  Its structural LTR step (step 4) most likely failed silently; see .command.log." >&2
        echo "  With kanghu/hite:3.3.3 this is the missing LtrDetector - set --hite_ltrdetector_dir." >&2
    fi

    # Final classified consensus library (drop-in for RepeatMasker -lib).
    LIB=\$(find "\$OUT" -name 'confident_TE.cons.fa' -print -quit 2>/dev/null || true)
    if [[ -n "\${LIB}" && -s "\${LIB}" ]]; then
        cp "\${LIB}" ${meta.id}.hite_library.fa
    else
        echo "WARN [${meta.id}]: HiTE produced no confident_TE.cons.fa" >&2
        touch ${meta.id}.hite_library.fa
    fi

    HITE_VER=\$(python /HiTE/main.py --version 2>/dev/null | grep -oiE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    [[ -z "\${HITE_VER}" ]] && HITE_VER=3.3.3
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hite: \${HITE_VER}
    END_VERSIONS
    """

    stub:
    """
    cat > ${meta.id}.hite_library.fa <<-EOF
    >${meta.id}_hite_stub#LTR/Gypsy
    GGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCT
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hite: "stub"
    END_VERSIONS
    """
}
