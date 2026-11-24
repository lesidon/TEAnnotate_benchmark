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
    path ltrdetector   // LtrDetector directory (HITE_LTRDETECTOR or --hite_ltrdetector_dir)

    output:
    tuple val(meta), path("${meta.id}.hite_library.fa"), emit: library
    tuple val(meta), path("${meta.id}.hite_status.tsv"), emit: status
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
    // with no LTR families at all. HITE_LTRDETECTOR builds the missing directory (or
    // --hite_ltrdetector_dir supplies it), and HiTE runs from hite_src/: links to /HiTE,
    // with that directory replaced and the few code files that locate it copied (see
    // below). This works under every container engine without a bind mount.
    // LtrDetector is checked before discovery starts, and the run log must show FiLTR
    // using the copy (else the task fails). Do NOT switch FiLTR off instead: the
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

    # HiTE exists only inside its image (/HiTE). Without it (no engine, or an offline host
    # where the image could not be pulled) there is nothing to fall back to.
    if [[ ! -f /HiTE/main.py ]]; then
        echo "ERROR [${meta.id}]: /HiTE/main.py not found: HiTE runs only from its container image" >&2
        echo "  (kanghu/hite). Pull or cache the image (--singularity_cache_dir), or rerun with --skip_hite." >&2
        exit 1
    fi

    # hite_src/: /HiTE with bin/FiLTR-main/bin/LtrDetector replaced by the staged one.
    link_except() {   # <source dir> <target dir> <entry to leave out>
        mkdir -p "\$2"
        for e in "\$1"/* "\$1"/.[!.]*; do
            [[ -e "\$e" && "\$(basename "\$e")" != "\$3" ]] && ln -s "\$e" "\$2/"
        done
        return 0
    }
    # The code that decides where FiLTR looks is copied, not linked (< 100 kB): Python
    # puts a linked script's REAL directory on sys.path, so linked entry scripts would
    # import their modules, and derive their paths, from /HiTE again. FiLTR also builds
    # FiLTR-main/configs/../bin/LtrDetector, and '..' after a link leads back into
    # /HiTE. HiTE calls its helpers by name from PATH (/HiTE/tools:/HiTE/module in the
    # image), so the copies go first on PATH.
    HS="\$(pwd -P)/hite_src"   # physical path, as Python reports it
    link_except /HiTE                       "\$HS"                   bin
    rm "\$HS/main.py" "\$HS/module" "\$HS/tools"
    cp -r --preserve=mode /HiTE/main.py /HiTE/module /HiTE/tools "\$HS/"
    link_except /HiTE/bin                   "\$HS/bin"               FiLTR-main
    link_except /HiTE/bin/FiLTR-main        "\$HS/bin/FiLTR-main"    ""
    rm "\$HS/bin/FiLTR-main/bin" "\$HS/bin/FiLTR-main/configs" "\$HS/bin/FiLTR-main/main.py"
    cp -r --preserve=mode /HiTE/bin/FiLTR-main/configs /HiTE/bin/FiLTR-main/main.py "\$HS/bin/FiLTR-main/"
    link_except /HiTE/bin/FiLTR-main/bin    "\$HS/bin/FiLTR-main/bin" LtrDetector
    ln -s "\$(readlink -f ${ltrdetector})" "\$HS/bin/FiLTR-main/bin/LtrDetector"
    export PATH="\$HS/tools:\$HS/module:\$PATH"

    # Preflight: a working LtrDetector prints its usage; a missing binary or library
    # (exit 126/127) would otherwise surface hours later as an LTR-less library. Its
    # usage call exits 1 even when it works, so test the output, not the status.
    LTRD="\$HS/bin/FiLTR-main/bin/LtrDetector/bin/LtrDetector"
    LTRD_OUT=\$("\$LTRD" -help 2>&1 || true)
    if ! grep -q -- '-fasta' <<< "\$LTRD_OUT"; then
        echo "ERROR [${meta.id}]: LtrDetector does not run in this image:" >&2
        head -3 <<< "\$LTRD_OUT" >&2
        echo "  Without it HiTE finds no LTR families. Check --hite_ltrdetector_dir, or rerun with --skip_hite." >&2
        exit 1
    fi

    OUT="${meta.id}_hite"
    mkdir -p "\$OUT"

    # --work_dir: HiTE's scratch defaults to /tmp, where coarse_boundary.py keeps its
    # all-vs-all BLAST of the longest repeats (~370 GB for the human genome). On nodes
    # whose /tmp is a RAM-backed tmpfs that fills the node, starves every other job on
    # it, and HiTE then exits 0 with no library. Keep the scratch in the task dir.
    mkdir -p hite_work
    python "\$HS/main.py" \\
        --genome ${genome} \\
        --out_dir "\$PWD/\$OUT" \\
        --work_dir "\$PWD/hite_work" \\
        --thread ${task.cpus} \\
        --plant ${is_plant} \\
        --annotate 0 2>&1 | tee hite.log
    rm -rf hite_work

    # HiTE does not check the status of its steps: a full disk only shows in the log.
    if grep -aq 'No space left on device' hite.log; then
        echo "ERROR [${meta.id}]: HiTE ran out of disk space (see hite.log); its library would be incomplete." >&2
        exit 1
    fi

    # FiLTR must have run LtrDetector from the copy: from /HiTE it finds nothing, and
    # HiTE would still exit 0. If the log does not show the call at all, LTR discovery
    # is unconfirmed: recorded, not fatal (a genome may also truly lack LTRs).
    LTRD_HOME=\$(grep -aoE -- '--LtrDetector_home [^ ]+' hite.log | head -1 | cut -d' ' -f2 || true)
    if [[ -n "\${LTRD_HOME}" && "\${LTRD_HOME}" != "\$HS"/* ]]; then
        echo "ERROR [${meta.id}]: HiTE ran LtrDetector from \${LTRD_HOME}, not from the copy in \$HS." >&2
        exit 1
    fi
    LTR_STEP=ran
    if [[ -z "\${LTRD_HOME}" ]]; then
        LTR_STEP=unconfirmed
        echo "WARN [${meta.id}]: HiTE's log does not show FiLTR running LtrDetector; its LTR discovery" >&2
        echo "  is unconfirmed (${meta.id}.hite_status.tsv)." >&2
    fi
    N_LTR=\$(grep -c '^>' "\$OUT/confident_ltr_cut.fa" 2>/dev/null || true); N_LTR=\${N_LTR:-0}

    # HiTE runs its structural LTR step without checking the exit status, so a broken
    # LTR toolchain still "succeeds" with zero LTR families. Make that visible.
    if [[ "\${N_LTR}" -eq 0 ]]; then
        echo "WARN [${meta.id}]: HiTE found no LTR families (\$OUT/confident_ltr_cut.fa empty or missing)." >&2
        echo "  LtrDetector passed the preflight, so its structural LTR step (step 4) found" >&2
        echo "  nothing or failed silently; see .command.log." >&2
    fi
    printf 'ltr_step\\t%s\\nltr_sequences\\t%s\\n' "\${LTR_STEP}" "\${N_LTR}" > ${meta.id}.hite_status.tsv

    # Final classified consensus library (drop-in for RepeatMasker -lib).
    LIB=\$(find "\$OUT" -name 'confident_TE.cons.fa' -print -quit 2>/dev/null || true)
    if [[ -n "\${LIB}" && -s "\${LIB}" ]]; then
        cp "\${LIB}" ${meta.id}.hite_library.fa
    else
        echo "WARN [${meta.id}]: HiTE produced no confident_TE.cons.fa" >&2
        touch ${meta.id}.hite_library.fa
    fi

    HITE_VER=\$(python "\$HS/main.py" --version 2>/dev/null | grep -oiE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    [[ -z "\${HITE_VER}" ]] && HITE_VER=3.3.3
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hite: \${HITE_VER}
    END_VERSIONS
    """

    stub:
    """
    printf 'ltr_step\\tstub\\nltr_sequences\\t0\\n' > ${meta.id}.hite_status.tsv
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
