// Starship element annotation for fungal genomes.
// Starship elements are giant (>100 kb) fungal transposons with a tyrosine-
// recombinase "captain" gene. This runs the real starfish single-genome
// workflow to extract Starship element sequences for the TE library:
//   1. annotate    - Gene Finder: locate captain (tyr) genes from the gene models
//   2. consolidate - merge captain genes into the input GFF
//   3. makeblastdb - self-BLAST db of the assembly (needed by insert)
//   4. insert      - Element Finder: define element boundaries around captains
//   5. summarize   - emit element metadata (.feat) and nucleotide sequences (.fna)
// Intermediate filenames vary across versions, so they are discovered with `find`.
// The process is non-fatal (label error_ignore): Starship annotation is a fungal
// bonus, so a failure (or a genome with no captains/elements) is skipped and the
// rest of the run continues with just the Earl Grey library.
// Runs via the egluckthaler Conda env (see conf/modules.config) — the ORAS image
// is incompatible with Nextflow's task wrapper (busybox shell). The YR db ships
// inside the env at $CONDA_PREFIX/db.
process STARFISH_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    input:
    tuple val(meta), path(fasta), path(gff)

    output:
    tuple val(meta), path("${meta.id}.starfish_library.fa"), emit: library
    path "${meta.id}_starfish/"                            , emit: output_dir
    path 'versions.yml'                                    , emit: versions

    script:
    """
    OUT="${meta.id}_starfish"
    # starfish requires the -o output directories to already exist (it does not
    # create them), so pre-create the geneFinder/elementFinder subdirectories.
    mkdir -p "\${OUT}/geneFinder" "\${OUT}/elementFinder"

    # Use the work-dir paths (bind-mounted in the container) for the control files.
    GENOME="\$(pwd)/${fasta}"
    GFF="\$(pwd)/${gff}"

    # ── Locate the bundled YR (tyrosine-recombinase) reference DB ──────────────
    # In the egluckthaler container the db ships at \$CONDA_PREFIX/db; fall back
    # to a search under /opt if CONDA_PREFIX is not exported in this shell.
    DB=""
    for cand in "\${CONDA_PREFIX:-}/db" "\${STARFISHDIR:-}/db"; do
        if [[ -f "\${cand}/YRsuperfamRefs.faa" ]]; then DB="\${cand}"; break; fi
    done
    if [[ -z "\${DB}" ]]; then
        DB=\$(dirname "\$(find "\${CONDA_PREFIX:-/opt}" /opt -name 'YRsuperfamRefs.faa' 2>/dev/null | head -1)")
    fi
    REFS="\${DB}/YRsuperfamRefs.faa"
    HMM=\$(find "\${DB}" -name 'YRsuperfams*.hmm' 2>/dev/null | head -1)
    if [[ ! -f "\${REFS}" || -z "\${HMM}" ]]; then
        echo "ERROR [${meta.id}]: could not locate starfish YR database (YRsuperfamRefs.faa / YRsuperfams*.hmm) under \${DB}" >&2
        exit 1
    fi
    echo "[${meta.id}] starfish DB: \${DB}" >&2

    # Emit versions early so the graceful "no Starships" exit paths still report it.
    SF_VER=\$(ls "\${CONDA_PREFIX:-}"/conda-meta/starfish-*.json /opt/conda/envs/*/conda-meta/starfish-*.json 2>/dev/null | sed -E 's@.*/starfish-([0-9][^-]*)-.*@\\1@' | head -1 || true)
    [[ -z "\${SF_VER}" ]] && SF_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        starfish: \${SF_VER}
    END_VERSIONS

    # ── Control files (ome -> assembly / gff) ──────────────────────────────────
    printf '%s\\t%s\\n' "${meta.id}" "\${GENOME}" > ome2assembly.txt
    printf '%s\\t%s\\n' "${meta.id}" "\${GFF}"     > ome2gff.txt

    # ── 1. Gene Finder: find captain (tyr) genes ───────────────────────────────
    # -n 'ID=' : AGAT/BRAKER store the feature name in the ID attribute, not the
    # default Name= field, so point starfish's nameField at ID= to read gene names.
    starfish annotate \\
        -T ${task.cpus} \\
        -x "${meta.id}_tyr" \\
        -a ome2assembly.txt \\
        -g ome2gff.txt \\
        -p "\${HMM}" \\
        -P "\${REFS}" \\
        -i tyr \\
        -n 'ID=' \\
        -o "\${OUT}/geneFinder/"

    CAPTAIN_GFF=\$(find "\${OUT}/geneFinder" -name '*.filt_intersect.gff' -print -quit 2>/dev/null || true)
    TYR_BED=\$(find "\${OUT}/geneFinder" -name '*.tyr.bed' -print -quit 2>/dev/null || true)

    # Genomes without captain (tyrosine-recombinase) genes have no Starships. Exit
    # gracefully with an empty library instead of cascading through consolidate/
    # insert/summarize with empty -G/--bed arguments.
    if [[ -z "\${CAPTAIN_GFF}" || ! -s "\${CAPTAIN_GFF}" || -z "\${TYR_BED}" || ! -s "\${TYR_BED}" ]]; then
        echo "INFO [${meta.id}]: No captain (tyrosine recombinase) genes found — no Starship elements" >&2
        touch "${meta.id}.starfish_library.fa"
        exit 0
    fi

    # ── 2. Consolidate captain genes into the input GFF ────────────────────────
    starfish consolidate -o ./ -g "\${GFF}" -G "\${CAPTAIN_GFF}"
    CONS_GFF=\$(find . -maxdepth 1 -name '*consolidated.gff' -print -quit)
    printf '%s\\t%s\\n' "${meta.id}" "\$(realpath "\${CONS_GFF}")" > ome2consolidatedGFF.txt

    # ── 3. Self-BLAST database of the assembly (required by insert) ────────────
    mkdir -p blastdb
    cut -f2 ome2assembly.txt | xargs cat > "blastdb/${meta.id}.assemblies.fna"
    makeblastdb -in "blastdb/${meta.id}.assemblies.fna" -out "blastdb/${meta.id}.assemblies" -parse_seqids -dbtype nucl

    # ── 4. Element Finder: define Starship boundaries around captains ──────────
    starfish insert \\
        -T ${task.cpus} \\
        -a ome2assembly.txt \\
        -d "blastdb/${meta.id}.assemblies" \\
        -b "\${TYR_BED}" \\
        -i tyr \\
        -x "${meta.id}" \\
        -o "\${OUT}/elementFinder/"

    INSERT_BED=\$(find "\${OUT}/elementFinder" -name '*.insert.bed'   -print -quit)
    INSERT_STATS=\$(find "\${OUT}/elementFinder" -name '*.insert.stats' -print -quit)

    # ── 5. Summarize: emit element nucleotide sequences (.fna) ─────────────────
    starfish summarize \\
        -a ome2assembly.txt \\
        -b "\${INSERT_BED}" \\
        -x "${meta.id}" \\
        -o "\${OUT}/elementFinder/" \\
        -S "\${INSERT_STATS}" \\
        -g ome2consolidatedGFF.txt

    # ── Collect element nucleotide sequences into the library ──────────────────
    ELEMENTS=\$(find "\${OUT}/elementFinder" -name '*.fna' -print -quit 2>/dev/null || true)
    if [[ -n "\${ELEMENTS}" ]]; then
        cp "\${ELEMENTS}" "${meta.id}.starfish_library.fa"
    else
        echo "INFO [${meta.id}]: No Starship elements detected" >&2
        touch "${meta.id}.starfish_library.fa"
    fi
    """

    stub:
    """
    mkdir -p "${meta.id}_starfish"
    cat > "${meta.id}.starfish_library.fa" <<-EOF
    >${meta.id}_stub_Starship#ClassII/Helitron
    GGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCT
    TAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCGGCTAGCG
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        starfish: "stub"
    END_VERSIONS
    """
}
