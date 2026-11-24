// Stage 3 (curation) — TEtrimmer.
// Automates manual curation of the de novo consensus library: for each consensus
// it BLASTs against the genome, extends/extracts hits, builds and cleans an MSA,
// trims boundaries, and emits a curated, de-duplicated consensus library.
// Classification is left to the dedicated downstream CLASSIFIER step, so
// --classify_all is intentionally NOT passed (it would invoke RepeatClassifier,
// which needs a configured Dfam database).
// --dedup (TEtrimmer's own CD-HIT on its input) is also NOT passed (default off):
// the input is already deduplicated class-aware in Stage 2, and a second, class-blind
// collapse would undo that. TEtrimmer's final merge cannot be switched off, but with
// --curation_class_aware its output is bypassed (TETRIMMER_POSTPROCESS reads the
// un-merged TEtrimmer_consensus.fasta).
process TETRIMMER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(library), path(genome)
    path pfam_dir

    output:
    tuple val(meta), path("${meta.id}.curated_library.fa"), emit: library
    tuple val(meta), path("${meta.id}_tetrimmer/")         , emit: output_dir, optional: true
    path 'versions.yml'                                    , emit: versions

    script:
    """
    # \$HOME (often a read-only HPC home) is not writable from inside the container,
    # so point matplotlib/fontconfig/numba caches at the (writable) task work dir.
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export FONTCONFIG_PATH="\$PWD/.fontconfig"
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME" "\$FONTCONFIG_PATH" "\$NUMBA_CACHE_DIR"

    # TE-Aid (TEtrimmer's terminal-repeat / full-length check) downloads RepeatPeps.lib
    # from a URL that now returns 404, so every TE-Aid call failed and TEtrimmer dropped
    # nearly all low-copy families ("check_low_copy is False"). The image directory is
    # not writable on every host (no overlay), so run TEtrimmer from a copy of its own
    # source (~350 kb) in the work dir: it finds TE-Aid relative to its own file, so
    # the copy's TE-Aid-master gets a db/ with the RepeatPeps.lib shipped with the
    # image's RepeatMasker (formatted here; the image's own index files are from mixed
    # BLAST builds). The copy is writable, so TEtrimmer's chmod of TE-Aid also works
    # without --writable-tmpfs.
    TETRIMMER_CMD=TEtrimmer
    TT_PY=\$(grep -oE '/[^ ]*/TEtrimmer\\.py' "\$(command -v TEtrimmer)" 2>/dev/null | head -1 || true)
    RM_LIBS="\$(dirname "\$(readlink -f "\$(command -v RepeatMasker)")" 2>/dev/null)/Libraries"
    if [[ -f "\${TT_PY}" ]]; then
        rm -rf tetrimmer_src
        cp -r "\$(dirname "\${TT_PY}")" tetrimmer_src && chmod -R u+w tetrimmer_src
        if [[ -s "\${RM_LIBS}/RepeatPeps.lib" && ! -s tetrimmer_src/TE-Aid-master/db/RepeatPeps.lib.phr ]]; then
            mkdir -p tetrimmer_src/TE-Aid-master/db
            makeblastdb -in "\${RM_LIBS}/RepeatPeps.lib" -dbtype prot \\
                -out tetrimmer_src/TE-Aid-master/db/RepeatPeps.lib > /dev/null || true
        fi
        TETRIMMER_CMD="python \$PWD/tetrimmer_src/\$(basename "\${TT_PY}")"
    fi
    [[ -s tetrimmer_src/TE-Aid-master/db/RepeatPeps.lib.phr ]] || \\
        echo "WARN [${meta.id}]: TE-Aid has no RepeatPeps.lib; TEtrimmer will drop most low-copy families" >&2

    # --pfam_dir points at a pre-downloaded, hmmpress-indexed Pfam dir (PFAM_FETCH).
    # PFAM is mandatory: TEtrimmer aborts curation without it. The dir is staged
    # read-only; since it is already indexed, TEtrimmer reads but never writes it.
    \${TETRIMMER_CMD} \\
        --input_file  ${library} \\
        --genome_file ${genome} \\
        --output_dir  ${meta.id}_tetrimmer \\
        --pfam_dir    ${pfam_dir} \\
        --num_threads ${task.cpus}

    # Final curated library is the de-duplicated consensus; fall back to the
    # pre-dedup consensus if the merged file is absent (older/edge versions).
    if [[ -s ${meta.id}_tetrimmer/TEtrimmer_consensus_merged.fasta ]]; then
        RAW_CONS=${meta.id}_tetrimmer/TEtrimmer_consensus_merged.fasta
    elif [[ -s ${meta.id}_tetrimmer/TEtrimmer_consensus.fasta ]]; then
        RAW_CONS=${meta.id}_tetrimmer/TEtrimmer_consensus.fasta
    else
        RAW_CONS=
    fi

    if [[ -n "\${RAW_CONS}" ]]; then
        # ── Drop low-confidence "Need_check" consensi before re-classification ────
        # TEtrimmer sorts each consensus into a confidence tier under
        # <out>/.../{Perfect,Good,Reco_check,Need_check}. Sequences that landed in a
        # Need_check folder are the ones flagged for manual inspection and must not
        # reach the classifier. Build the set of Need_check ids from those folders,
        # then keep a FASTA record only if its id is NOT in the set AND its header
        # does not itself carry a need_check tag (belt-and-braces for versions that
        # annotate the status inline). The awk record filter is multi-line-safe.
        : > need_check.ids
        for d in \$(find ${meta.id}_tetrimmer -type d -iname 'need_check' 2>/dev/null); do
            find "\${d}" -maxdepth 1 -type f -printf '%f\\n' 2>/dev/null
        done | sed -E 's/\\.[^.]+\$//' | sort -u > need_check.ids || true

        # Keep every record whose id is NOT a Need_check id and whose header carries
        # no inline need_check tag. Read the id list in BEGIN (getline) rather than the
        # NR==FNR two-file idiom: when need_check.ids is EMPTY, NR==FNR stays true into
        # the FASTA and swallows the whole library (empty output). getline is robust to
        # an empty/absent id list — nothing to drop => keep everything.
        awk -v badfile=need_check.ids '
             BEGIN { while ((getline line < badfile) > 0) if (line != "") bad[line]=1 }
             /^>/ { id=\$0; sub(/^>/,"",id); sub(/[[:space:]#].*\$/,"",id);
                    keep = !(id in bad) && (tolower(\$0) !~ /need_?check/) }
             keep' "\${RAW_CONS}" > ${meta.id}.curated_library.fa

        DROPPED=\$( (wc -l < need_check.ids) 2>/dev/null || echo 0)
        echo "INFO [${meta.id}]: removed \${DROPPED} Need_check consensus id(s) before classification" >&2
        [[ -s ${meta.id}.curated_library.fa ]] || \\
            echo "WARN [${meta.id}]: curated library empty after Need_check filtering" >&2
    else
        echo "WARN [${meta.id}]: TEtrimmer produced no consensus library" >&2
        touch ${meta.id}.curated_library.fa
    fi

    TT_VER=\$(ls "\${CONDA_PREFIX:-}"/conda-meta/tetrimmer-*.json /usr/local/conda-meta/tetrimmer-*.json 2>/dev/null | sed -E 's@.*/tetrimmer-([0-9][^-]*)-.*@\\1@' | head -1 || true)
    [[ -z "\${TT_VER}" ]] && TT_VER=unknown
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetrimmer: \${TT_VER}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_tetrimmer
    cp ${library} ${meta.id}.curated_library.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetrimmer: "stub"
    END_VERSIONS
    """
}
