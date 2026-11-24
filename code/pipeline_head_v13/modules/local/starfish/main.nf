// Starship captain report for fungal genomes (on by default; off with --skip_starfish).
// Starship elements are giant (15-800 kb) fungal transposons carried by a tyrosine-
// recombinase "captain". This runs the single-genome part of the starfish workflow:
//   1. annotate - captain genes predicted de novo (metaeuk + the YR profile HMMs and
//                 reference captains bundled with starfish); no gene annotation needed
//   2. sketch   - genomic neighbourhoods around the captains
//   3. insert   - element boundaries. Starfish places them by aligning a captain region
//                 to an empty insertion site, normally in ANOTHER genome, so a single
//                 genome usually yields captains but no elements; any it does find are
//                 run through flank + summarize and reported too.
// The result is a REPORT only (<id>.starfish_captains.tsv and <id>_starfish/): it never
// enters the TE library, curation, classification or masking.
// Starfish splits genome and feature IDs on '_', so contigs are renamed <code>_<n> for
// the run (contig_map.tsv); the captain table gives the original contig names.
// The process is non-fatal (label error_ignore), and so is its environment: Starfish
// (egluckthaler channel; no image that works with Nextflow's task wrapper) runs from a
// Conda env that the task itself builds once (micromamba, mamba or conda on the host)
// and caches under --conda_cachedir (default <workDir>/conda), or from a prebuilt env
// given with --starfish_env. A Conda env declared to Nextflow would abort the whole
// run when it cannot be built; here a failed build only skips this report. The YR
// database ships inside the env at <env>/db.
process STARFISH_ANNOTATE {
    tag "$meta.id"
    label 'process_low'
    label 'error_ignore'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.starfish_captains.tsv"), emit: captains
    path "${meta.id}_starfish/"                              , emit: output_dir
    path 'versions.yml'                                      , emit: versions

    script:
    // genome code without the separator (and anything else Starfish could split on)
    def ome        = meta.id.replaceAll(/[^A-Za-z0-9]/, '') ?: 'genome'
    def env_yml    = "${projectDir}/assets/environments/starfish.yml"
    def cache_root = params.conda_cachedir ?: "${workflow.workDir}/conda"
    """
    OUT="${meta.id}_starfish"
    mkdir -p "\${OUT}" geneFinder elementFinder blastdb
    printf 'captain\\tcontig\\tstart\\tend\\tstrand\\tmetaeuk_score\\tbest_reference_captain\\tneighbourhood\\n' \\
        > ${meta.id}.starfish_captains.tsv
    printf '"%s":\\n    starfish: "not run"\\n' "${task.process}" > versions.yml
    skip_report() {
        echo "WARN [${meta.id}]: \$1; Starship report skipped, the run continues" >&2
        exit 0
    }

    # ── Starfish environment: prebuilt (--starfish_env) or built here once ─────
    SF_ENV="${params.starfish_env ?: ''}"
    if [[ -n "\${SF_ENV}" ]]; then
        [[ -x "\${SF_ENV}/bin/starfish" ]] || skip_report "--starfish_env \${SF_ENV} has no bin/starfish"
    else
        SF_ENV="${cache_root}/starfish-\$(md5sum < ${env_yml} | cut -c1-12)"
        if [[ ! -f "\${SF_ENV}/.teannotate_ok" ]]; then
            SOLVER=\$(command -v micromamba || command -v mamba || command -v conda || true)
            [[ -n "\${SOLVER}" ]] || skip_report "no micromamba, mamba or conda on PATH to build the Starfish environment"
            mkdir -p "${cache_root}"
            build_env() {
                [[ -f "\${SF_ENV}/.teannotate_ok" ]] && return 0      # built meanwhile by another task
                rm -rf "\${SF_ENV}"
                echo "INFO [${meta.id}]: building the Starfish environment \${SF_ENV} with \${SOLVER}" >&2
                if [[ "\$(basename "\${SOLVER}")" == micromamba ]]; then
                    "\${SOLVER}" create -y -p "\${SF_ENV}" -f ${env_yml} > env_build.log 2>&1
                else
                    "\${SOLVER}" env create -p "\${SF_ENV}" -f ${env_yml} > env_build.log 2>&1
                fi
            }
            # one build at a time when several fungal samples start together
            if command -v flock > /dev/null; then
                ( flock -w 10800 9 && build_env ) 9> "\${SF_ENV}.lock" && BUILT=1 || BUILT=0
            else
                build_env && BUILT=1 || BUILT=0
            fi
            if [[ \${BUILT} == 1 && -x "\${SF_ENV}/bin/starfish" ]]; then
                touch "\${SF_ENV}/.teannotate_ok"
            else
                rm -rf "\${SF_ENV}"
                tail -n 20 env_build.log >&2 2>/dev/null || true
                skip_report "building the Starfish environment failed (env_build.log)"
            fi
        fi
    fi
    export PATH="\${SF_ENV}/bin:\${PATH}" CONDA_PREFIX="\${SF_ENV}"

    # ── YR (tyrosine-recombinase) reference DB bundled with starfish ───────────
    DB="\${SF_ENV}/db"
    REFS="\${DB}/YRsuperfamRefs.faa"
    HMM=\$(find "\${DB}" -name 'YRsuperfams*.hmm' 2>/dev/null | head -1 || true)
    if [[ ! -f "\${REFS}" || -z "\${HMM}" ]]; then
        skip_report "no starfish YR database (YRsuperfamRefs.faa / YRsuperfams*.hmm) under \${DB}"
    fi

    SF_VER=\$(ls "\${SF_ENV}"/conda-meta/starfish-*.json 2>/dev/null | sed -E 's@.*/starfish-([0-9][^-]*)-.*@\\1@' | head -1 || true)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        starfish: \${SF_VER:-unknown}
    END_VERSIONS

    # ── Contigs renamed <code>_<n> (map kept for the report) ──────────────────
    zcat -f ${fasta} | awk -v o="${ome}" '
        /^>/ { n++; id = substr(\$1, 2); print id "\\t" o "_" n > "contig_map.tsv"; print ">" o "_" n; next }
        { print }' > ${ome}.fna
    printf '%s\\t%s\\n' "${ome}" "\$(pwd)/${ome}.fna" > ome2assembly.txt
    cp contig_map.tsv "\${OUT}/"

    # ── 1. Captain genes, de novo ──────────────────────────────────────────────
    starfish annotate -T ${task.cpus} -x ${ome}_tyr -a ome2assembly.txt \\
        -p "\${HMM}" -P "\${REFS}" -i tyr -o geneFinder/

    CAP_GFF=geneFinder/${ome}_tyr.filt.gff
    CAP_IDS=geneFinder/${ome}_tyr.filt.ids
    if [[ ! -s "\${CAP_IDS}" ]]; then
        echo "INFO [${meta.id}]: no captain (tyrosine recombinase) genes - no Starships" >&2
        exit 0
    fi
    cp "\${CAP_GFF}" "\${OUT}/${meta.id}.captains.gff"
    cp geneFinder/${ome}_tyr.filt.fas "\${OUT}/${meta.id}.captains.faa" 2>/dev/null || true

    # ── 2. Neighbourhoods around the captains ─────────────────────────────────
    printf '%s\\t%s\\n' "${ome}" "\$(pwd)/\${CAP_GFF}" > ome2gff.txt
    starfish sketch -m 10000 -q "\${CAP_IDS}" -g ome2gff.txt -i s -x ${ome} -o geneFinder/
    awk -F'\\t' '\$5 == "tyr"' geneFinder/${ome}.bed > geneFinder/${ome}_tyr.bed
    cp geneFinder/${ome}.bed "\${OUT}/${meta.id}.neighbourhoods.bed"

    # Captain table with the original contig names
    awk -F'\\t' -v OFS='\\t' '
        FILENAME == ARGV[1] { orig[\$2] = \$1; next }
        FILENAME == ARGV[2] { hood[\$4] = \$7; next }
        \$3 == "mRNA" {
            name = ""; target = ""
            k = split(\$9, a, ";")
            for (i = 1; i <= k; i++) {
                if (a[i] ~ /^Name=/)      name   = substr(a[i], 6)
                if (a[i] ~ /^Target_ID=/) target = substr(a[i], 11)
            }
            print name, ((\$1 in orig) ? orig[\$1] : \$1), \$4, \$5, \$7, \$6, target, ((name in hood) ? hood[name] : ".")
        }' contig_map.tsv geneFinder/${ome}_tyr.bed "\${CAP_GFF}" >> ${meta.id}.starfish_captains.tsv

    # ── 3. Element boundaries (usually none from a single genome) ─────────────
    makeblastdb -in ${ome}.fna -out blastdb/${ome}.assemblies -parse_seqids -dbtype nucl > /dev/null
    starfish insert -T ${task.cpus} -a ome2assembly.txt -d blastdb/${ome}.assemblies \\
        -b geneFinder/${ome}_tyr.bed -i tyr -x ${ome} -o elementFinder/
    if [[ -s elementFinder/${ome}.insert.bed ]]; then
        starfish flank -a ome2assembly.txt -b elementFinder/${ome}.insert.bed -x ${ome} -o elementFinder/
        starfish summarize -a ome2assembly.txt -b elementFinder/${ome}.flank.bed -x ${ome} -o elementFinder/ \\
            -S elementFinder/${ome}.insert.stats -f elementFinder/${ome}.flank.singleDR.stats \\
            -g ome2gff.txt -t "\${CAP_IDS}"
        mkdir -p "\${OUT}/elements"
        find elementFinder -maxdepth 1 -type f -exec cp {} "\${OUT}/elements/" \\;
    else
        echo "INFO [${meta.id}]: \$(( \$(wc -l < ${meta.id}.starfish_captains.tsv) - 1 )) captain(s); no element boundaries (single genome: starfish needs an empty insertion site, normally in another genome)" >&2
    fi
    """

    stub:
    """
    mkdir -p "${meta.id}_starfish"
    printf 'captain\\tcontig\\tstart\\tend\\tstrand\\tmetaeuk_score\\tbest_reference_captain\\tneighbourhood\\n' \\
        > ${meta.id}.starfish_captains.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        starfish: "stub"
    END_VERSIONS
    """
}
