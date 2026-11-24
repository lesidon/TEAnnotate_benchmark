// Stage 3 (masking) — single-pass RepeatMasker over the ORIGINAL genome.
//
// One RepeatMasker run with the curated + classified consensus library (-lib) and
// RepeatMasker's built-in simple/low-complexity/satellite detection LEFT ON (no
// -nolow, no -noint). That is what makes the summary .tbl carry BOTH the interspersed
// TEs AND the tandem/simple/satellite repeats with their genome percentages — the
// earlier two-round + ProcessRepeats concatenation dropped the simple rows from the
// combined table. Flags:
//   -xsmall : soft-mask (lowercase) every masked base -> the *.masked deliverable.
//   -a      : write the alignment file (*.align), consumed by the Kimura landscape.
//   -no_is  : skip the bacterial insertion-sequence check (eukaryote-irrelevant, faster).
// The *.out table is split by class into three GFF3 views (simple-only, TE-only,
// combined) plus a combined BED. The complete *.tbl and *.align pass straight through.
//
// Biocontainer preferred, with bioconda::repeatmasker as the Conda fallback.
// --writable-tmpfs: RepeatMasker writes config/library files inside its (read-only)
// install dir on Singularity/Apptainer images.
process REPEATMASKER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(genome, stageAs: 'genome_in.fa'), path(te_library)

    output:
    tuple val(meta), path("${meta.id}.softmasked.fa")  , emit: masked
    tuple val(meta), path("${meta.id}.repeats.gff3")   , emit: gff        // combined (all classes)
    tuple val(meta), path("${meta.id}.simple.gff3")    , emit: simple_gff // simple / low-complexity / satellite
    tuple val(meta), path("${meta.id}.te.gff3")        , emit: te_gff     // interspersed TEs
    tuple val(meta), path("${meta.id}.repeats.bed")    , emit: bed
    tuple val(meta), path("${meta.id}.repeats.align")  , emit: align
    tuple val(meta), path("${meta.id}.repeats.out")    , emit: out
    tuple val(meta), path("${meta.id}.repeats.tbl")    , emit: tbl
    path 'versions.yml'                                , emit: versions

    script:
    // RepeatMasker -pa forks N parallel jobs that each use ~4 cores with the rmblast
    // engine, so divide cpus by 4 (min 1) to avoid 4x oversubscription.
    def pa = Math.max(1, (task.cpus as int).intdiv(4))
    """
    # \$HOME on this HPC is read-only from inside the container; redirect HOME/caches.
    export HOME="\$PWD"
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    mkdir -p "\$MPLCONFIGDIR" rmout

    # An empty TE library says nothing about simple repeats: RepeatMasker then runs in
    # simple/low-complexity-only mode (-noint) instead of being skipped.
    if [[ -s ${te_library} ]]; then
        LIB_ARGS=(-lib ${te_library})
    else
        echo "WARN [${meta.id}]: empty TE library; masking simple / low-complexity repeats only (-noint)" >&2
        LIB_ARGS=(-noint)
    fi

    # RepeatMasker builds its general libraries (is.lib etc.) with makeblastdb in
    # <Libraries>/general.working. Inside the image that directory lives in the
    # --writable-tmpfs overlay, and on hosts where Apptainer falls back to
    # fuse-overlayfs makeblastdb refuses to write there ("You do not have write
    # permissions on 'general.working'" -> "createLib(): Error invoking makeblastdb").
    # Point -libdir at a task-local directory of symlinks to the installed Libraries,
    # so that build happens in the work dir on every host.
    RM_LIBSRC="\$(dirname "\$(readlink -f "\$(command -v RepeatMasker)")")/Libraries"
    mkdir -p rmlib
    for f in "\${RM_LIBSRC}"/*; do ln -s "\$f" rmlib/; done

    # Default sensitivity; simple/low-complexity detection ON (no -nolow) so the .tbl
    # reports tandem + interspersed together. -xsmall soft-masks; -a writes the .align
    # for the divergence landscape.
    RM_RC=0
    RepeatMasker \\
        -libdir "\$PWD/rmlib" \\
        "\${LIB_ARGS[@]}" \\
        -xsmall -a -no_is \\
        -pa ${pa} \\
        -dir rmout \\
        genome_in.fa || RM_RC=\$?

    # A RepeatMasker failure must fail the task (Nextflow then retries it, see
    # conf/modules.config), never silently pass the genome through: the
    # makeblastdb/general.working failure above once produced an empty annotation
    # that was cached as a success and broke every downstream QC step.
    # RepeatMasker writes a .out even when it finds no repeats, so a missing .out
    # always means it did not run to completion.
    RMOUT=\$(find rmout -name '*.out' -print -quit 2>/dev/null || true)
    if [[ \${RM_RC} -ne 0 || -z "\${RMOUT}" ]]; then
        echo "ERROR [${meta.id}]: RepeatMasker failed (exit \${RM_RC}); no annotation produced." >&2
        exit 1
    fi

    MASKED=\$(find rmout -name '*.masked' -print -quit 2>/dev/null || true)
    if [[ -n "\${MASKED}" && -s "\${MASKED}" ]]; then
        cp "\${MASKED}" ${meta.id}.softmasked.fa
    else
        echo "WARN [${meta.id}]: RepeatMasker masked nothing; passing genome through" >&2
        cp -L genome_in.fa ${meta.id}.softmasked.fa
    fi

    ALIGN=\$(find rmout -name '*.align' -print -quit 2>/dev/null || true)
    [[ -n "\${ALIGN}" && -s "\${ALIGN}" ]] && cp "\${ALIGN}" ${meta.id}.repeats.align || true

    TBL=\$(find rmout -name '*.tbl' -print -quit 2>/dev/null || true)
    [[ -n "\${TBL}" && -s "\${TBL}" ]] && cp "\${TBL}" ${meta.id}.repeats.tbl || true

    OUT=\$(find rmout -name '*.out' -print -quit 2>/dev/null || true)
    [[ -n "\${OUT}" && -s "\${OUT}" ]] && cp "\${OUT}" ${meta.id}.repeats.out || true

    # ── Split the .out into simple / TE / combined GFF3 + a combined BED ──────────
    # .out: 3 header lines, then whitespace columns —
    #   1=SW score  5=query  6=begin  7=end  9=strand(+/C)  10=repeat  11=class/family.
    # A class starting Simple_repeat / Low_complexity / Satellite is a "simple" repeat;
    # everything else (LTR/LINE/SINE/DNA/RC/Unknown/…) is an interspersed TE.
    # Features are annotated regions, not alignments, so the repeat goes into Name and a
    # lower-case repeat_class attribute (GFF3 reserves Target for "id start end"); IDs
    # are unique per file and reserved characters in values are percent-encoded.
    write_gff () {  # args: bedfile  gff_out  feature_type
        {
            echo "##gff-version 3"
            awk -v ft="\$3" 'BEGIN{FS=OFS="\\t"}
                function esc(s) { gsub(/%/, "%25", s); gsub(/;/, "%3B", s); gsub(/=/, "%3D", s);
                                  gsub(/&/, "%26", s); gsub(/,/, "%2C", s); return s }
                { split(\$4, a, "#"); rep=esc(a[1]); cls=esc(a[2]);
                  print \$1,"RepeatMasker",ft,\$2+1,\$3,\$5,\$6,".",
                        "ID=" ft "_" NR ";Name=" rep ";repeat_class=" cls }' "\$1"
        } > "\$2"
    }

    if [[ -s ${meta.id}.repeats.out ]]; then
        # Full table -> tab file: chrom start0 end repeat#class score strand class
        awk 'BEGIN{OFS="\\t"} NR>3 && NF>=11 { s=(\$9=="C")?"-":"+";
            print \$5, \$6-1, \$7, \$10"#"\$11, \$1, s, \$11 }' \\
            ${meta.id}.repeats.out > repeats.all.tab || true
    else
        echo "WARN [${meta.id}]: no RepeatMasker .out table (no repeats?)" >&2
        : > repeats.all.tab
    fi

    # Combined BED (drop the trailing class column used only for splitting).
    cut -f1-6 repeats.all.tab > ${meta.id}.repeats.bed || : > ${meta.id}.repeats.bed

    # Simple / TE partitions (class is the 7th column of repeats.all.tab).
    awk 'BEGIN{FS=OFS="\\t"} \$7 ~ /^(Simple_repeat|Low_complexity|Satellite)/ {print \$1,\$2,\$3,\$4,\$5,\$6}' \\
        repeats.all.tab > simple.bed || : > simple.bed
    awk 'BEGIN{FS=OFS="\\t"} \$7 !~ /^(Simple_repeat|Low_complexity|Satellite)/ {print \$1,\$2,\$3,\$4,\$5,\$6}' \\
        repeats.all.tab > te.bed || : > te.bed

    write_gff ${meta.id}.repeats.bed ${meta.id}.repeats.gff3 repeat_region
    write_gff simple.bed             ${meta.id}.simple.gff3  tandem_repeat
    write_gff te.bed                 ${meta.id}.te.gff3      dispersed_repeat

    # Guarantee the optional artefacts exist (empty placeholder when no repeats), so
    # the visualization joins never silently drop a sample.
    for f in ${meta.id}.repeats.align ${meta.id}.repeats.out ${meta.id}.repeats.tbl; do
        [[ -f "\$f" ]] || : > "\$f"
    done

    RM_VER=\$(RepeatMasker -v 2>&1 | grep -oE '[0-9]+\\.[0-9.]+' | head -1 || true)
    [[ -z "\${RM_VER}" ]] && RM_VER=4.1.5
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \${RM_VER}
    END_VERSIONS
    """

    stub:
    """
    cp -L genome_in.fa ${meta.id}.softmasked.fa
    : > ${meta.id}.repeats.bed
    for g in repeats simple te; do echo "##gff-version 3" > ${meta.id}.\$g.gff3; done
    for f in ${meta.id}.repeats.align ${meta.id}.repeats.out ${meta.id}.repeats.tbl; do : > "\$f"; done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: "stub"
    END_VERSIONS
    """
}
