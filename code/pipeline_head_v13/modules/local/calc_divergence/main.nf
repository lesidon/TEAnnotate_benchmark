// Stage 4 (visualization) — class-by-divergence matrix for the Kimura landscape.
//
// Primary: CpG-adjusted Kimura distances via calcDivergenceFromAlign.pl (ships with
// RepeatMasker) on the .align, extracting its "coverage for each repeat class and
// divergence" matrix. Fallback: if that yields nothing (util missing, empty/absent
// .align), bin the RepeatMasker .out perc-divergence column directly into the same
// class x divergence (bp) matrix — so a landscape is produced whenever repeats exist.
// That matrix's first column is "PercDiv" instead of "Div", so the plot labels its
// axis as plain divergence, not CpG-adjusted Kimura distance.
//
// Same RepeatMasker biocontainer as masking (bioconda fallback); --writable-tmpfs
// because the utility touches files inside the read-only install dir.
process CALC_DIVERGENCE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(rm_out), path(align)

    output:
    tuple val(meta), path("${meta.id}.kimura_matrix.tsv"), emit: matrix
    tuple val(meta), path("${meta.id}.divsum")           , emit: divsum, optional: true
    path 'versions.yml'                                  , emit: versions

    script:
    """
    export HOME="\$PWD"

    # ── Primary: CpG-adjusted Kimura matrix from the .align ──────────────────────
    CDA=\$(command -v calcDivergenceFromAlign.pl || true)
    [[ -z "\${CDA}" ]] && CDA=\$(find / -name calcDivergenceFromAlign.pl 2>/dev/null | head -1 || true)

    # calcDivergenceFromAlign.pl uses RepeatMasker's own Perl modules (SearchResult.pm,
    # etc.), but some images (the biocontainer, some conda builds) don't put the
    # RepeatMasker install dir on @INC -> "Can't locate SearchResult.pm". Find the module
    # and add its directory to PERL5LIB so the util loads regardless of the container.
    SRPM=\$(find / -name SearchResult.pm -path '*RepeatMasker*' 2>/dev/null | head -1 || true)
    [[ -n "\${SRPM}" ]] && export PERL5LIB="\$(dirname "\${SRPM}"):\${PERL5LIB:-}"

    if [[ -n "\${CDA}" && -s ${align} ]]; then
        perl "\${CDA}" -s ${meta.id}.divsum ${align} 1>&2 || true
        if [[ -s ${meta.id}.divsum ]]; then
            awk 'f && NF>0 { print }
                 /Coverage for each repeat class and divergence/ { f=1 }' \\
                ${meta.id}.divsum > ${meta.id}.kimura_matrix.tsv || true
        fi
    fi

    # ── Fallback: bin the .out perc-div column into a class x divergence matrix ───
    # .out cols: 2=perc div, 6=begin, 7=end, 11=class/family. Top-level class before '/'.
    if [[ ! -s ${meta.id}.kimura_matrix.tsv && -s ${rm_out} ]]; then
        echo "INFO [${meta.id}]: Kimura matrix from .align empty; binning .out perc-div" >&2
        awk 'BEGIN{OFS="\\t"}
             NR>3 && NF>=11 { d=int(\$2+0.5); split(\$11,c,"/"); cl=c[1]; if(cl=="")cl="Unknown";
                              bp=\$7-\$6+1; if(bp<0)bp=-bp; cov[d SUBSEP cl]+=bp;
                              cls[cl]=1; if(d>maxd)maxd=d }
             END{ n=0; for(x in cls){ n++; name[n]=x }
                  printf "PercDiv"; for(i=1;i<=n;i++) printf "\\t%s", name[i]; printf "\\n";
                  for(d=0; d<=maxd; d++){ printf "%d", d;
                      for(i=1;i<=n;i++){ k=d SUBSEP name[i]; printf "\\t%d", (k in cov)?cov[k]:0 }
                      printf "\\n" } }' \\
            ${rm_out} > ${meta.id}.kimura_matrix.tsv || true
    fi
    [[ -s ${meta.id}.kimura_matrix.tsv ]] || : > ${meta.id}.kimura_matrix.tsv

    RM_VER=\$(RepeatMasker -v 2>&1 | grep -oE '[0-9]+\\.[0-9.]+' | head -1 || true)
    [[ -z "\${RM_VER}" ]] && RM_VER=4.1.5
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \${RM_VER}
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.kimura_matrix.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: "stub"
    END_VERSIONS
    """
}
