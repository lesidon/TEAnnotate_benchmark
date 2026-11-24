include { INPUT_VALIDATION                      } from '../subworkflows/local/input_validation/main'
include { DENOVO_TE_ANIMAL                      } from '../subworkflows/local/denovo_te_animal/main'
include { DENOVO_TE_FUNGUS                      } from '../subworkflows/local/denovo_te_fungus/main'
include { DENOVO_TE_PLANT                       } from '../subworkflows/local/denovo_te_plant/main'
include { LIBRARY_DEDUP                         } from '../subworkflows/local/library_dedup/main'
include { TE_EVIDENCE                           } from '../modules/local/te_evidence/main'
include { PFAM_FETCH                            } from '../modules/local/pfam_fetch/main'
include { HITE_LTRDETECTOR                      } from '../modules/local/hite_ltrdetector/main'
include { CURATION                             } from '../subworkflows/local/curation/main'
include { MASKING                              } from '../subworkflows/local/masking/main'
include { VISUALIZATION                        } from '../subworkflows/local/visualization/main'
include { FINAL_LIBRARY                        } from '../modules/local/final_library/main'
include { validateParameters; samplesheetToList } from 'plugin/nf-schema'
include { suppliedAssets                        } from '../modules/local/utils/asset_key'

workflow TEANNOTATE {

    main:
    ch_versions  = Channel.empty()
    ch_manifests = Channel.empty()   // source.tsv of every fetched asset

    // Validate all params against nextflow_schema.json before touching channels.
    validateParameters()

    // HiTE has no Conda recipe: without a container engine it cannot run, so stop now
    // rather than after hours of other discovery work. A -stub run never uses the image.
    if (!params.skip_hite && !workflow.containerEngine && !workflow.stubRun) {
        error "HiTE runs only from its container image (kanghu/hite), and this run has no " +
              "container engine (-profile conda / mamba). Use -profile singularity, apptainer, " +
              "docker or podman, or add --skip_hite to run without HiTE (the discovery set changes)."
    }

    // Samplesheet paths: absolute, or relative to the launch directory, or relative to the
    // samplesheet itself (so a samplesheet and its genomes can be moved together).
    def sheet_dir = file(params.input).parent
    def resolveInput = { String p, String what, String id ->
        if (!p) return null
        def f = file(p)
        if (!f.exists() && !p.startsWith('/') && !(p ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/)) {
            def alt = file("${sheet_dir}/${p}")
            if (alt.exists()) f = alt
        }
        if (!f.exists()) {
            error "Sample '${id}': ${what} '${p}' not found (looked in the launch directory and next to the samplesheet ${params.input})."
        }
        f
    }

    // ── Parse and schema-validate samplesheet (nf-schema) ─────────────────────
    Channel
        .fromList(
            samplesheetToList(
                params.input,
                "${projectDir}/assets/samplesheet_schema.json"
            )
        )
        .map { row ->
            def meta = [
                id            : row[0].id,
                organism_type : (row[0].organism_type ?: 'animal').toLowerCase(),
            ]
            [ meta, resolveInput(row[1].toString(), 'genome', meta.id) ]
        }
        .set { ch_input }

    // ── Stage 0: Input validation + genome statistics (QC report) ─────────────
    INPUT_VALIDATION(ch_input)
    ch_versions = ch_versions.mix(INPUT_VALIDATION.out.versions)

    // ── Stage 1: Organism-aware de novo TE discovery ──────────────────────────
    // Branch the validated genome channel by organism type. Anything that is not
    // explicitly fungus or plant (including an unspecified type) takes the animal route.
    INPUT_VALIDATION.out.genome
        .branch { meta, fasta ->
            fungus: meta.organism_type == 'fungus'
            plant : meta.organism_type == 'plant'
            animal: true
        }
        .set { ch_by_organism }

    // HiTE's LtrDetector (missing from its image): built once and cached, or a
    // prebuilt directory via --hite_ltrdetector_dir. Shared by every route.
    ch_ltrdetector = Channel.value([])
    if (!params.skip_hite) {
        if (params.hite_ltrdetector_dir) {
            ch_ltrdetector = Channel.value(file(params.hite_ltrdetector_dir, type: 'dir', checkIfExists: true))
        } else {
            HITE_LTRDETECTOR()
            ch_ltrdetector = HITE_LTRDETECTOR.out.dir.first()
            ch_versions    = ch_versions.mix(HITE_LTRDETECTOR.out.versions)
            ch_manifests   = ch_manifests.mix(HITE_LTRDETECTOR.out.manifest)
        }
    }

    // Animal (default route): Earl Grey + HiTE
    DENOVO_TE_ANIMAL(ch_by_organism.animal, ch_ltrdetector)
    ch_versions = ch_versions.mix(DENOVO_TE_ANIMAL.out.versions)

    // Fungus: Earl Grey + HiTE (+ the opt-in Starfish captain report)
    DENOVO_TE_FUNGUS(ch_by_organism.fungus, ch_ltrdetector)
    ch_versions = ch_versions.mix(DENOVO_TE_FUNGUS.out.versions)

    // Plant: EDTA (+ HELIANO) + HiTE
    DENOVO_TE_PLANT(ch_by_organism.plant, ch_ltrdetector)
    ch_versions = ch_versions.mix(DENOVO_TE_PLANT.out.versions)

    // Merge raw libraries from all implemented routes
    ch_raw_library = Channel.empty()
        .mix(DENOVO_TE_ANIMAL.out.raw_library)
        .mix(DENOVO_TE_FUNGUS.out.raw_library)
        .mix(DENOVO_TE_PLANT.out.raw_library)

    ch_raw_library
        .view { meta, lib ->
            "RAW_LIBRARY | sample=${meta.id} | organism=${meta.organism_type} | ${lib}"
        }

    // Which tool each raw sequence came from (cross-tool support in TE_EVIDENCE)
    ch_provenance = Channel.empty()
        .mix(DENOVO_TE_ANIMAL.out.provenance)
        .mix(DENOVO_TE_FUNGUS.out.provenance)
        .mix(DENOVO_TE_PLANT.out.provenance)

    // TE_EVIDENCE, TEtrimmer and the dedup representative choice work against the
    // genome: re-attach the validated genome by sample id.
    ch_genome_by_id = INPUT_VALIDATION.out.genome.map { meta, fasta -> [ meta.id, fasta ] }

    // ── Stage 2: Library deduplication & merge (CD-HIT-EST + user library) ─────
    LIBRARY_DEDUP(ch_raw_library, ch_genome_by_id)
    ch_versions = ch_versions.mix(LIBRARY_DEDUP.out.versions)

    ch_dedup_library = LIBRARY_DEDUP.out.library
    ch_dedup_library
        .view { meta, lib ->
            "DEDUP_LIBRARY | sample=${meta.id} | organism=${meta.organism_type} | ${lib}"
        }

    // Pfam-A: mandatory for TEtrimmer, and the TE-domain search of TE_EVIDENCE.
    // Fetched + indexed once (cached), or --pfam_db for a local copy; a value channel
    // so the single DB broadcasts to every sample.
    def need_pfam = !params.skip_tetrimmer || (params.evidence_gate != 'off' && params.evidence_pfam)
    ch_pfam = Channel.value([])
    if (need_pfam) {
        if (params.pfam_db) {
            ch_pfam = Channel.value(file(params.pfam_db, type: 'dir', checkIfExists: true))
        } else {
            PFAM_FETCH()
            ch_pfam      = PFAM_FETCH.out.db.first()
            ch_versions  = ch_versions.mix(PFAM_FETCH.out.versions)
            ch_manifests = ch_manifests.mix(PFAM_FETCH.out.manifest)
        }
    }

    // ── Stage 2 (gate): keep consensi with evidence of being TEs ──────────────
    // TE protein domains, the same element in another tool's library, or >=2
    // full-length genomic copies; non-TE repeats go to a separate file. Off with
    // --evidence_gate off.
    ch_evidence = Channel.empty()
    if (params.evidence_gate != 'off') {
        ch_dedup_library
            .map { meta, lib -> [ meta.id, meta, lib ] }
            .join(ch_raw_library.map { meta, raw -> [ meta.id, raw ] }, by: 0)
            .join(ch_genome_by_id, by: 0)
            .join(LIBRARY_DEDUP.out.lineage.map { meta, lin -> [ meta.id, lin ] }, by: 0)
            .join(ch_provenance.map { meta, prov -> [ meta.id, prov ] }, by: 0, remainder: true)
            .filter { it[2] != null }
            .map { id, meta, lib, raw, genome, lin, prov -> [ meta, lib, raw, genome, prov ?: [], lin ] }
            .set { ch_evidence_input }

        TE_EVIDENCE(ch_evidence_input, ch_pfam)
        ch_versions      = ch_versions.mix(TE_EVIDENCE.out.versions.first())
        ch_dedup_library = TE_EVIDENCE.out.library
        ch_evidence      = TE_EVIDENCE.out.evidence
    }

    // ── Stage 3: Curation (TEtrimmer) ─────────────────────────────────────────
    ch_dedup_library
        .map { meta, lib -> [ meta.id, meta, lib ] }
        .join(ch_genome_by_id, by: 0)
        .map { id, meta, lib, genome -> [ meta, lib, genome ] }
        .set { ch_curation_input }

    CURATION(ch_curation_input, ch_evidence, LIBRARY_DEDUP.out.lineage, ch_pfam)
    ch_versions  = ch_versions.mix(CURATION.out.versions)
    ch_manifests = ch_manifests.mix(CURATION.out.manifests)

    ch_classified_library = CURATION.out.library

    // The library masking uses, published as <id>.te_library.fa on every route.
    FINAL_LIBRARY(ch_classified_library)
    ch_versions = ch_versions.mix(FINAL_LIBRARY.out.versions.first())

    // ── Stage 3 (final): Genome masking ───────────────────────────────────────
    // Single-pass RepeatMasker soft-masking: one run with the classified library and
    // built-in simple/low-complexity detection ON, so the .tbl reports simple +
    // interspersed together; -a writes the .align used by the visualization stage.
    // Re-attach the validated genome to the classified library by id.
    // INPUT_VALIDATION.out.genome is a broadcast process output (safe to reuse).
    if (!params.skip_masking) {
        ch_classified_library
            .map { meta, lib -> [ meta.id, meta, lib ] }
            .join(INPUT_VALIDATION.out.genome.map { meta, fasta -> [ meta.id, fasta ] }, by: 0)
            .map { id, meta, lib, genome -> [ meta, genome, lib ] }
            .set { ch_masking_input }

        MASKING(ch_masking_input)
        ch_versions = ch_versions.mix(MASKING.out.versions)

        MASKING.out.genome_masked
            .view { meta, fasta ->
                "MASKED_GENOME | sample=${meta.id} | organism=${meta.organism_type} | ${fasta}"
            }

        // ── Stage 4: Visualization & QC (Nature-styled) ───────────────────────
        // Kimura divergence landscape (from the RepeatMasker .align), TE-density
        // karyoplot (from the combined BED + genome), and a repeat detection / masking
        // QC figure (from the .tbl + BED + classified library). Each figure ships the
        // R script that drew it. Re-attach the validated genome by sample id.
        if (!params.skip_visualization) {
            VISUALIZATION(
                MASKING.out.repeats_align,
                MASKING.out.repeats_out,
                INPUT_VALIDATION.out.genome,
                MASKING.out.repeats_bed,
                MASKING.out.repeats_tbl,
                ch_classified_library
            )
            ch_versions = ch_versions.mix(VISUALIZATION.out.versions)
        }
    }

    // One software-version manifest for the run (per-process fragments, de-duplicated).
    ch_versions
        .map { it.text.trim() }
        .unique()
        .collectFile(
            name     : 'software_versions.yml',
            storeDir : "${params.outdir}/pipeline_info",
            newLine  : true,
            sort     : true,
            seed     : "\"Workflow\":\n    nextflow: ${workflow.nextflow.version}\n" +
                       "    teannotate: ${workflow.manifest.version ?: 'dev'} (${workflow.commitId ?: 'local checkout'})"
        )

    // Assets: every fetched one (its cache directory, then its source.tsv rows:
    // requested and resolved identity, checksum) and every one given as a local path.
    ch_manifests
        .map { f -> f.readLines().findAll { it }.collect { "${f.parent}\t${it}" }.join('\n') }
        .mix(Channel.fromList(suppliedAssets()))
        .collectFile(name: 'asset_manifest.tsv', storeDir: "${params.outdir}/pipeline_info",
                     newLine: true, sort: true, seed: 'asset\tfield\tvalue')

    emit:
    versions = ch_versions
}
