include { INPUT_VALIDATION                      } from '../subworkflows/local/input_validation/main'
include { DENOVO_TE_ANIMAL                      } from '../subworkflows/local/denovo_te_animal/main'
include { DENOVO_TE_FUNGUS                      } from '../subworkflows/local/denovo_te_fungus/main'
include { DENOVO_TE_PLANT                       } from '../subworkflows/local/denovo_te_plant/main'
include { LIBRARY_DEDUP                         } from '../subworkflows/local/library_dedup/main'
include { TE_EVIDENCE                           } from '../modules/local/te_evidence/main'
include { PFAM_FETCH                            } from '../modules/local/pfam_fetch/main'
include { CURATION                             } from '../subworkflows/local/curation/main'
include { MASKING                              } from '../subworkflows/local/masking/main'
include { VISUALIZATION                        } from '../subworkflows/local/visualization/main'
include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

workflow TEANNOTATE {

    main:
    ch_versions = Channel.empty()

    // Validate all params against nextflow_schema.json before touching channels.
    validateParameters()

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
                id              : row[0].id,
                organism_type   : (row[0].organism_type ?: 'animal').toLowerCase(),
                assembly_quality: params.assembly_quality,
                genome_size_mb  : null,
                num_seqs        : null,
                n50_mb          : null,
                gc_percent      : null,
            ]
            // Optional GFF3 gene annotation (last non-meta column); [] when absent.
            def gff = (row.size() > 2 && row[2]) ? file(row[2]) : []
            [ meta, file(row[1]), gff ]
        }
        .set { ch_samples }

    // Genome channel for validation / discovery; gene-annotation channel keyed by id.
    // The gene annotation is only used by the fungal Starfish route, which handles
    // any GTF/GFF -> GFF3 normalization itself (only when needed).
    ch_samples.map { meta, genome, gff -> [ meta, genome ] }.set { ch_input }
    ch_samples.map { meta, genome, gff -> [ meta.id, gff ] }.set { ch_gene_annotation }

    // ── Stage 0: Input validation + genome statistics ─────────────────────────
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

    // Animal (default route): Earl Grey + HiTE
    DENOVO_TE_ANIMAL(ch_by_organism.animal)
    ch_versions = ch_versions.mix(DENOVO_TE_ANIMAL.out.versions)

    // Fungus: Earl Grey + HiTE + Starfish (Starfish needs a GFF; attach it).
    ch_by_organism.fungus
        .map { meta, fasta -> [ meta.id, meta, fasta ] }
        .join(ch_gene_annotation, by: 0)
        .map { id, meta, fasta, gff -> [ meta, fasta, gff ] }
        .set { ch_fungus_input }

    DENOVO_TE_FUNGUS(ch_fungus_input)
    ch_versions = ch_versions.mix(DENOVO_TE_FUNGUS.out.versions)

    // Plant: EDTA (+ HELIANO) + HiTE
    DENOVO_TE_PLANT(ch_by_organism.plant)
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
            ch_pfam     = PFAM_FETCH.out.db.first()
            ch_versions = ch_versions.mix(PFAM_FETCH.out.versions)
        }
    }

    // ── Stage 2 (gate): keep consensi with evidence of being TEs ──────────────
    // TE protein domains, the same element in another tool's library, or >=2
    // full-length genomic copies; non-TE repeats go to a separate file. Off with
    // --evidence_gate off.
    ch_evidence = Channel.empty()
    if (params.evidence_gate != 'off') {
        def ch_trusted = params.user_te_library ? file(params.user_te_library, checkIfExists: true) : []
        ch_dedup_library
            .map { meta, lib -> [ meta.id, meta, lib ] }
            .join(ch_raw_library.map { meta, raw -> [ meta.id, raw ] }, by: 0)
            .join(ch_genome_by_id, by: 0)
            .join(ch_provenance.map { meta, prov -> [ meta.id, prov ] }, by: 0, remainder: true)
            .filter { it[2] != null }
            .map { id, meta, lib, raw, genome, prov -> [ meta, lib, raw, genome, prov ?: [] ] }
            .set { ch_evidence_input }

        TE_EVIDENCE(ch_evidence_input, ch_trusted, ch_pfam)
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

    CURATION(ch_curation_input, ch_evidence, ch_pfam)
    ch_versions = ch_versions.mix(CURATION.out.versions)

    ch_classified_library = CURATION.out.library

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

    emit:
    versions = ch_versions
}
