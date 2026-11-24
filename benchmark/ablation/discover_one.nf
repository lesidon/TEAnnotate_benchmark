#!/usr/bin/env nextflow
// discover_one.nf - run ONE de novo discovery tool of the pipeline on one genome, with the
// pipeline's own module (pipeline_head_v11, unchanged: all run-time fixes included), so
// every tool x genome gets an identical, separately timed run. Used for the discovery-tool
// ablation (which tools per kingdom, accuracy vs compute) and as the standalone baselines.
//   nextflow run discover_one.nf -c $H/nextflow.config -c ablation.config -profile singularity \
//       --tool hite|earlgrey|edta|heliano --genome G --id SP --organism_type plant|animal|fungus \
//       --hite_ltrdetector_dir DIR --outdir OUT
nextflow.enable.dsl = 2

def H = '/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/pipeline_head_v11'
include { HITE     } from '/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/pipeline_head_v11/modules/local/hite/main'
include { EARLGREY } from '/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/pipeline_head_v11/modules/local/earlgrey/main'
include { EDTA     } from '/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/pipeline_head_v11/modules/local/edta/main'
include { HELIANO  } from '/trinity/home/leonid.sidorov/tmp/teannotate_benchmark/whole_pipeline_benchmark/pipeline_head_v11/modules/local/heliano/main'

workflow {
    def meta = [ id: params.id, organism_type: params.organism_type ]
    ch = Channel.of([ meta, file(params.genome, checkIfExists: true) ])
    switch (params.tool) {
        case 'hite':
            HITE(ch, Channel.value(file(params.hite_ltrdetector_dir, type: 'dir', checkIfExists: true))); break
        case 'earlgrey': EARLGREY(ch); break
        case 'edta':     EDTA(ch);     break
        case 'heliano':  HELIANO(ch);  break
        default: error "--tool must be hite, earlgrey, edta or heliano"
    }
}
