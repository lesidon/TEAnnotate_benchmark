#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TEANNOTATE } from './workflows/teannotate'

workflow {
    TEANNOTATE ()
}

workflow.onComplete {
    log.info (
        workflow.success
            ? "\nPipeline complete! Results: ${params.outdir}\n"
            : "\nPipeline failed — check .nextflow.log for details.\n"
    )
}
