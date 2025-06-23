#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for rna seq
 */
include { salmon_index } from './index_workflows.nf'
include { salmon_index as salmon_index_transcript_level } from './index_workflows.nf'
include { star_index } from './index_workflows.nf'
include { fastP } from './modules/fastp/qc.nf'

workflow {
    /*
    salmon_genome_index = salmon_index(file(params.genome), file(params.transcriptome), "salmon_gene_level")
    if (params.transcript_level) {
        salmon_transcript_index = salmon_index_transcript_level(file(params.genome), file(params.transcript_level_transcriptome), "salmon_transcript_level")
    }
    else {
        salmon_transcript_index = null
    }
    */
    samples = Channel.fromPath(params.samplesheet).splitCsv()
    fastp_results = fastP(samples)

}
