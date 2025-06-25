#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for rna seq
 */
include { STARindex } from './modules/STAR/star.nf'
include { STARalign } from './modules/STAR/star.nf'
include { qc_samples } from './workflows/qc_workflow.nf'
include { salmon_quant } from './modules/salmon/salmon.nf'

workflow {
    samples=Channel.fromPath(params.samplesheet).splitCsv().map {
        def sample = it[0]
        def r1 = file(it[1])
        def r2 = file(it[2])
        return [sample, r1, r2]
    } | groupTuple
    qc = qc_samples(samples)
    star_index = STARindex(file(params.genome), file(params.gtf), params.readlength)
    align = STARalign(qc.trimmed, star_index.index, file(params.gtf))
    salmon_quant = salmon_quant(align.bam, file(params.gtf), params.gibbs_sampling, 
                                params.seq_bias, params.gc_bias, params.dump_eq)
}