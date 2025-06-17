#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for rna seq
 */
include { STARindex } from './modules/STAR/star.nf'
include { salmonIndex } from './modules/salmon/salmon.nf'


workflow {
    gtf = file(params.gtf)
    genome = file(params.genome)
    transcriptome = file(params.transcriptome)
    STARindex(genome, gtf)
    salmonIndex(genome, transcriptome)
}
