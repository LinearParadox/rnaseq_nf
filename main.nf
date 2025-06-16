#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for rna seq
 */
include { STARindex } from './modules/STAR/star.nf'

workflow {
    // Generate a STAR index
    STARindex(file(params.genome), file(params.gtf))
}

