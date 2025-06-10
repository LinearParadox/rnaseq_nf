#! /usr/bin/env nextflow
params.outdir = 'results'
/*
 * Run fastp on a set of fastq files.
 */
process fastP{
    tag "FASTP on $sample"
    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(sample), path(reads)
    output:
    
}
