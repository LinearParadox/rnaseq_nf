#!/usr/bin/env nextflow

/*  
 * Workflows for indexing
 */
include { STARindex } from './modules/STAR/star.nf'
include { salmonIndex as salmonIndex } from './modules/salmon/salmon.nf'


workflow salmon_index {
    take:
        genome
        transcriptome
        out_name
    main:
        salmon_index=salmonIndex(genome, transcriptome, out_name)

    emit: 
        salmon_index= salmon_index
        }
workflow star_index {
    take:
        genome
        gtf

    main:
        gtf = gtf
        genome = genome
        STAR_index=STARindex(genome, gtf)

    emit: 
        star_index= STAR_index
}