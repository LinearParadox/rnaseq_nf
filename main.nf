#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for rna seq
 */
include { STARindex } from './modules/STAR/star.nf'
include { salmonIndex as salmonIndex } from './modules/salmon/salmon.nf'
include {salmonIndex as salmonTranscriptome} from './modules/salmon/salmon.nf'


workflow index {
    gtf = file(params.gtf)
    genome = file(params.genome)
    transcriptome = file(params.transcriptome)
    STAR_index=STARindex(genome, gtf)
    salmon_index=salmonIndex(genome, transcriptome, "salmon_genome")
    if ( params.transcript_level ){
        transcript_level_transcriptome = file(params.transcript_level_transcriptome)
        full_transcript = salmonTranscriptome(genome, transcript_level_transcriptome, "salmon_transcriptome")
    }
    }
workflow {
    samples = Channel.fromPath(params.samplesheet).splitCsv().view { row -> "${row[0]} - ${row[1]} - ${row[2]}" }
}
