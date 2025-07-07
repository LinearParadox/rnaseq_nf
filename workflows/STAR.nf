include { STARindex } from '../modules/STAR/star.nf'
include { STARalign } from '../modules/STAR/star.nf'

workflow star {
    take:
        samples 
        star_index 
        gtf 
        genome
        read_length
    main:
        if (file(star_index).exists()) {
        star_index = channel.fromPath(star_index).collect()
        } else {
        star_index = STARindex(file(genome), file(gtf), read_length).index.collect()
        }
        align = STARalign(samples, star_index, file(gtf))
    emit:
        bam = align.bam
        starlog = align.log
        }
