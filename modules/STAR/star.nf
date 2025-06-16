 #!/usr/bin/env nextflow

/*
 * Star processes for creating an index
 */
 
 process STARindex {

    publishDir "${params.outdir}/ref/", mode: 'copy'

    input: 
    path genome
    path gtf

    output:
        path 'STAR_index', emit: index
    script:
    """
    mkdir -p STAR_index
    STAR --runMode genomeGenerate \\
            --genomeDir STAR_index \\
            --genomeFastaFiles ${genome} \\
            --sjdbGTFfile ${gtf} \\
            --runThreadN 4
    """
 }

