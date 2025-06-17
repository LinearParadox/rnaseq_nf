 #!/usr/bin/env nextflow


 /*
 * Decompress files if they are gzipped
 */



/*
 * Star processes for creating an index
 */
 
 process STARindex {

    publishDir "${params.outdir}/ref/", mode: 'copy'

    input: 
    file genome
    file gtf
    output:
    path 'STAR_index', emit: index
    script:
    decompressed_genome = genome.name.replaceAll(/\.gz$/, '')
    decompressed_gtf = gtf.name.replaceAll(/\.gz$/, '')
    """
    mkdir -p STAR_index
    if [[ "${genome}" == *.gz ]]; then
        gunzip -f ${genome}
    fi
    if [[ "${gtf}" == *.gz ]]; then
        gunzip -f ${gtf}
    fi
    STAR --runMode genomeGenerate \\
            --genomeDir STAR_index \\
            --genomeFastaFiles "${decompressed_genome}" \\
            --sjdbGTFfile "${decompressed_gtf}" \\
            --runThreadN 4
    """
 }

