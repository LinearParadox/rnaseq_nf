 #!/usr/bin/env nextflow


/*
 * Star processes for creating an index
 */
 
 process STARindex {
    tag "STAR index"
    label 'star'

    publishDir "${params.outdir}/ref/", mode: 'copy'

    cpus 16
    memory 64.GB

    input: 
    path genome
    path gtf
    val read_length

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
            --runThreadN ${task.cpus} \\
            --genomeDir STAR_index \\
            --genomeFastaFiles "${decompressed_genome}" \\
            --sjdbGTFfile "${decompressed_gtf}" \\
            --sjdbOverhang ${read_length}-1 \\
    """
 }

process STARalign {
    cpus 16
    memory 64.GB
    label 'star'
    tag "STAR align on $sample"

    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    path gtf
    output:
    tuple val(sample), path("Aligned.out.bam"), emit: bam
    path "Log.final.out", emit: log

    script:
    decompressed_gtf = gtf.name.replaceAll(/\.gz$/, '')
    """
    if [[ "${gtf}" == *.gz ]]; then
        gunzip -f ${gtf}
    fi
    STAR --genomeDir ${index} \\
         --readFilesIn ${r1} ${r2} \\
         --runThreadN ${task.cpus} \\
         --twoPassMode Basic \\
         --readFilesCommand zcat \\
         --sjdbGTFfile ${decompressed_gtf} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFilterType BySJout    //reduces the number of "spurious" junctions
         --outFilterMultimapNmax 20         //max number of multiple alignments allowed for a read: if exceeded, the read is consi>
         --alignSJoverhangMin 8          //min overhang for unannotated junctions
         --alignSJDBoverhangMin 1          //min overhang for annotated junctions
         --outFilterMismatchNmax 999        //max number of mismatches per pair (absolute)
         --outFilterMismatchNoverLmax 0.06       //max number of mismatches per pair relative to read length: for 2x100b, max number of>
         --alignIntronMin 20         //min intron
         --alignIntronMax 1000000    //max intron
         --alignMatesGapMax 1000000    //max genomic distance between pairs
    """
}
