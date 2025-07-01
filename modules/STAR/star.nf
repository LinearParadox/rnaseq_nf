 #!/usr/bin/env nextflow


/*
 * Star processes for creating an index
 */
 
 process STARindex {
    tag "STAR index"
    label 'star'

    publishDir "${params.outdir}/ref/", mode: 'copy'

    cpus 16
    memory 72.GB

    input: 
    file genome
    file gtf
    val read_length
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
            --sjdbOverhang \$(( "${read_length}-1" ))
    """
 }

process STARalign {
    cpus 16
    memory 64.GB
    label 'star'
    tag "STAR align"

    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*.bam"
    publishDir "${params.outdir}/pipeline_info/", mode: "copy", pattern: "STAR_version.txt"

    input:
    tuple val(sample), path(r1), path(r2)
    path index 
    path gtf
    output:
    tuple val(sample), path("*.bam"), emit: bam
    path "${sample}_Log.final.out", emit: log
    path "STAR_version.txt", emit: star_version

    script:
    decompressed_gtf = gtf.name.replaceAll(/\.gz$/, '')
    """
    STAR | grep "^STAR version" > STAR_version.txt
    if [[ "${gtf}" == *.gz ]]; then
        gunzip -f ${gtf}
    fi
    STAR --genomeDir ${index} \\
         --readFilesIn ${r1} ${r2} \\
         --runThreadN ${task.cpus} \\
         --twopassMode Basic \\
         --readFilesCommand zcat \\
         --sjdbGTFfile ${decompressed_gtf} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFilterType BySJout \\
         --outFilterMultimapNmax 20 \\
         --alignSJoverhangMin 8 \\
         --alignSJDBoverhangMin 1 \\
         --outFilterMismatchNmax 999 \\
         --outFilterMismatchNoverLmax 0.06 \\
         --alignIntronMin 20 \\
         --alignIntronMax 1000000 \\
         --alignMatesGapMax 1000000
    mv  Log.final.out ${sample}_Log.final.out
    """
}
