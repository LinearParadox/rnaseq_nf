#! /usr/bin/env nextflow
/*
 * Run fastp on a set of fastq files.
 */

 process merge_lanes{
    cpus 2
    memory '8 GB'
    tag "Merging lanes for $sample"
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    tuple val(sample), path("R1_merged.fastq.gz"), path("R2_merged.fastq.gz"), emit: reads
    script:
    """
    cat ${r1} > R1_merged.fastq.gz
    cat ${r2} > R2_merged.fastq.gz
    """



 }
process fastP{
    cpus fastp_threads
    tag "FASTP on $sample"
    publishDir "${params.outdir}/per-sample-outs/${sample}/fastp.html", mode: 'copy', pattern: "*.html"
    input:
    tuple val(sample), path(r1), path(r2)
    val fastp_threads 
    output:
    path "*fastp.{json,html}", emit: fastp_results
    tuple val(sample), path("R1_trimmed.fastq.gz"), path("R2_trimmed.fastq.gz"), emit: reads
    script:
    """
    fastp -i ${r1} -I ${r2} -o R1_trimmed.fastq.gz -O R2_trimmed.fastq.gz -j fastp_${sample}.json -h ${sample}_fastp.html --detect_adapter_for_pe -w ${fastp_threads}
    """
}
