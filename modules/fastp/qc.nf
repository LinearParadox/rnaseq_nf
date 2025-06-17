#! /usr/bin/env nextflow
params.outdir = 'results'
/*
 * Run fastp on a set of fastq files.
 */
process fastP{
    tag "FASTP on $sample"
    publishDir "${params.outdir}/qc/${sample}/fastp.html", mode: 'copy', pattern: "*.html"
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    tuple val(sample), path("fastp_${sample}.json"), emit: json
    tuple val(sample), path("${r1.simpleName}_trimmed.fastq.gz"), path(${r2.simpleName}_trimmed.fastq.gz) emit: reads
    script:
    """
    fastp -i ${r1} -I ${r2} -o ${r1.simpleName}_trimmed.fastq.gz -O ${r2.simpleName}_trimmed.fastq.gz -j fastp_${sample}.json -h ${sample}_fastp.html --detect_adapter_for_pe --length_required 20 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --max_len1 200 --max_len2 200
    """
}
