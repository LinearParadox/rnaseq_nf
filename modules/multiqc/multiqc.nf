process multiqc{
    cpus 1
    memory 4.GB
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    input:
    path fastp_logs
    path salmon_logs
    path star_logs
    path "csv"
    output:
    path "multiqc_report.html", emit: report
    script:
    """
    multiqc .
    """

}