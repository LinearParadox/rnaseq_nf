process multiqc{
    cpus 1
    memory 4.GB
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    input:
    path "figs"
    path "csv"
    output:
    path "multiqc_report.html", emit: report
    script:
    """
    multiqc ${figs} ${csv} -o multiqc_report.html
    """

}