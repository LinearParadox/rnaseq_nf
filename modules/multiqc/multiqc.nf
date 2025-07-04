process multiqc{
    cpus 1
    memory 4.GB
    tag "MultiQC"
    publishDir "${params.outputDir}/", mode: 'copy'
    input:
    path fastp_logs
    path salmon_logs
    path star_logs
    output:
    path "multiqc_report.html", emit: report
    script:
    """
    echo 'use_filename_as_sample_name:
        - fastp' > multiqc_config.yaml
    multiqc -c multiqc_config.yaml .

    """

}