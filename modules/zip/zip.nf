process zip_files {
    cpus 4
    memory 16.GB
    publishDir "${params.outdir}/", mode: 'copy', pattern: "output.zip"
    input:
    path de_csv, stageAs: "differential_genes/csv"
    path de_figs, stageAs: "differential_genes/figs"
    path dt_csv, stageAs: "differential_transcripts/csv"
    path dt_figs, stageAs: "differential_transcripts/figs"
    path multiqc, stageAs: "multiqc_report.html"
    output:
    path "output.zip"
    script:
    """
    zip -r output.zip differential_genes/ differential_transcripts/ multiqc_report.html
    """
}