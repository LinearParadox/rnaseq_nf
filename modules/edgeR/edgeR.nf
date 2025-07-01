

 /*
 Gene level differential expression analysis using edgeR
 */
 
 process differential_expression{
    cpus 2
    memory 8.GB
    tag "Differential Expression"
    publishDir "${params.outdir}/differential_genes", mode: 'copy', pattern: "csv"
    publishDir "${params.outdir}/differential_genes/", mode: 'copy', pattern: "figs"
    publishDir "${params.outdir}/pipeline_info/", mode: 'copy', pattern: "sessionInfo_deg.txt"
    input:
    file salmon_files
    file script
    file index
    file samplesheet
    file transcriptome
    file gtf
    val organism
    file design
    file contrast_matrix
    output:
    path "figs", emit: figures
    path "csv", emit: csv
    path "sessionInfo.txt", emit: session_info
    script:
    """
    mkdir cache
    Rscript ${script} ${samplesheet} ${transcriptome} ${gtf} ${organism} ${design} ${contrast_matrix}

    """

 }

  process differential_transcripts{
    cpus 2
    memory 8.GB
    tag "Differential Expression"
    publishDir "${params.outdir}/differential_transcripts", mode: 'copy', pattern: "csv"
    publishDir "${params.outdir}/differential_transcripts/", mode: 'copy', pattern: "figs"
    publishDir "${params.outdir}/pipeline_info/", mode: 'copy', pattern: "sessionInfo_det.txt"
    input:
    file salmon_files
    file script
    file samplesheet
    val organism
    file design
    file contrast_matrix
    output:
    path "figs", emit: figures
    path "csv", emit: csv
    path "sessionInfo_det.txt", emit: session_info
    script:
    """
    mkdir cache
    Rscript ${script} ${samplesheet} ${organism} ${design} ${contrast_matrix}

    """

 }