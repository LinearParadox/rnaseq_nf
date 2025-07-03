

 /*
 Gene level differential expression analysis using edgeR
 */
 
 process differential_expression{
    cpus 2
    memory 16.GB
    tag "Differential Expression"
    label "edgeR"
    publishDir "${params.outdir}/differential_genes", mode: 'copy', pattern: "csv"
    publishDir "${params.outdir}/differential_genes/", mode: 'copy', pattern: "figs"
    publishDir "${params.outdir}/pipeline_info/", mode: 'copy', pattern: "sessionInfo_deg.txt"
    input:
    file salmon_files
    file index
    file transcriptome
    file gtf
    val organism
    file design
    file contrast_matrix
    val min_size
    val max_size
    output:
    path "figs", emit: figures
    path "csv", emit: csv
    path "sessionInfo_deg.txt", emit: session_info
    script:
    """
    mkdir cache
    differential_genes.R ${transcriptome} ${gtf} ${organism} ${design} ${contrast_matrix} ${min_size} ${max_size}

    """

 }

  process differential_transcripts{
    cpus 2
    memory 16.GB
    tag "Differential Expression"
    label "edgeR"
    publishDir "${params.outdir}/differential_transcripts", mode: 'copy', pattern: "csv"
    publishDir "${params.outdir}/differential_transcripts/", mode: 'copy', pattern: "figs"
    publishDir "${params.outdir}/pipeline_info/", mode: 'copy', pattern: "sessionInfo_det.txt"
    input:
    file salmon_files
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
    differential_transcripts.R  ${organism} ${design} ${contrast_matrix}

    """

 }
