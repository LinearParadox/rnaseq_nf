 
 /*
 Gene level differential expression analysis using edgeR
 */
 
 process differential_expression{
    cpus 2
    memory 8.GB
    tag "Differential Expression"
    publishDir "${params.outdir}/deg", mode: 'copy'
    input:
    path salmon_files
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
    script:
    """
    mkdir cache
    Rscript ${script} ${samplesheet} ${transcriptome} ${gtf} ${organism} ${design} ${contrast_matrix}

    """

 }