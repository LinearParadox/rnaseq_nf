
 /*
Zipping outputs for easy download!
 */

 process zip_outputs{
    cpus 1
    memory 8.GB
    publishDir "${params.outdir}/outs.zip", mode: 'copy', pattern: "zip"
    input:
    path differential_transcripts, name: 'differential_transcripts/csv'
    path differential_genes, name: 'differential_genes/csv'
    output:
    file "outs.zip"
    script:
    """
    zip -r outs.zip ./*
    """
 }