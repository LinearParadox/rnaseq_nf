
 /*
Zipping outputs for easy download!
 */

 process zip_outputs{
    cpus 1
    memory 8.GB
    publishDir "${params.outdir}/", mode: 'copy'
    input:
    path differential_transcripts, name: 'differential_transcripts/csv'
    path differential_genes, name: 'differential_genes/csv'
    output:
    path "out_files.zip", emit: zip
    script:
    """
    zip -r out_files.zip ./*
    """
 }