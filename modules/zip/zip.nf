
 /*
Zipping outputs for easy download!
 */

 process zip_outputs{
    cpus 1
    memory 8.GB
    publishDir "${params.outdir}/outs.zip", mode: 'copy', pattern: "zip"
    input:
    tuple val(sample), path(r1), path(r2)
    output:
    path differential_transcripts, name: 'differential_transcripts/csv'
    path differential_genes, name: 'differential_genes/csv'
    script:
    """
    zip -r outs.zip .
    """
 }