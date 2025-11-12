process STARfusion {
    container 'trinityctat/starfusion:1.15.1'
    cpus 16
    memory 96.GB
    label 'fusion'

    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "fusion"

    input:
    tuple val(sample), path(r1), path(r2)
    path index 
    output:
    tuple val(sample), path("fusion"), emit: bam

    script:
    """
    mkdir fusion
    STAR-Fusion --genome_lib_dir ${index} \\
        --left_fq ${r1} \\
        --right_fq ${r2} \\
        --output_dir fusion \\
        --FusionInspector 'validate' \\
        --examine_coding_effect \\
        --denovo_reconstruct \\
        --extract_fusion_reads \\
        --CPU ${task.cpus}
    """
}