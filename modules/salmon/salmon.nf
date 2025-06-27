 #!/usr/bin/env nextflow

process salmon_quant{
    label 'salmon'
    tag "Salmon quant"
    cpus 4
    memory 12.GB
    publishDir "${params.outdir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "*.sf"
    publishDir "${params.outdir}/per-sample-outs/${sample}/equiv_classes/", mode: 'copy', pattern: "equiv_classes/*", when: params.dump_eq
    input:
    tuple val(sample), path(bam)
    file transcriptome
    file gtf
    val gibbs_sampling
    val seq_bias
    val gc_bias
    val dump_eq
    output:
    path "*.sf", emit: quant

    script:
    """
    command="salmon quant -l A \
        -a ${bam} \
        -t ${transcriptome} \
        -p ${task.cpus} \
        --numGibbsSamples ${gibbs_sampling} \
        -g ${gtf} \
        --gencode \
        -o salmon_quant_${sample}"
    if [ ${seq_bias} = true ]; then
        command+=" --seqBias"
    fi
    if [ ${gc_bias} = true ]; then
        command+=" --gcBias"
    fi
    if [ ${dump_eq} = true ]; then
        command+=" --dumpEq --auxDir equiv_classes"
    fi
    eval "\${command}"
    """
}