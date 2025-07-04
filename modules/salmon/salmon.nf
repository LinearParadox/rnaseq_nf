 #!/usr/bin/env nextflow

process salmon_quant{
    label 'salmon'
    tag "Salmon quant"
    cpus 16
    memory 24.GB
    publishDir "${params.outputDir}/per-sample-outs/${sample}/", mode: 'copy', pattern: "${sample}", saveAs: {file -> "salmon_quant"}
    publishDir "${params.outputDir}/pipeline_info/", mode: "copy", pattern: "salmon_version.txt"
    input:
    tuple val(sample), path(r1), path(r2)
    path salmon_index
    file gtf
    val library_type
    val gibbs_sampling
    val seq_bias
    val gc_bias
    val pos_bias
    val dump_eq
    output:
    tuple val(sample), path("${sample}"), emit: salmon_output
    path ("${sample}"), emit: salmon_file
    path "salmon_version.txt", emit: salmon_version



    script:
    """
    salmon -v > salmon_version.txt
    command="salmon quant -l A \
        -i ${salmon_index} \
        -1 ${r1} \
        -2 ${r2} \
        --validateMappings \
        -p ${task.cpus} \
        --numGibbsSamples ${gibbs_sampling} \
        -g ${gtf} \
        -o ${sample}"
    if [ ${seq_bias} = true ]; then
        command+=" --seqBias"
    fi
    if [ ${gc_bias} = true ]; then
        command+=" --gcBias"
    fi
    if [ ${dump_eq} = true ]; then
        command+=" --dumpEq"
    fi
    if [ ${pos_bias} = true ]; then
        command+=" --posBias"
    fi
    eval "\${command}"
    """
}

process salmon_index{
    label 'salmon'
    tag "Salmon index"
    publishDir "${params.outputDir}/ref/", mode: 'copy'
    input:
    file transcriptome
    file genome
    val kmer_size
    output:
    path "salmon_index", emit: index
    script:
    """
    grep "^>" <(gunzip -cf ${genome}) | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat ${transcriptome} ${genome} > salmon_genome.fa.gz
    salmon index -t salmon_genome.fa.gz -d decoys.txt -p ${task.cpus} -i salmon_index --gencode
    """
} 
