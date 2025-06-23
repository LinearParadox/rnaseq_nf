 #!/usr/bin/env nextflow

/*
 * Star processes for creating an index
 */
 
 process salmonIndex {
    cpus 8
    publishDir "${params.outdir}/ref/", mode: 'copy'

    input: 
    path full_genome
    path transcriptome
    val output_name
    output:
        path '${output_name}', emit: index
    script:
    """
    grep "^>" <(gunzip -c ${full_genome}) | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat ${transcriptome} ${full_genome} > salmon_combined_genome.fa
    salmon index -t salmon_combined_genome.fa \\
                 -d decoys.txt \\
                 -i ${output_name} \\
                 --gencode \\
                 -p ${task.cpus}
    """
 }
