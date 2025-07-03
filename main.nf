#!/usr/bin/env nextflow

/*  
 * Basic nextflow pipeline for rna seq
 */

include { qc_samples } from './workflows/qc_workflow.nf'
include { salmon_quant } from './modules/salmon/salmon.nf'
include { salmon_index } from './modules/salmon/salmon.nf'
include { star } from './workflows/STAR.nf'
include { differential_expression } from './modules/edgeR/main.nf'
include { differential_transcripts } from './modules/edgeR/main.nf'
include { multiqc } from './modules/multiqc/multiqc.nf'
nextflow.enable.moduleBinaries = true

workflow {
    if (!params.gtf | !params.samplesheet){
        error "A gtf file and a samplesheet must be provided for the pipeline to run."
    }
    samples=Channel.fromPath(params.samplesheet).splitCsv().map {
        def sample = it[0]
        def r1 = file(it[1])
        def r2 = file(it[2])
        return [sample, r1, r2]
    } | groupTuple()
    qc_samples(samples)
    star_logs = channel.empty()
    if (params.star_align) {
        star(qc_samples.out.trimmed, params.star_index, params.gtf, params.genome, params.readlength)
        star_logs = star.out.starlog.collect()
    }
    if (!file(params.salmon_index).exists()) {
        salmon_index = salmon_index(file(params.salmon_transcriptome), file(params.genome), params.kmer_size).index.collect()
    } else {
        salmon_index = channel.fromPath(params.salmon_index).collect()
    }
    salmon_quant = salmon_quant(qc_samples.out.trimmed, salmon_index, file(params.gtf), params.library_type, params.gibbs_sampling,
                                params.seq_bias, params.gc_bias, params.pos_bias, params.dump_eq)
    salmon_files = salmon_quant.salmon_file.collect()
    differential_expression(
        salmon_files,
        salmon_index,
        file(params.salmon_transcriptome),
        file(params.gtf),
        params.organism,
        file(params.design),
        file(params.contrast_matrix),
        params.min_gs_size,
        params.max_gs_size
    )
    differential_transcripts(
        salmon_files,
        params.organism,
        file(params.design),
        file(params.contrast_matrix)
    )
    multiqc(
        qc_samples.out.multiqc.collect(),
        salmon_quant.salmon_file.collect(),
        star_logs.ifEmpty([])
    )
}
