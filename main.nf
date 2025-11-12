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
include { zip_outputs } from './modules/zip/zip.nf'
include { STARfusion } from './modules/STAR-fusion/fusion.nf'

workflow {
    if (!params.gtf | !params.samplesheet){
        error "A gtf file and a samplesheet must be provided for the pipeline to run."
    }
    if (!params.deg_only) {
        samples=channel.fromPath(params.samplesheet).splitCsv().map { fields ->
            def sample = fields[0]
            def r1 = file(fields[1])
            def r2 = file(fields[2])
            return [sample, r1, r2]
        } | groupTuple()
        qc_samples(samples)
        star_logs = channel.empty()
        if (params.run_star) {
            star(qc_samples.out.trimmed, params.star_index, params.gtf, params.genome, params.readlength)
            star_logs = star.out.starlog.collect()
        }
        if (params.run_salmon) {
            if (!file(params.salmon_index).exists()) {
                salmon_index = salmon_index(file(params.salmon_transcriptome), file(params.genome), params.kmer_size).index.collect()
            } else {
                salmon_index = channel.fromPath(params.salmon_index).collect()
            }
            salmon_quant = salmon_quant(qc_samples.out.trimmed, salmon_index, file(params.gtf), params.library_type, params.gibbs_sampling,
                                        params.seq_bias, params.gc_bias, params.pos_bias, params.dump_eq)
            salmon_files = salmon_quant.salmon_file.collect()
        }
        if (params.run_fusion) {
            STARfusion(qc_samples.out.trimmed, file(params.star_fusion_index))
        }
        if (params.run_multiqc) {
        multiqc(
            qc_samples.out.multiqc.collect(),
            salmon_files.ifEmpty([]),
            star_logs.ifEmpty([])
        )
        }
    } else {
        salmon_files = channel.fromPath(params.samplesheet).splitCsv().map { fields ->
            def sample = fields[0]
            def quantFile = file(fields[1])
            return [sample, quantFile]
        }.collect()
    }
    if (params.run_deg) {
    differential_expression(
        salmon_files,
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
    zip_outputs(differential_transcripts.out.csv, differential_expression.out.csv)
    }
}
