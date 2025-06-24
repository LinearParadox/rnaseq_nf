include { fastP } from '../modules/fastp/qc.nf'
include { merge_lanes } from '../modules/fastp/qc.nf'

workflow qc_samples {
    take:
        samples
        fastp_threads
    main:
        fastp_results = fastP(merge_lanes(samples), fastp_threads)
    emit:
        trimmed = fastp_results.reads
        multiqc = fastp_results.fastp_results
    }