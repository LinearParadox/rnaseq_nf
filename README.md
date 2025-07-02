A pipeline for analyzing RNA seq data written in nextflow

All parameters can be passed through a config file:

## Mandatory parameters:  
 -  samplesheet -- A samplesheet, in the format sample_name,fq1,fq2 with no header
 - gtf - gtf for salmon and star (Note even if a salmon index is provided this is needed for how DEG is currently implemented)
 - salmon_transcriptome - transcriptome fasta for salmon (Note even if a salmon index is provided this is needed for how DEG is currently implemented)
 - organism - must be one of "human" or "mouse"
 - design_matrix - design matrix for edgeR in csv format. Generated through model.matrix in R, and then writing to a csv. Must have the same rownames as the samplesheet
 - contrasts - matrix of contrasts in csv format. 

## STAR Parameters
 - star_align - whether to run STAR
 - genome - genome fasta file for STAR
 - star_index - path to a precomputed star index. Will skip indexing if present
 - readlength - readlength for STAR index. Leaving at 100 is generally fine. (Note that 1 is subtracted internally)
## Salmon Parameters
 - gibbs_sampling -  number of gibbs samples for transcript level quantification. 200 is usually fine.
 - gc_bias - whether to use the gc bias correction model in salmon
 - seq_bias - whether to use the seq bias correction model in salmon
 - pos_bias - whether to use the position bias correction model in salmon 
 - dump_eq - whether to dump equivalence classes in salmon
 - salmon_index - path to a precomputed salmon index. This will skip indexing
 - library_type - salmon library type. Autodetect by default
 - kmer_size - kmer size for salmon index generation. 
 ## DEG Parameters
 - min_gs_size - min size for gene sets
 - max_gs_size - max size for gene sets