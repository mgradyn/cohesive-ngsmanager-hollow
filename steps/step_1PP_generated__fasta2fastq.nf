nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;flattenPath } from '../functions/common.nf'
include { getSingleInput } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '1PP_generated'
def METHOD = 'fasta2fastq' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_1PP_generated__fasta2fastq {
    take: 
      reads
    main:
      fasta2fastq(reads)
}

workflow {
    step_1PP_generated__fasta2fastq(getSingleInput())
}