nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata } from '../functions/common.nf'
include { getSingleInput;param;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '1PP_downsampling'
def METHOD = 'bbnorm' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_1PP_downsampling__bbnorm {
    take: 
      reads
      k
      target
    main:
      bbnorm(reads, k, target)
}

workflow {
    step_1PP_downsampling__bbnorm(getSingleInput(), param('k'), param('target'))
}