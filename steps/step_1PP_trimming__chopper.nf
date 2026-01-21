nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getInput;isCompatibleWithSeqType;paramWrap;optWrap } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '1PP_trimming'
def METHOD = 'chopper' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_1PP_trimming__chopper {
    take: 
      rawreads
    main:
      chopper_out = chopper(rawreads)
      nanoplot(chopper_out.trimmed)
    emit:
      chopper_out.trimmed      
}

workflow {
    step_1PP_trimming__chopper(getInput())
}