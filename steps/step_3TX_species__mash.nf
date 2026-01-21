nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '3TX_species'
def METHOD = 'mash' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"


workflow step_3TX_species__mash {
    take: 
      reads
    main:
      mash_out = mash(reads)
}

workflow {
    step_3TX_species__mash(getInput())
}