nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;param } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '4AN_AMR'
def METHOD = 'filtering'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_4AN_AMR__filtering {
    take: 
      data
      coverage
      identity
    main:
      filtering_out = abricate_filtering(data, coverage, identity)
}

workflow {
    step_4AN_AMR__filtering(getSingleInput(), param('coverage'), param('identity'))
}