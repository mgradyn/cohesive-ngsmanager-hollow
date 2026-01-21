nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4AN_AMR'
def METHOD = 'abricate'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_4AN_AMR__abricate {
    take: data
    main:
      abricate_out = abricate(data)
}

workflow {
    step_4AN_AMR__abricate(getInput())
}