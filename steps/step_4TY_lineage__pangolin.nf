nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4TY_lineage'
def METHOD = 'pangolin' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_4TY_lineage__pangolin {
    take: 
      consensus
    main:
      pangolin_out = pangolin(consensus)
}

workflow {
    step_4TY_lineage__pangolin(getInput())
}