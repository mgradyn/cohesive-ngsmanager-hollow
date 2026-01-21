nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getInput;param } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4TY_MLST'
def METHOD = 'mlst' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_4TY_MLST__mlst {
    take: 
      assembly
    main:
      mlst_out = mlst(assembly)
}

workflow {
  step_4TY_MLST__mlst(getInput())
}