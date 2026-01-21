nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getSingleInput;param } from '../functions/parameters.nf'
include { stepInputs;flattenPath } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4AN_AMR'
def METHOD = 'staramr' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_ALLOWED = [
  'campylobacter'
]

workflow step_4AN_AMR__staramr {
    take: 
      assembly
      genus_species
    main:
      staramr_out = staramr(assembly, genus_species)
}

workflow {
    step_4AN_AMR__staramr(getSingleInput(), param('genus_species'))
}