nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;isSpeciesSupported } from '../functions/common.nf'
include { getInput;param } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4TY_flaA'
def METHOD = 'flaA' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def MLST_SCHEMA_NAME = 'flaA'
def GENUS_ALLOWED = [
  'campylobacter'
]

workflow step_4TY_flaA__flaA {
    take: 
      assembly
      genus_species
    main:
      mlst_out = mlst_flaa(assembly, genus_species)
}

workflow {
    step_4TY_flaA__flaA(getInput(), param('genus_species'))
}