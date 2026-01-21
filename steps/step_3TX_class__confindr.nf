nextflow.enable.dsl=2

include { stepInputs;parseMetadataFromFileName;executionMetadata;isSpeciesSupported;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '3TX_class'
def METHOD = 'confindr' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def GENUS_ALLOWED = [
  'escherichia',
  'salmonella',
  'listeria' 
]

workflow step_3TX_class__confindr {
    take: 
      reads
      genus_species
    main:
      confindr(reads, genus_species)
}

workflow {
    step_3TX_class__confindr(getSingleInput(),param('genus_species'))
}