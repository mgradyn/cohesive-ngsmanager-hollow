nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4AN_AMR'
def METHOD = 'resfinder'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_4AN_AMR__resfinder {
    take: 
      reads
      genus_species
    main:
      resfinder_out = resfinder(reads, genus_species)
}

workflow {
    step_4AN_AMR__resfinder(getSingleInput(),param('genus_species'))
}