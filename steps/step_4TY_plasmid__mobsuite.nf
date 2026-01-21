nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;param;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4TY_plasmid'
def METHOD = 'mobsuite'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_4TY_plasmid__mobsuite {
    take: 
      reads
    main:
      mobsuite_out = mobsuite(reads)
    emit:
      plasmids = mobsuite_out.plasmids
}

workflow {
    step_4TY_plasmid__mobsuite(getSingleInput())
}