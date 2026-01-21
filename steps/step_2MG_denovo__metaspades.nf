nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2MG_denovo'
def METHOD = 'metaspades' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_2MG_denovo__metaspades{
    take: data
    main:
      spades_out = metaspades(data)
      filter_out = assembly_filter(spades_out.scaffolds)
      quast(filter_out.fasta)
    emit:
      assembled = filter_out.fasta
}

workflow {
    step_2MG_denovo__metaspades(getSingleInput())
}