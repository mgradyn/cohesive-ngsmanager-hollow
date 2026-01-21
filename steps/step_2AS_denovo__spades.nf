nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2AS_denovo'
def METHOD = 'spades' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_2AS_denovo__spades {
    take: data
    main:
      denovo_out = denovo(data)
      filter_out = assembly_filter(denovo_out.scaffolds)
      quast(filter_out.fasta)
    emit:
      assembled = filter_out.fasta
}

workflow {
    step_2AS_denovo__spades(getSingleInput())
}