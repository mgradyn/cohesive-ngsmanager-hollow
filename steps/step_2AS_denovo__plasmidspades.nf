nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2AS_denovo'
def METHOD = 'plasmidspades' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_2AS_denovo__plasmidspades {
    take: data
    main:
      spades_out = plasmid_spades(data)
      filter_out = assembly_filter(spades_out.scaffolds)
      quast(filter_out.fasta)
    emit:
      assembled = filter_out.fasta
}

workflow {
    step_2AS_denovo__plasmidspades(getSingleInput())
}