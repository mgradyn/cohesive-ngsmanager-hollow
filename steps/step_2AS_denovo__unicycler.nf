nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2AS_denovo'
def METHOD = 'unicycler' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_2AS_denovo__unicycler {
    take: data
    main:
      unicycler_out = unicycler(data)
      filter_out = assembly_filter(unicycler_out.scaffolds)
      quast(filter_out.fasta)
    emit:
      assembled = filter_out.fasta       
}

workflow {
    step_2AS_denovo__unicycler(getSingleInput())
}