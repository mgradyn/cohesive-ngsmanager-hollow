nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;getLongReads;param } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2AS_hybrid'
def METHOD = 'unicycler' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_2AS_hybrid__unicycler {
    take: 
      short_reads
      long_reads
    main:
      unicycler_out = unicycler(short_reads, long_reads)
      quast(unicycler_out.scaffolds)
    emit:
      scaffolds = unicycler_out.scaffolds
}

workflow {  
    step_2AS_hybrid__unicycler(getSingleInput(), getLongReads())
}