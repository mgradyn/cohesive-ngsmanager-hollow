nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;stepInputs;extractKey } from '../functions/common.nf'
include { getSingleInput;isCompatibleWithSeqType;getHostReference;getRisCd } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '1PP_hostdepl'
def METHOD = 'minimap2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow prepare_inputs {
    main:
      getSingleInput().cross(getHostReference()) { extractKey(it) }
        .multiMap { 
            reads: it[0]
            host:  it[1]
        }.set { input }
    emit:
      reads = input.reads
      host = input.host
}

workflow step_1PP_hostdepl__minimap2 {
    take: 
      reads
      host
    main:
      minimap_out = minimap2(reads, host)
      samtools_out = samtools(minimap_out.sam)
    emit:
      samtools_out.depleted 
}

workflow {
    inputs = prepare_inputs()
    step_1PP_hostdepl__minimap2(inputs.reads, inputs.host)
}