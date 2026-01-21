nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;stepInputs;extractKey } from '../functions/common.nf'
include { getSingleInput;isCompatibleWithSeqType;getReference;getRisCd} from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '1PP_filtering'
def METHOD = 'minimap2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow prepare_inputs {
    main:
      getSingleInput()
        .cross(getReference('fa')) { extractKey(it) }
        .multiMap { 
            reads: it[0]
            refs:  it[1][1..3]
        }
        .set { prepared_data }
    emit:
      reads = prepared_data.reads
      refs = prepared_data.refs
}

workflow step_1PP_filtering__minimap2 {
    take: 
      reads
      reference
    main:
      minimap_out = minimap2(reads, reference)
      samtools_out = samtools(minimap_out.sam)
    emit:
      samtools_out.filtered  
}

workflow {
    inputs = prepare_inputs()
    step_1PP_filtering__minimap2(inputs.reads, inputs.refs)
}