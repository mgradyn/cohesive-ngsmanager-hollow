nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory;stepInputs } from '../functions/common.nf'
include { getSingleInput;getReference;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '2AS_mapping'
def METHOD = 'snippy' 
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

workflow step_2AS_mapping__snippy {
  take: 
    reads
    reference 
  main:
    snippy(reads, reference)
}

workflow {
    inputs = prepare_inputs()
    step_2AS_mapping__snippy(inputs.reads, inputs.refs)
}