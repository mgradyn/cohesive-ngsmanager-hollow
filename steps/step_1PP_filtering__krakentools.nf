nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;extractKey } from '../functions/common'
include { getTrimmedReads;getParamTaxaId;getParamIncludeChildren;getParamIncludeParents;getKrakenResults } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '1PP_filtering'
def METHOD = 'krakentools' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow prepare_inputs {
    main:
      getKrakenResults().cross(getTrimmedReads(false)) { extractKey(it) }.multiMap { 
        kraken: it[0]
        trimmed: it[1]
      }.set { krakenAndTrimmed }
    
    emit:
      kraken = krakenAndTrimmed.kraken
      trimmed = krakenAndTrimmed.trimmed
}

workflow step_1PP_filtering__krakentools {
    take: 
      kraken
      trimmed
      taxaid
      include_children
      include_parents
    main:
      krakentools(kraken, trimmed, taxaid, include_children, include_parents)
}

workflow {
    inputs = prepare_inputs()
    step_1PP_filtering__krakentools(inputs.kraken, inputs.trimmed, getParamTaxaId(), getParamIncludeChildren(), getParamIncludeParents())
}