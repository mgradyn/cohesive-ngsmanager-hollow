nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'
include { param } from '../functions/parameters.nf'

FILTERABLE_REFERENCES_PATH = param('step_2AS_filtering__seqio')

def ex = executionMetadata()
def STEP = '2AS_denovo'
def METHOD = 'spades' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def isReferenceFilterable(refCode, _path) {
    try {
      def referencePath = (_path instanceof java.util.Collection) ? _path.flatten()[0] : _path
      if (referencePath.empty) return false
      return referencePath.toRealPath().toString().contains(FILTERABLE_REFERENCES_PATH)
    } catch(Throwable t) { return false } 
}

workflow step_2AS_filtering__seqio {
    take: 
      calls
      assembly
      reference
    main:
      filter(calls, assembly, reference)
}