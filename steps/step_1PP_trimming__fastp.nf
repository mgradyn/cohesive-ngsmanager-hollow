nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getInput;isIonTorrent;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd;extractKey } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '1PP_trimming'
def METHOD = 'fastp' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_1PP_trimming__fastp {
    take: 
      rawreads
    main:
      trimmed = fastp(rawreads).trimmed

      fastqc(trimmed)
      readsCheckInput = rawreads.cross(trimmed) { extractKey(it) }.multiMap { 
        rawreads: it[0]
        trimmed: it[1]
      }      
      sample_reads_check(readsCheckInput.rawreads, readsCheckInput.trimmed)

    emit:
      trimmed      
}

workflow {
    step_1PP_trimming__fastp(getInput())
}