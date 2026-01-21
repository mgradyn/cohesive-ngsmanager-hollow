nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;extractKey } from '../functions/common'
include { getInput;isIonTorrent;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '1PP_trimming'
def METHOD = 'trimmomatic' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_1PP_trimming__trimmomatic {
    take: rawreads
    main:
      trimmed = trimmomatic(rawreads).fastq;
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
    step_1PP_trimming__trimmomatic(getInput())
}