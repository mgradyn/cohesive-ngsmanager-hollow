nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common'
include { getInput;isCompatibleWithSeqType;isIlluminaPaired } from '../functions/parameters.nf'
include { stepInputs;parseRISCD } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '0SQ_rawreads'
def METHOD = 'fastq'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_0SQ_rawreads__fastq {
    take: data
    main:
      fastqc(data)

      nanoplot_out = nanoplot(data)
      nanopore_reads_check(nanoplot_out.stats)
}

workflow {
  step_0SQ_rawreads__fastq(getInput())
}