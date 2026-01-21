nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory;extractKey } from '../functions/common'
include { stepInputs;getRisCd } from '../functions/common.nf'
include { getSingleInput;getReference;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '1PP_filtering'
def METHOD = 'bowtie' 

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

workflow step_1PP_filtering__bowtie {
    take: 
      reads
      reference 
    main:
      bowtie2(reads, reference) 
      samtools(bowtie2.out.sam)
    emit:
      samtools.out.filtered      
}

workflow {
    data = prepare_inputs()
    step_1PP_filtering__bowtie(data.reads, data.refs)
}