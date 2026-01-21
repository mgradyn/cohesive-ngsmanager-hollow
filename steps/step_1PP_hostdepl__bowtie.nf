nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common'
include { param;getSingleInput;getHostUnkeyed;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def refDir = param('hosts_dir')
def ex = executionMetadata()
def STEP = '1PP_hostdepl'
def METHOD = 'bowtie' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow prepare_inputs {
    main:
      ch_input = getSingleInput()
      ch_host = getHostUnkeyed()
      
      ch_combined = ch_input.combine(ch_host)
    emit:
      combined_inputs = ch_combined
}

workflow step_1PP_hostdepl__bowtie {
    take: 
      trimmedAndHost
    main:
      bowtie_out = bowtie2(trimmedAndHost)
      samtools_out = samtools(bowtie_out.sam)
    emit:
      samtools_out.depleted      
}

workflow {    
    inputs = prepare_inputs()
    step_1PP_hostdepl__bowtie(inputs.combined_inputs)    
}