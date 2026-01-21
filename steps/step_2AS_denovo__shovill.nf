nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2AS_denovo'
def METHOD = 'shovill' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_2AS_denovo__shovill {
    take: rawreads
    main:
      shovill_out = shovill(rawreads)
      shovill_se_out = shovill_se(rawreads)
      
      contigs = shovill_out.assembly.mix(shovill_se_out.assembly)
      
      quast(contigs)
      
      if (!params.skip_checkm) {
        checkm(contigs)
      }
    emit:
      assembly = contigs
}

workflow {
    step_2AS_denovo__shovill(getSingleInput())
}