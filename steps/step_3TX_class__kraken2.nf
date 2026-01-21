nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { isSarsCov2;isPositiveControlSarsCov2;isNegativeControlSarsCov2 } from '../functions/sampletypes'
include { param;getSingleInput;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters'
include { stepInputs;getRisCd } from '../functions/common.nf'

def KRAKEN2_DB = param('step_3TX_class__kraken2__db')

def ex = executionMetadata()
def STEP = '3TX_class'
def METHOD = 'kraken2' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_3TX_class__kraken2 {
    take: reads
    main:
      kraken_out = kraken2(reads)
      braken_out = braken2(kraken_out.report)
     emit:
       genus_report = braken_out.genus_report
}

workflow {
  step_3TX_class__kraken2(getSingleInput())
}