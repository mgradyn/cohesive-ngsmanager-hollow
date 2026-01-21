nextflow.enable.dsl=2

include { getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { isSarsCov2;isPositiveControlSarsCov2;isNegativeControlSarsCov2;isNGSMG16S } from '../functions/sampletypes'
include { param;isFullOutput;getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters'
include { stepInputs;getRisCd } from '../functions/common.nf'

def db_kraken=param('step_3TX_class__kraken__db_kraken')
def db_bracken=param('step_3TX_class__kraken__db_bracken')
def ex = executionMetadata()
def STEP = '3TX_class'
def METHOD = 'kraken' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_3TX_class__kraken {
    take: reads
    main:
      kraken_out = kraken(reads) 
      braken_out = braken(kraken_out.report)
     emit:
       genus_report = braken_out.genus_report
}

workflow {
  step_3TX_class__kraken(getSingleInput())
}