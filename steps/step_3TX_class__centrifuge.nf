nextflow.enable.dsl=2

include { getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common.nf'
include { param;isFullOutput;getSingleInput;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters'
include { stepInputs;getRisCd } from '../functions/common.nf'

def DB_PATH=param('step_3TX_class__centrifuge__db_path')
def DB_NAME=param('step_3TX_class__centrifuge__db_name')
def ex = executionMetadata()
def STEP = '3TX_class'
def METHOD = 'centrifuge' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow step_3TX_class__centrifuge {
    take: reads
    main: 
      centrifuge_out = centrifuge(reads)
      import_taxa(centrifuge_out.reports)
}

workflow {
  step_3TX_class__centrifuge(getSingleInput())
}