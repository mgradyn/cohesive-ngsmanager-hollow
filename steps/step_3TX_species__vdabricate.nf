nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata; taskMemory } from '../functions/common.nf'
include { param;getInput } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def abricateDBDir = param('step_3TX_species__vdabricate__db')
def ex = executionMetadata()
def STEP = '3TX_species'
def METHOD = 'vdabricate' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow prepare_inputs {
    main:
      ch_combined = getInput().combine(['viruses_TREF'])
    emit:
      scaffoldsAbricatedb = ch_combined
}

workflow step_3TX_species__vdabricate {
    take: data
    main:
      abricate_out = abricate(data)
    emit:
      calls = abricate_out.calls
}

workflow {
    inputs = prepare_inputs()
    step_3TX_species__vdabricate(inputs.scaffoldsAbricatedb)
}