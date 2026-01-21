nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata;taskMemory;extractKey } from '../functions/common.nf'
include { getInput;getKingdom;getReferenceOptional;checkEnum } from '../functions/parameters.nf'
include { stepInputs } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4AN_genes'
def METHOD = 'prokka'
def ENTRYPOINT = "step_${STEP}__${METHOD}"

enum KINGDOM {
  Viruses, Bacteria, Archaea, Mitochondria
}

workflow prepare_inputs {
    main:
      getInput()
        .cross(getKingdom()) { extractKey(it) }
        .cross(getReferenceOptional('gb')) { extractKey(it) }
        .map { it.flatten() }  // [ riscd assembly ds kingdom ds riscd_ref refid refpath]
        .map { 
            [ it[0], it[1], it[3], it[5], it[6], it[7] ]
        }
        .set { input }
    emit:
      input
}

workflow step_4AN_genes__prokka {
    take: 
      data
    main:
      prokka_out = prokka(data)
}

workflow {
  input = prepare_inputs()
  step_4AN_genes__prokka(input)
}