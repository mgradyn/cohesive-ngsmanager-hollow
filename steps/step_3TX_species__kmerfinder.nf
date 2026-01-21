nextflow.enable.dsl=2

include { flattenPath; parseMetadataFromFileName; executionMetadata; csv2map; extractKey;taskMemory } from '../functions/common.nf'
include { param;getInput } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

KMERFINDER_SPECIES_DIR = param('step_3TX_species__kmerfinder__db')
KMERFINDER_REFERENCE_DIR = "${KMERFINDER_SPECIES_DIR}/Bacteria/Fasta/"

def ex = executionMetadata()
def STEP = '3TX_species'
def METHOD = 'kmerfinder' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getBacterialReferencePath(checkFile) {
    try {
      def bacterialReferenceId = csv2map(checkFile, "\\t").assembly_accBacteria
      return [ '-', bacterialReferenceId, file("${KMERFINDER_REFERENCE_DIR}/${bacterialReferenceId}*.fa") ] 
    } catch(Throwable t) { exit 1, "Error: ${t.asString()}" }
}

def getCalculatedSpecies(checkFile) {
    try {
      return csv2map(checkFile, "\\t").speciesAssigned
    } catch(Throwable t) { exit 1, "Error: ${t.asString()}" }
}

workflow step_3TX_species__kmerfinder {
    take: data
    main:
      kmerfinder_out = kmerfinder(data)
      
      assigned_species = kmerfinder_out.check.map { 
          [ it[0], getCalculatedSpecies(it[1]), getBacterialReferencePath(it[1]) ] 
      }
    emit:
      assigned_species = assigned_species
}

workflow {
    step_3TX_species__kmerfinder(getInput())
}