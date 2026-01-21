nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getInput;param;optionalOrDefault;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

SPECIES_SCHEMA = [
  listeria_monocytogenes : ['l_mono_chewie_1748_220623'],
  escherichia_coli : ['e_coli_chewie_2360_210531'],
  salmonella_enterica : ['s_enterica_chewie_3255_210531']
]

SCHEMAS = [
  l_mono_chewie_1748_220623 : "/schemas/Listeria_monocytogenes_Pasteur_cgMLST_2022-06-23T18_03_54.613576.zip",
  e_coli_chewie_2360_210531 : "/schemas/Escherichia_coli_INNUENDO_wgMLST_2021-05-31T14_24_05.304225.zip",
  s_enterica_chewie_3255_210531 : "/schemas/Salmonella_enterica_INNUENDO_cgMLST_2021-05-31T20_28_21.350919.zip"
]

CHEWBBACA_SINGLE_END_PARAMS = [
  'l_mono_chewie_1748_220623': ' --minimum-length 144 --st 0.1 --bsr 0.6 ',
  'e_coli_chewie_2360_210531': ' --minimum-length 0 --st 0.01 --bsr 0.6 --genes-list /schemas/Escherichia_coli_INNUENDO_cgMLST_EFSA_filterlist.txt ',
  's_enterica_chewie_3255_210531': ' --minimum-length 0 --st 0.01 --bsr 0.6 --genes-list /schemas/Salmonella_enterica_INNUENDO_cgMLST_EFSA_filterlist.txt '
]

CHEWBBACA_PAIRED_END_PARAMS = [
  'l_mono_chewie_1748_220623': ' --minimum-length 144 ',
  'e_coli_chewie_2360_210531': ' --minimum-length 0 --genes-list /schemas/Escherichia_coli_INNUENDO_cgMLST_EFSA_filterlist.txt ',
  's_enterica_chewie_3255_210531': ' --minimum-length 0 --genes-list /schemas/Salmonella_enterica_INNUENDO_cgMLST_EFSA_filterlist.txt '
]

def getExtraParams(schema) {
  try {   
    if (isIonTorrent(null)) {
      if (CHEWBBACA_SINGLE_END_PARAMS.containsKey(schema)) return CHEWBBACA_SINGLE_END_PARAMS[schema]
    } else {
      if (CHEWBBACA_PAIRED_END_PARAMS.containsKey(schema)) return CHEWBBACA_PAIRED_END_PARAMS[schema]
    }        
    return ''
  } catch(Throwable t) { exit 1, "Error: ${t.asString()}" } 
}

def ex = executionMetadata()
def STEP = '4TY_cgMLST'
def METHOD = 'chewbbaca' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getSchema(gsp, schema) {
  try {  
    def genus_species = gsp ? gsp.toLowerCase() : ''
    def (genus, species) = genus_species.contains("_") ? genus_species.split('_') : [ genus_species, null ]
    def allowedSchemas = []
    if (SPECIES_SCHEMA.containsKey(genus_species)) {
      allowedSchemas =  SPECIES_SCHEMA.get(genus_species) 
    } else if (SPECIES_SCHEMA.containsKey(genus)) {
      allowedSchemas = SPECIES_SCHEMA.get(genus)
    }
    if (!allowedSchemas) return null
    if (!schema) return allowedSchemas[0] 
    if (allowedSchemas.contains(schema)) return schema
    return null;      
  } catch(Throwable t) { exit 1, "Error: ${t.asString()}" } 
}

workflow step_4TY_cgMLST__chewbbaca {
    take: 
      assembly
      genus_species
      schema
    main:
      chewbbaca_out = chewbbaca(assembly, genus_species, schema)      
      hashing(chewbbaca_out.alleles)
      chewbbaca_check(chewbbaca_out.stats)
}

workflow {
    step_4TY_cgMLST__chewbbaca(getInput(), param('genus_species'), optionalOrDefault('schema', ''))
}