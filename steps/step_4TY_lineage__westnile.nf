nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory } from '../functions/common.nf'
include { getSingleInput;_getSingleReference } from '../functions/parameters.nf'
include { stepInputs;csv2map;getRisCd;parseRISCD } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '4TY_lineage'
def METHOD = 'westnile' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

STEP_ASSETS_DIR = "${params.assets_dir}/${ENTRYPOINT}"
def HCOV_THRESHOLD = params.step_4TY_lineage__westnile___threshold
def WESTNILE_LINEAGE_REFERENCES_PATH = "${STEP_ASSETS_DIR}/westnile_lineage_references.json"
WESTNILE_LINEAGE_REFERENCES = new groovy.json.JsonSlurper().parseText(file(WESTNILE_LINEAGE_REFERENCES_PATH).text)

def getReferenceForLineage(lineageFile) {
   try {        
      def lineage = csv2map(lineageFile, ",").lineage
      def refData = WESTNILE_LINEAGE_REFERENCES.find { it.lineage == lineage }
      return [ refData.ref_riscd, refData.ref_code, "${STEP_ASSETS_DIR}/${refData.ref_path}" ]
   } catch(Throwable t) { exit 1, "Error: ${t.asString()}" }   
}

def isValidLineage(lineageResult) {
   try {        
      def lineage = csv2map(lineageResult[1], ",").lineage   
      def notAssigned =  lineage == 'NA'
      def notDetermined = lineage == 'ND'
      return !notAssigned && !notDetermined
   } catch(Throwable t) { exit 1, "Error: ${t.asString()}" }   
}
 
def getWNVReferences() {
   try {        
      return Channel.fromList(WESTNILE_LINEAGE_REFERENCES.collect { [ it.ref_riscd, it.ref_code, "${STEP_ASSETS_DIR}/${it.ref_path}" ] } )
   } catch(Throwable t) { exit 1, "Error: ${t.asString()}" }   
}

workflow step_4TY_lineage__westnile {
    take: 
      reads
    main:
      references = getWNVReferences()
      comb = reads.combine(references)
      
      snippy_out = snippy(comb)
      samtools_out = samtools_depth(snippy_out.bam)
      
      grouped_cov = samtools_out.coverage.groupTuple()
      
      lineage_out = lineage_selection(grouped_cov)
    emit:
      lineage = lineage_out.lineage.filter { isValidLineage(it) }
}

workflow {  
    step_4TY_lineage__westnile(getSingleInput())
}