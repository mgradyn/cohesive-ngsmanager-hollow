nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;extractKey;taskMemory;stepInputs;getRisCd;extractDsRef } from '../functions/common.nf'
include { getSingleInput;getReference;isCompatibleWithSeqType } from '../functions/parameters.nf'

def ex = executionMetadata()
def STEP = '2AS_mapping'
def METHOD = 'medaka' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

workflow prepare_inputs {
    main:
      getSingleInput()
        .cross(getReference('fa')) { extractKey(it) }
        .multiMap { 
            reads: it[0]
            refs:  it[1][1..3]
        }
        .set { prepared_data }
    emit:
      reads = prepared_data.reads
      refs = prepared_data.refs
}

workflow step_2AS_mapping__medaka {
  take: 
    reads
    reference 
  main:
    medaka_out = medaka(reads, reference)
    minmax_out = coverage_minmax(medaka_out.bam, METHOD)
    coverage_plot(minmax_out.coverage_depth)
    depth_out = samtools_depth(medaka_out.bam, METHOD)
    ch_cov_keyed = depth_out.coverage.map { [extractDsRef(it), it] }
    ch_con_keyed = medaka_out.consensus.map { [extractDsRef(it), it] }
    ch_cov_con_crossed = ch_cov_keyed.cross(ch_con_keyed)
    ch_cov_con_crossed.map { 
        def cov = it[0][1] 
        def con = it[1][1] 
        return [ cov[0], con[1], cov[1] ]
    }.set { coverageRefAndConsensus }
    
    check_out = coverage_check(coverageRefAndConsensus, METHOD)
    coverageBasic = check_out.coverage_basic
    ch_extra_keyed = minmax_out.coverage_extra.map { [it[0] + "-" + it[1], it] }
    ch_basic_keyed = coverageBasic.map { [it[0] + "-" + it[1], it] }
    ch_checks_crossed = ch_extra_keyed.cross(ch_basic_keyed)
    ch_checks_crossed.map { 
        def extra = it[0][1]
        def basic = it[1][1]
        return [ extra[0], extra[1], extra[2], basic[2] ]
    }.set { crossedChecks }
    coverage_check_merge(crossedChecks, METHOD)
  emit:
    consensus = medaka_out.consensus
}

workflow {
    inputs = prepare_inputs()
    step_2AS_mapping__medaka(inputs.reads, inputs.refs)
}