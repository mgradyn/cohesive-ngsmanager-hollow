nextflow.enable.dsl=2

include { extractDsRef;getEmpty;flattenPath; parseMetadataFromFileName; executionMetadata; extractKey;taskMemory } from '../functions/common.nf'
include { getSingleInput;getReferences;getReferenceCodes;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent;isSegmentedMapping } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()
def STEP = '2AS_mapping'
def METHOD = 'ivar' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

def getExDt(reference, ex) {
    try {        
      def reflist = getReferenceCodes()
      if (reflist.size() < 2 || !reference || !reflist.contains(reference)) return ex.dt
      return ex.dt + reflist.indexOf(reference)
    } catch(Throwable t) { return ex.dt }
}

def canBeAggregated(actual) {
    try {        
      def expected = getReferenceCodes()
      if (expected.size() < 2) return false
      if (actual.size() != expected.size()) return false
      if (actual.any {!expected.contains(it)}) return false
      return true
    } catch(Throwable t) { return false }
}

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

workflow step_2AS_mapping__ivar {
  take: 
    reads
    reference 
  main:
    snippy_out = snippy(reads, reference)
    pileup_out = samtools_pileup(snippy_out.bam)
    ivar_out = ivar(pileup_out.pileup)
    consensus = ivar_out.consensus
    minmax_out = coverage_minmax(snippy_out.bam, 'vdsnippy')
    coverage_plot(minmax_out.coverage_depth)
    depth_out = samtools_depth(snippy_out.bam, 'vdsnippy')
    ch_cov_keyed = depth_out.coverage.map { [extractDsRef(it), it] }
    ch_con_keyed = consensus.map { [extractDsRef(it), it] }
    ch_cov_con_crossed = ch_cov_keyed.cross(ch_con_keyed)
    ch_cov_con_crossed.map { 
        def cov = it[0][1] 
        def con = it[1][1] 
        return [ cov[0], con[2], cov[1] ]
    }.set { coverageRefAndConsensus }
    check_out = coverage_check(coverageRefAndConsensus, 'ivar')
    coverageBasic = check_out.coverage_basic
    ch_extra_keyed = minmax_out.coverage_extra.map { [it[0] + "-" + it[1], it] }
    ch_basic_keyed = coverageBasic.map { [it[0] + "-" + it[1], it] }
    ch_checks_crossed = ch_extra_keyed.cross(ch_basic_keyed)
    ch_checks_crossed.map { 
        def extra = it[0][1]
        def basic = it[1][1]
        return [ extra[0], extra[1], extra[2], basic[2] ]
    }.set { crossedChecks }
    merge_out = coverage_check_merge(crossedChecks, 'vdsnippy')
    grouped_merge = merge_out.coverage_merged.groupTuple(by: 0)
    coverage_check_group(grouped_merge, 'vdsnippy')
    ch_reads_keyed = reads.map { [extractKey(it), it] }
    ch_con_keyed_for_agg = consensus.map { [extractKey(it), it] }
    ch_con_grouped = ch_con_keyed_for_agg.groupTuple()
    ch_agg_crossed = ch_reads_keyed.cross(ch_con_grouped)
    ch_agg_crossed.map {
         def r = it[0][1]
         def c_group = it[1][1]
         def refs = c_group.collect { it[1] }
         def paths = c_group.collect { it[2] }
         return [ r[0], refs, paths ]
    }.set { consensus_to_aggregate }
    aggregate(consensus_to_aggregate)
  emit:
    consensus = consensus.map { it[0,2] }
    coverage_depth = minmax_out.coverage_depth
}

workflow {
    inputs = prepare_inputs()
    step_2AS_mapping__ivar(inputs.reads, inputs.refs)
}