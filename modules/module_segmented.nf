nextflow.enable.dsl=2

include { extractKey } from '../functions/common.nf'
include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'
include { getSingleInput;getReferences } from '../functions/parameters.nf'

workflow module_segmented {
    take: 
        reads
        reference 
    main:
        ivar_results = step_2AS_mapping__ivar(reads, reference)
    emit:
        consensus = ivar_results.consensus
}

workflow prepare_inputs {
    take:
        raw_reads
        raw_refs
    main:
        raw_reads.cross(raw_refs) { extractKey(it) }
            .multiMap { 
                reads: it[0] 
                refs:  it[1][1..3] 
            }
            .set { ch_prepared }
    emit:
        reads = ch_prepared.reads
        refs = ch_prepared.refs
}

workflow {
    ch_in = getSingleInput()
    ch_ref = getReferences('any')
    
    ch_ready = prepare_inputs(ch_in, ch_ref)
    
    module_segmented(ch_ready.reads, ch_ready.refs)
}