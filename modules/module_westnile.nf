nextflow.enable.dsl=2

include { step_4TY_lineage__westnile; getReferenceForLineage } from '../steps/step_4TY_lineage__westnile'
include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'
include { extractKey } from '../functions/common.nf'
include { getSingleInput } from '../functions/parameters.nf'

workflow module_westnile {
    take: 
        reads        
    main:
        ch_lineage = step_4TY_lineage__westnile(reads)

        reads.cross(ch_lineage) { extractKey(it) }
            .multiMap { 
                reads: it[0]
                reference: getReferenceForLineage(it[1][1])
            }
            .set { ch_ready }
        
        ivar_out = step_2AS_mapping__ivar(ch_ready.reads, ch_ready.reference)

    emit:
        consensus = ivar_out.consensus
}

workflow {
    ch_input = getSingleInput()
    module_westnile(ch_input)
}