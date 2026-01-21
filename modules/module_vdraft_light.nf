nextflow.enable.dsl=2

include { step_1PP_hostdepl__bowtie } from '../steps/step_1PP_hostdepl__bowtie'
include { step_2AS_mapping__bowtie } from '../steps/step_2AS_mapping__bowtie'
include { extractKey } from '../functions/common.nf'
include { getSingleInput;getHostOptional;getReference } from '../functions/parameters.nf'

workflow module_vdraft_light {
    take: 
        trimmedReads
        host
        reference
    main:
        trimmedReads.cross(host) { extractKey(it) }
            .map { [ it[0][0], it[0][1], it[1][1] ] }
            .branch {
                with_host: it[1][1]
                without_host: true
            }
            .set { ch_branched }

        ch_depleted = step_1PP_hostdepl__bowtie(ch_branched.with_host)

        ch_branched.without_host
            .mix(ch_depleted)
            .map { it[0,1] }
            .set { ch_ready }

        ch_ready.cross(reference) { extractKey(it) }
            .multiMap { 
                reads: it[0]
                refs:  it[1][1..3]
            }
            .set { ch_final }

        step_2AS_mapping__bowtie(ch_final.reads, ch_final.refs)
}

workflow {
    ch_in = getSingleInput()
    ch_host = getHostOptional()
    ch_ref = getReference('fa')
    module_vdraft_light(ch_in, ch_host, ch_ref)
}