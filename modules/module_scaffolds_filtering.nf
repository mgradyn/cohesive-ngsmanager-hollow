nextflow.enable.dsl=2

include { step_2AS_filtering__seqio } from '../steps/step_2AS_filtering__seqio'
include { step_3TX_species__vdabricate } from '../steps/step_3TX_species__vdabricate'
include { extractKey } from '../functions/common.nf'
include { getInput;getReference;getDS } from '../functions/parameters.nf'

workflow module_scaffolds_filtering {    
    take: 
        assembled
        reference 
        abricatedatabase 
    main:
        assembled.cross(abricatedatabase) { extractKey(it) }
            .map { 
                [ it[0][0], it[0][1], it[1][1] ] 
            }
            .set { ch_scaffolds_db }

        ch_calls = step_3TX_species__vdabricate(ch_scaffolds_db)

        ch_calls.cross(assembled) { extractKey(it) }
            .cross(reference) { extractKey(it) }
            .multiMap { 
                calls: it[0][0]
                assembly: it[0][1]
                reference: it[1][1..3]
            }
            .set { ch_filt }

        step_2AS_filtering__seqio(ch_filt.calls, ch_filt.assembly, ch_filt.reference)
}

workflow {
    ch_input = getInput()
    ch_ref = getReference('fa')
    ch_db = Channel.of([ getDS(), 'viruses_TREF' ])
    module_scaffolds_filtering(ch_input, ch_ref, ch_db)
}