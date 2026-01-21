nextflow.enable.dsl=2

include { module_denovo } from '../modules/module_denovo'
include { module_scaffolds_filtering } from '../modules/module_scaffolds_filtering'
include { module_draft_genome } from '../modules/module_draft_genome'
include { extractKey } from '../functions/common.nf'
include { getSingleInput;getHost;getReference;getReferenceOptional;getDS } from '../functions/parameters.nf'

workflow module_vdraft {
    take: 
        reads
        host
        reference
        referenceGB
        abricateDatabase
    main:
        denovo_out = module_denovo(reads, host)

        denovo_out.assembled
            .cross(reference) { extractKey(it) }
            .cross(abricateDatabase) { extractKey(it) }
            .multiMap { 
                assembly: it[0][0][0..1]
                reference: it[0][1]
                abricateDatabase: it[1]
            }
            .set { ch_scaff_filt }
            
        module_scaffolds_filtering(ch_scaff_filt.assembly, ch_scaff_filt.reference, ch_scaff_filt.abricateDatabase)
        
        denovo_out.depleted
            .cross(reference) { extractKey(it) }
            .cross(referenceGB) { extractKey(it) }
            .multiMap {
                depleted: it[0][0][0..1]
                reference: it[0][1]
                referenceGB: it[1]
            }
            .set { ch_draft }
            
        module_draft_genome(ch_draft.depleted, ch_draft.reference, ch_draft.referenceGB)
}

workflow {
    ch_in = getSingleInput()
    ch_host = getHost()
    ch_ref_fa = getReference('fa')
    ch_ref_gb = getReferenceOptional('gb')
    ch_db = Channel.of([ getDS(), 'viruses_TREF' ])
    
    module_vdraft(ch_in, ch_host, ch_ref_fa, ch_ref_gb, ch_db)
}