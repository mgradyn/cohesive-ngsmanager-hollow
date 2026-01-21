nextflow.enable.dsl=2

include { step_2AS_mapping__bowtie } from '../steps/step_2AS_mapping__bowtie'
include { step_3TX_species__kmerfinder; getBacterialReferencePath } from '../steps/step_3TX_species__kmerfinder'
include { step_4AN_genes__prokka } from '../steps/step_4AN_genes__prokka'
include { step_4AN_AMR__abricate } from '../steps/step_4AN_AMR__abricate'
include { step_4AN_AMR__staramr } from '../steps/step_4AN_AMR__staramr'
include { step_4TY_cgMLST__chewbbaca } from '../steps/step_4TY_cgMLST__chewbbaca'
include { step_4TY_MLST__mlst } from '../steps/step_4TY_MLST__mlst'
include { step_4TY_flaA__flaA } from '../steps/step_4TY_flaA__flaA'
include { csv2map; extractKey; getEmpty } from '../functions/common.nf'
include { getTrimmedReads; getAssembly } from '../functions/parameters.nf'

workflow module_typing_bacteria {
    take: 
        trimmed
        assembly
    main:
        kmer_results = step_3TX_species__kmerfinder(assembly)
        ch_assigned = kmer_results.assigned_species
      
        if (!params.skip_bestref_mapping) {
            trimmed.cross(ch_assigned) { extractKey(it) }
                .multiMap { 
                    trimmed: it[0]
                    species: it[1][1]
                    referencePath: it[1][2]
                }
                .set { ch_mapping }
            
            step_2AS_mapping__bowtie(ch_mapping.trimmed, ch_mapping.referencePath)
        } 

        step_4AN_AMR__abricate(assembly)

        ch_prokka_in = assembly.map { [ it[0], it[1], 'Bacteria', '-', '-', getEmpty() ] }
        step_4AN_genes__prokka(ch_prokka_in)

        assembly.cross(ch_assigned) { extractKey(it) }
            .multiMap { 
                assembly: it[0]
                species: it[1][1]
            }
            .set { ch_typing }

        step_4AN_AMR__staramr(ch_typing.assembly, ch_typing.species)
        step_4TY_MLST__mlst(ch_typing.assembly)
        step_4TY_flaA__flaA(ch_typing.assembly, ch_typing.species)
        step_4TY_cgMLST__chewbbaca(ch_typing.assembly, ch_typing.species, '')
    
    emit:
        genus_species = ch_assigned
}

workflow {
    ch_trimmed = getTrimmedReads(true)
    ch_assembly = getAssembly()
    module_typing_bacteria(ch_trimmed, ch_assembly)
}