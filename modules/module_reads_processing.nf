nextflow.enable.dsl=2

include { step_0SQ_rawreads__fastq } from '../steps/step_0SQ_rawreads__fastq'
include { step_1PP_trimming__trimmomatic } from '../steps/step_1PP_trimming__trimmomatic'
include { step_1PP_trimming__fastp } from '../steps/step_1PP_trimming__fastp'
include { step_3TX_class__kraken } from '../steps/step_3TX_class__kraken'
include { getInput;hasFastqData;hasEnoughFastqData;isIlluminaPaired;isIonTorrent;isNanopore } from '../functions/parameters.nf'
include { isBacterium } from '../functions/sampletypes.nf'

workflow module_reads_processing {
    take: 
        rawReads
    main:
        rawReads.branch {
            with_data: hasFastqData(it[1])
            no_reads: true
        }.set { ch_raw }

        step_0SQ_rawreads__fastq(ch_raw.with_data)        

        ch_raw.with_data.branch {
            illumina: isIlluminaPaired(it[1])
            ion: isIonTorrent(it[1])
            nanopore: isNanopore(it[1])
            other: true 
        }.set { ch_tech }

        ch_tech.illumina.branch {
            bacteria: isBacterium(it)
            other: true 
        }.set { ch_illumina }

        proc_trimmomatic = step_1PP_trimming__trimmomatic(ch_illumina.other)
        
        ch_fastp_in = ch_tech.ion.mix(ch_illumina.bacteria)
        proc_fastp = step_1PP_trimming__fastp(ch_fastp_in)

        ch_trimmed = proc_trimmomatic.trimmed.mix(proc_fastp.trimmed)

        ch_trimmed.branch {
            with_data: hasEnoughFastqData(it[1])
            insufficient: true
        }.set { ch_final }

        step_3TX_class__kraken(ch_final.with_data)

    emit:
        no_reads = ch_raw.no_reads
        trimmed_with_data = ch_final.with_data
        insufficient_number_of_reads = ch_final.insufficient
}

workflow {
    module_reads_processing(getInput())
}