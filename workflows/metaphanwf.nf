include { concatFastq } from '../modules/concatfastq'
include { doMetaphlan } from '../modules/metaphlan'
include { mergeMetaphlan } from '../modules/merge_metaphlan'


workflow METAPHLAN {
take:
    ch_fastq_filtered

main:
    //ch_fastq_filtered.view{ "Humann3 input: $it" }
    concatFastq(ch_fastq_filtered)
    ch_concat_fastq = concatFastq.out
        .view{ "concat fastq output: $it" }

    doMetaphlan(
            params.doMetaphlan.bowtie2db,
            params.doMetaphlan.metaphlan_index, 
            ch_concat_fastq
    )
    ch_metaphlan = doMetaphlan.out
        .view{ "Metaphlan output: $it" }

    ch_metaphlan_join = ch_metaphlan
        .map{it -> tuple( it[1])}
        .flatten()
        .collect()
        .view{ "Metaphlan output join: $it" }
    mergeMetaphlan(ch_metaphlan_join)
    ch_metaphlan_merged = mergeMetaphlan.out
        .view{ "Metaphlan merged output: $it" }


emit:
    ch_metaphlan
    ch_metaphlan_merged
}