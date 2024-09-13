include { concatFastq } from '../modules/concatfastq'
include { doHumann3 } from '../modules/humann3'
include { mergeHumann } from '../modules/merge_humann'
include { translateHumann } from '../modules/translate_humann'

workflow HUMANN3 {
take:
ch_fastq_filtered

main:
//ch_fastq_filtered.view{ "Humann3 input: $it" }
concatFastq(ch_fastq_filtered)
ch_concat_fastq = concatFastq.out
    //.view{ "concat fastq output: $it" }
doHumann3(
        params.doHumann3.bowtie2db,
        params.doHumann3.metaphlan_index, 
        ch_concat_fastq
)
ch_humann3 = doHumann3.out
    //.view{ "Humann3 output: $it" }
    .flatten()
    .collect()
    //.view{ "Humann3 output flat: $it" }

mergeHumann(ch_humann3)
ch_humann3_merged = mergeHumann.out
    //.view{ "Humann3 output merged: $it" }

translateHumann(ch_humann3_merged, params.translateHumann3.cazy_db)
ch_humann3_translated = translateHumann.out
    //.view{ "Humann3 output translated: $it" }

emit:
ch_humann3
ch_humann3_merged
ch_humann3_translated
} 