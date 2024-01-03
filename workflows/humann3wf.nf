include { concatFastq } from '../modules/concatfastq'
include { doHumann3 } from '../modules/humann3'
//include { mergeHumann } from '../modules/merge_humann'

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
    .view{ "Humann3 output: $it" }

ch_humann3_2merge = ch_humann3.collect{ it->"${it}/*.tsv" }
    .view{ "Humann3 output by file: $it" }
    .toList()
    .map{ it->it.fromPath() }
    .view{ "Humann3 output fromPath: $it" }
    .view{ it.getClass() }

/* ch_humann3_2merge = Channel.fromPath(ch_humann3_2merge)
    .view{ "Humann3 merge: $it" } */
/* mergeHumann(ch_humann3_2merge)
ch_humann3_merged = mergeHumann.out */

emit:
ch_humann3
//ch_humann3_merged
}