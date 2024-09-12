include { concatFastq } from '../modules/concatfastq'
include { doCentrifuge } from '../modules/centrifuge'



workflow CENTRIFUGE {
    take:
        ch_fastq_filtered

    main:
        concatFastq(ch_fastq_filtered)
        ch_concat_fastq = concatFastq.out
        .view{ "concat fastq output: $it" }

        doCentrifuge(
            params.doCentrifuge.centrifuge_db, 
            params.doCentrifuge.centrifuge_index,
            ch_fastq_filtered
        )
        ch_centrifuge = doCentrifuge.out
            .view { "Centrifuge output: $it" }

    emit:
        ch_centrifuge
}
