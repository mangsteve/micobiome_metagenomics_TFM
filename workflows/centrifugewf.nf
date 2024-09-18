include { doCentrifuge } from '../modules/centrifuge'



workflow CENTRIFUGE {
    take:
        ch_fastq_filtered

    main:

        doCentrifuge( 
            params.doCentrifuge.centrifuge_index,
            ch_fastq_filtered
        )
        ch_centrifuge = doCentrifuge.out
            .view { "Centrifuge output: $it" }

    emit:
        ch_centrifuge
}
