include { doCLARK } from '../modules/clark'

workflow CLARK {
    take:
        ch_fastq_paired

    main:
        doCLARK(
            params.doCLARK.clark_db,  
            ch_fastq_paired             
        )
        ch_clark_reports = doCLARK.out
        ch_clark_reports.view { "CLARK output: $it" }

    emit:
        ch_clark_reports
}

