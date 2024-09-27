workflow CLARKK {
    take:
        ch_fastq_paired

    main:
        doCLARKK(
            params.doCLARKK.clarkk_db,  // Base de datos CLARKK
            ch_fastq_paired             // Pares de archivos FASTQ
        )
        ch_clarkk_reports = doCLARKK.out
        ch_clarkk_reports.view { "CLARKK output: $it" }

    emit:
        ch_clarkk_reports
}

