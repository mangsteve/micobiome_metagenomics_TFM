workflow MASH_KALLISTO {
    take:
    ch_fastq_paired
    ch_mash_output

    main:
    // Ejecutar Mash Kallisto Pipeline
    doMashKallistoPipeline(
        ch_mash_output,  // Archivo de salida de Mash
        params.doMashKallistoPipeline.top_strains,  // Número de cepas a seleccionar
        ch_fastq_paired  // Secuencias paired-end
    )

    emit:
    doMashKallistoPipeline.out
}
