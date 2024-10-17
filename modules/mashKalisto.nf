process doMashKallistoPipeline {
    label 'mg20_mash_kallisto_pipeline'
    conda params.doMashKallistoPipeline.conda
    cpus params.resources.doMashKallistoPipeline.cpus
    memory params.resources.doMashKallistoPipeline.mem
    queue params.resources.doMashKallistoPipeline.queue
    clusterOptions params.resources.doMashKallistoPipeline.clusterOptions
    errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
    maxRetries 10
    publishDir "$results_dir/mg20_mash_kallisto_pipeline", mode: 'symlink'

    input:
        path mash_output, val(top_strains)  // Archivo de salida de Mash y número de cepas
        tuple(val(illumina_id), path(fastq_paired))  // Secuencias paired-end

    output:
        path 'mash_kallisto_output/*'  // Ruta de salida

    script:
    """
    # Crear un directorio para la salida de Kallisto
    mkdir -p mash_kallisto_output

    # Ejecutar Mash Kallisto Pipeline con las secuencias paired-end
    python3 mash_kallisto_pipeline.py !{mash_output} !{top_strains} \
        --directory mash_kallisto_output \
        --reads1 !{fastq_paired[0]} \
        --reads2 !{fastq_paired[1]}
    """
}
