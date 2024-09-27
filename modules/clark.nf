process doCLARKK {
  label 'mg19_clarkk'
  conda params.doCLARKK.conda
  cpus params.resources.doCLARKK.cpus
  memory params.resources.doCLARKK.mem
  queue params.resources.doCLARKK.queue
  clusterOptions params.resources.doCLARKK.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg18_clarkk", mode: 'symlink'

  input:
    val clarkk_db  // Índice de CLARKK
    tuple(val(illumina_id), path(fastq_paired))

  output:
    tuple(val(illumina_id), path('*_clarkk_report.txt'))

  shell:
  '''
  output=!{illumina_id}_clarkk_report.txt

  clark -k !{clarkk_db} \
    --seq !{fastq_paired[0]} \
    --seq !{fastq_paired[1]} \
    --threads !{params.resources.doCLARKK.cpus} \
    > $output
  '''

  stub:
  """
  touch !{illumina_id}_clarkk_report.txt
  """
}
