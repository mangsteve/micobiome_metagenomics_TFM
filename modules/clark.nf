process doCLARK {
  label 'mg19_clark'
  conda params.doCLARK.conda
  cpus params.resources.doCLARK.cpus
  memory params.resources.doCLARK.mem
  queue params.resources.doCLARK.queue
  clusterOptions params.resources.doCLARK.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg19_clark", mode: 'symlink'

  input:
    val clark_db  // Índice de CLARKK
    tuple(val(illumina_id), path(fastq_paired))

  output:
    tuple(val(illumina_id), path('*_clark_report.txt'))

  shell:
  '''
  output=!{illumina_id}_clark_report.txt

  clark -k !{clark_db} \
    --seq !{fastq_paired[0]} \
    --seq !{fastq_paired[1]} \
    --threads !{params.resources.doCLARK.cpus} \
    > $output
  '''

  stub:
  """
  touch !{illumina_id}_clark_report.txt
  """
}
