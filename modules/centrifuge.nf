process doCentrifuge {
  label 'mg18_centrifuge'
  conda params.doCentrifuge.conda
  cpus params.resources.doCentrifuge.cpus
  memory params.resources.doCentrifuge.mem
  queue params.resources.doCentrifuge.queue
  clusterOptions params.resources.doCentrifuge.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg18_centrifuge", mode: 'symlink'
  
  input:
    path centrifuge_db
    val centrifuge_index
    tuple(val(illumina_id), path(fastq_paired))

  output:
    tuple(val(illumina_id), path('*_centrifuge_report.txt'), path('*_centrifuge_classification.txt'))
  
  shell:
  '''
  report=!{illumina_id}_centrifuge_report.txt
  classification=!{illumina_id}_centrifuge_classification.txt

  centrifuge -x !{centrifuge_db}/!{centrifuge_index} \
    -1 !{fastq_paired[0]} \
    -2 !{fastq_paired[1]} \
    -p !{params.resources.doCentrifuge.cpus} \
    --report-file $report > $classification
  '''

  stub:
  """
  touch !{illumina_id}_centrifuge_report.txt
  touch !{illumina_id}_centrifuge_classification.txt
  """
}

