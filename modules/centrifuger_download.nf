process CentrifugerDownload{
  label 'mg19_centrifuger_download'
  conda params.CentrifugerDownload.conda
  cpus params.resources.CentrifugerDownload.cpus
  memory params.resources.CentrifugerDownload.mem
  queue params.resources.CentrifugerDownload.queue
  //array params.resources.array_size
  clusterOptions params.resources.CentrifugerDownload.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg19_centrifuger_download", mode: 'symlink'
  input:
    tuple(val(is_taxonomy), val(database), val(arguments), val(output_name))

  output:
    tuple(val(is_taxonomy), val(output_name), path(output_name), path('*seqid2taxid.map'))
  
  script:
  if(is_taxonomy == true)
  """
  centrifuger-download -o $output_name taxonomy
  touch tx_seqid2taxid.map
  """
  else
  """
  centrifuger-download -o $output_name $arguments $database > $output_name'_seqid2taxid.map'
  """
  }