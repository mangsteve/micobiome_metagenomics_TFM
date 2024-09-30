process CentrifugerMakeFileList{
  label 'mg20_centrifuger_filelist'
  conda params.CentrifugerMakeFileList.conda
  cpus params.resources.CentrifugerMakeFileList.cpus
  memory params.resources.CentrifugerMakeFileList.mem
  queue params.resources.CentrifugerMakeFileList.queue
  //array params.resources.array_size
  clusterOptions params.resources.CentrifugerMakeFileList.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg20_centrifuger_filelist", mode: 'symlink'
  input:
    tuple(val(ref_name), path(ref_directory))

  output:
    path('*.txt')
  
  shell:
  '''
  for p in !{ref_directory}; do
    ls $p'/*/*.fna.gz' >> !{ref_name}_file_list .txt
  done
  '''
  }