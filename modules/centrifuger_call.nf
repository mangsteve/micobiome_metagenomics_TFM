process CentrifugerCall{
  label 'mg22_centrifuger_call'
  conda params.CentrifugerCall.conda
  cpus params.resources.CentrifugerCall.cpus
  memory params.resources.CentrifugerCall.mem
  queue params.resources.CentrifugerCall.queue
  //array params.resources.array_size
  clusterOptions params.resources.CentrifugerCall.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg22_centrifuger_call", mode: 'symlink'
  input:
    tuple(val(ref_name), path(ref_dir), val(sample_id), path(fastq))

  output:
    tuple(path("*cfgr.out"), path("*cfgr.err"))
  
  shell:
  '''
  centrifuger -1 !{fastq[1]} -2 !{fastq[1]} \
    -x !{ref_dir}/!{ref_name} !{params.CentrifugerCall.options} \
    > !{sample_id}_!{ref_name}.cfgr.out 2> !{sample_id}_!{ref_name}.cfgr.err
  '''
  }