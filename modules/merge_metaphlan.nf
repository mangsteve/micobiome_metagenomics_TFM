process mergeMetaphlan{
  label 'mg17_mergemetaphlan'
  conda params.mergeMetaphlan.conda
  cpus params.resources.mergeMetaphlan.cpus
  memory params.resources.mergeMetaphlan.mem
  queue params.resources.mergeMetaphlan.queue
  clusterOptions params.resources.mergeMetaphlan.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg17_mergemetaphlan", mode: 'symlink'
  input:
    path metaphlan_results_list
  output:
  path("metaphlan_merged.tsv")
  
  shell:
  '''
  merge_metaphlan_tables.py *.txt -o metaphlan_merged.tsv  
  '''

  stub:
  """
  echo $humann_results_list
  touch humann3_merged.tsv
  """
  }