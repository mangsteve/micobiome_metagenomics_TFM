process getFastQCIllumina{
  label 'mg02_fastqc'
  conda params.getFastQCIllumina.conda
  cpus params.resources.getFastQCIllumina.cpus
  memory params.resources.getFastQCIllumina.mem
  queue params.resources.getFastQCIllumina.queue
  clusterOptions params.resources.getFastQCIllumina.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_fastqc", mode: 'symlink'
  input:
  path fastq
  
  output:
  tuple(path('*.html'), path('*.zip'))

  shell:
  '''
  fastqc -q !{fastq} --threads !{params.resources.getFastQCIllumina.cpus}
  '''
}