process braken2mpa{
  label 'mg08_mpa'
  conda params.braken2mpa.conda
  cpus params.resources.braken2mpa.cpus
  memory params.resources.braken2mpa.mem
  queue params.resources.braken2mpa.queue
  //array params.resources.array_size
  clusterOptions params.resources.braken2mpa.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg08_mpa", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(taxonomy_level_name), path(bracken_report))

  output:
  tuple(val(illumina_id), val(taxonomy_level_name), path("*.mpa.txt"))
  

  shell:
  '''
  outfile=$(basename -s .txt !{bracken_report}).mpa.txt
  kreport2mpa.py -r !{bracken_report} -o $outfile --display-header
  '''
}