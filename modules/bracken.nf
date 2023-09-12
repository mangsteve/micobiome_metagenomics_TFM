process callBracken{
  label 'mg07_Bracken'
  conda params.callBracken.conda
  cpus params.resources.callBracken.cpus
  memory params.resources.callBracken.mem
  queue params.resources.callBracken.queue
  clusterOptions params.resources.callBracken.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg07_Bracken", mode: 'symlink'
  input:
  path k2database
  val threshold
  val readlen
  tuple(val(illumina_id), path(k2report), val(taxonomy_level_name), val(taxonomy_level_param))

  output:
  tuple(val(illumina_id), val(taxonomy_level_name), path("*.bracken.txt"), path("*.bracken.report.txt"), path("*bracken.err"))
  

  shell:
  '''
  outfile=!{illumina_id}.!{taxonomy_level_name}.bracken.txt
  report=!{illumina_id}.!{taxonomy_level_name}.bracken.report.txt
  summary=!{illumina_id}.!{taxonomy_level_name}.bracken.err

  bracken -d !{k2database} \
        -i !{k2report} \
        -o $outfile \
        -w $report \
        -r !{readlen}  \
        -t !{threshold} \
        -l !{taxonomy_level_param} > $summary

  '''
}