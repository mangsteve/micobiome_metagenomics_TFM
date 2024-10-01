process callKronaFromKraken2{
  label 'mg10_krona'
  conda params.callKronaFromKraken2.conda
  cpus params.resources.callKronaFromKraken2.cpus
  memory params.resources.callKronaFromKraken2.mem
  queue params.resources.callKronaFromKraken2.queue
  //array params.resources.array_size
  clusterOptions params.resources.callKronaFromKraken2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg10_krona", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(kraken_report))

  output:
  tuple(val(illumina_id), path("*.krona.html"))
  

  shell:
  '''
  outfile=!{illumina_id}.krona.html

  ktImportTaxonomy -m 3 -t 5 !{kraken_report} -o $outfile

  '''
}