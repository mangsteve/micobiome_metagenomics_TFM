process callKraken2{
  label 'mg06_kraken2'
  conda params.callKraken2.conda
  cpus params.resources.callKraken2.cpus
  memory params.resources.callKraken2.mem
  queue params.resources.callKraken2.queue
  clusterOptions params.resources.callKraken2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg06_kraken2", mode: 'symlink'
  input:
  path k2database
  val confidence
  tuple(val(illumina_id), val(fastq))

  output:
  tuple(val(illumina_id), path("*standard.kraken2.gz"), path("*.standard.kraken2.report") ,path("*tx.fastq.gz"), path("*kraken2.err"))
  

  shell:
  '''
  outfile=!{illumina_id}.standard.kraken2
  report=!{illumina_id}.standard.kraken2.report
  unclassified=!{illumina_id}.unclassified
  summary=!{illumina_id}.standard.kraken2.err

  kraken2 --db !{k2database} \
        --confidence !{confidence} \
        --threads !{params.resources.callKraken2.cpus} \
        --unclassified-out $unclassified#.tx.fastq \
        --paired !{fastq[0]} !{fastq[1]} \
        --output $outfile \
        --report $report 2> $summary
  pigz -p 4 $unclassified'_1.tx.fastq' $unclassified'_2.tx.fastq'
  pigz -p 4 $outfile
  '''
}