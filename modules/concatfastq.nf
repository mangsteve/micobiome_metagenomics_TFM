process concatFastq{
  label 'mg12_concat_fastq'
  conda params.concatFastq.conda
  cpus params.resources.concatFastq.cpus
  memory params.resources.concatFastq.mem
  queue params.resources.concatFastq.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg12_concat_fastq", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(fastq))

  output:
  tuple(val(illumina_id), path("*_merged.fastq.gz"))
  
  shell:
  '''
  merged_fastq=!{illumina_id}_merged.fastq.gz
  zcat !{fastq[0]} !{fastq[1]} | pigz -p !{params.resources.doHumann3.cpus} > $merged_fastq
  '''
}