

process alignBowtie2{
  label 'mg03_bowtie'
  conda params.alignBowtie2.conda
  cpus params.resources.alignBowtie2.cpus
  memory params.resources.alignBowtie2.mem
  queue params.resources.alignBowtie2.queue 
  array params.resources.array_size
  clusterOptions params.resources.alignBowtie2.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_bowtie", mode: 'symlink'
  input:
  val index
  val options
  tuple(val(illumina_id), path(fastq))
  
  output:
  tuple(val(illumina_id), path('*.bam'), path('*.err'))
  
  shell:
  '''
  ## bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
  outname=!{illumina_id}'.bam'
  logname=!{illumina_id}'.bowtie2.err'
  bowtie2 -x !{index} --threads !{params.resources.alignBowtie2.cpus} !{options} -1 !{fastq[0]} -2 !{fastq[1]} 2> $logname | \
    samtools view -bSh - -o $outname 
  '''
}