process samtoolsSort{
  label 'mg04_sortbam'
  conda params.samtoolsSort.conda
  cpus params.resources.samtoolsSort.cpus
  memory params.resources.samtoolsSort.mem
  queue params.resources.samtoolsSort.queue
  //array params.resources.array_size
  clusterOptions params.resources.samtoolsSort.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg04_sortbam", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(bam))
  
  output:
  tuple(val(illumina_id), path ('*.bam'), path('*.bam.bai'))
  
  shell:
  '''
  outname=$(basename -s .bam !{bam}).sorted.bam
  samtools sort -@ !{params.resources.samtoolsSort.cpus} !{bam} -T samtools.sort.tmp -O bam -o $outname
  samtools index $outname
  '''
}