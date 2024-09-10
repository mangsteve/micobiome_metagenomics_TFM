process doMetaphlan{
  label 'mg16_metaphlan'
  conda params.doMetaphlan.conda
  cpus params.resources.doMetaphlan.cpus
  memory params.resources.doMetaphlan.mem
  queue params.resources.doMetaphlan.queue
  clusterOptions params.resources.doMetaphlan.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg16_metaphlan", mode: 'symlink'
  input:
    path bowtie2db
    val metaphlan_index
    tuple(val(illumina_id), path(fastq_merged))

  output:
    tuple(val(illumina_id), path('*_metaphlan.txt'), path('*_bowtie2.out'))
  
  shell:
  '''
  bwtout=!{illumina_id}_bowtie2.out
  output=!{illumina_id}_metaphlan.txt

  metaphlan !{fastq_merged} \
    --input_type fastq \
    --nproc !{params.resources.doMetaphlan.cpus} \
    --bowtie2out $bwtout \
    --bowtie2db !{bowtie2db} > $output
  '''

  stub:
  '''
  bwtout=!{illumina_id}_bowtie2.out
  output=!{illumina_id}_metaphlan.txt
  
  touch $bwtout
  touch $output
  '''
  }