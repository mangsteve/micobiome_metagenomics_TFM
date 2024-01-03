process doHumann3{
  label 'mg13_humann3'
  conda params.doHumann3.conda
  cpus params.resources.doHumann3.cpus
  memory params.resources.doHumann3.mem
  queue params.resources.doHumann3.queue
  clusterOptions params.resources.doHumann3.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg13_humann3", mode: 'symlink'
  input:
    path bowtie2db
    val metaphlan_index
    tuple(val(illumina_id), path(fastq_merged))

  output:
    path("${illumina_id}_humann3results/*.tsv")
  
  shell:
  '''
  outdir=!{illumina_id}_humann3results
  humann --input !{fastq_merged} \
    --output  $outdir \
    --threads !{params.resources.doHumann3.cpus}  \
    --input-format "fastq.gz"  \
    --remove-temp-output \
    --metaphlan-options "--input_type fastq --nproc !{params.resources.doHumann3.cpus} --index !{metaphlan_index} --bowtie2db !{bowtie2db}" \
    --resume


  '''

  stub:
  """
  resdir=$illumina_id'_humann3results'
  mkdir \$resdir
  touch \$resdir/$illumina_id'_merged_genefamilies.tsv'
  touch \$resdir/$illumina_id'_merged_pathabundance.tsv'
  touch \$resdir/$illumina_id'_merged_pathcoverage.tsv'
  """
}