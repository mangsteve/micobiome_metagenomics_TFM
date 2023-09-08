process doHumann3{
  label 'mg13_humann3'
  conda params.doHumann3.conda
  cpus params.resources.doHumann3.cpus
  memory params.resources.doHumann3.mem
  queue params.resources.doHumann3.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg13_humann3", mode: 'symlink'
  input:
    path bowtie2db
    path metaphlan_index
    tuple(val(illumina_id), path(fastq_merged))

  output:
  path("*_humann3results")
  
  shell:
  '''
  outdir=!{illumina_id}_humann3results
  humann --input !{fastq_merged} \
    --output  $outdir \
    --threads !{params.resources.doHumann3.cpus}  \
    --input-format "fastq"  \
    --metaphlan-options "--input_type fastq --nproc !{params.resources.doHumann3.cpus} --index !{metaphlan_index}  --bowtie2db !{bowtie2db}" \
    --resume
  '''
}