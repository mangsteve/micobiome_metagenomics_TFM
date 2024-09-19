process removeHumanReads{
  label 'mg05_filtreads'
  conda params.removeHumanReads.conda
  cpus params.resources.removeHumanReads.cpus
  memory params.resources.removeHumanReads.mem
  queue params.resources.removeHumanReads.queue
  array params.resources.array_size
  clusterOptions params.resources.removeHumanReads.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg05_filtreads", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(bam), path(bai), path(fastq))
  
  output:
  tuple(val(illumina_id), path ('*.readids.txt'), path('*.filt.fastq.gz'))
  

  shell:
  '''
  # 1) Get the IDs of reads in bam file (reads mapping human), in different files for R1 and R2

  readsfile=$(basename -s .bam !{bam}).readids.txt
  readsfile_r1=$(basename -s .bam !{bam}).R1.readids.txt
  readsfile_r2=$(basename -s .bam !{bam}).R2.readids.txt
  samtools view !{bam} | cut -f 1 | sort -u > $readsfile
  awk '{print $1"/1"}' $readsfile > $readsfile_r1
  awk '{print $1"/2"}' $readsfile > $readsfile_r2

  # 2) Generate fastq 1 without human reads using seqkit, compress with pigz (parallel gzip)

  oname1=!{illumina_id}_1.filt.fastq.gz
  seqkit grep -v -f $readsfile_r1 !{fastq[0]} | gzip > $oname1 #pigz -p !{params.resources.removeHumanReads.cpus} 
  
  # 3) Generate fastq 2 without human reads using seqkit, compress with pigz (parallel gzip)

  oname2=!{illumina_id}_2.filt.fastq.gz
  seqkit grep -v -f $readsfile_r2 !{fastq[1]} | gzip > $oname2
  '''
}