process doTrimmomatic{
  label 'fq01_trimmomatic'
  conda params.doTrimmomatic.conda
  cpus params.resources.doTrimmomatic.cpus
  memory params.resources.doTrimmomatic.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/fq01_trimmomatic", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(fastq))
  val illuminaclip
  val slidingwindow
  val minlen
  
  output:
  tuple(val(illumina_id), path('*.trimmed.fastq.gz'), path('*.trimmed.single.fastq.gz'), path('*.trimlog'))

  shell:
  '''
  trimmomatic PE -threads !{params.resources.doTrimmomatic.cpus} -phred33 \
        -trimlog !{illumina_id}.trimlog \
        !{fastq[0]} !{fastq[1]} \
        !{illumina_id}_1.trimmed.fastq.gz \
        !{illumina_id}_1.trimmed.single.fastq.gz \
        !{illumina_id}_2.trimmed.fastq.gz \
        !{illumina_id}_2.trimmed.single.fastq.gz \
        ILLUMINACLIP:!{illuminaclip} \
        SLIDINGWINDOW:!{slidingwindow} \
        MINLEN:!{minlen}
  '''
}

process getFastQCIllumina{
  label 'fq02_fastqc'
  conda params.getFastQCIllumina.conda
  cpus params.resources.getFastQCIllumina.cpus
  memory params.resources.getFastQCIllumina.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/fq02_fastqc", mode: 'symlink'
  input:
  path fastq
  
  output:
  tuple(path('*.html'), path('*.zip'))

  shell:
  '''
  fastqc -q !{fastq} --threads !{params.resources.getFastQCIllumina.cpus}
  '''
}

process alignBowtie2{
  label 'fq03_bowtie'
  conda params.alignBowtie2.conda
  cpus params.resources.alignBowtie2.cpus
  memory params.resources.alignBowtie2.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/fq03_bowtie", mode: 'symlink'
  input:
  val index
  val options
  tuple(val(illumina_id), val(fastq))
  
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

process samtoolsSort{
  label 'fq04_sortbam'
  conda params.samtoolsSort.conda
  cpus params.resources.samtoolsSort.cpus
  memory params.resources.samtoolsSort.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/fq04_sortbam", mode: 'symlink'
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

process removeHumanReads{
  label 'fq05_filtreads'
  conda params.removeHumanReads.conda
  cpus params.resources.removeHumanReads.cpus
  memory params.resources.removeHumanReads.mem
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/fq05_filtreads", mode: 'symlink'
  input:
  tuple(val(illumina_id), path(bam), path(bai), val(fastq))
  
  output:
  tuple(val(illumina_id), path ('*.readids.txt'), path('*.filt.fastq.gz'))
  

  shell:
  '''
  # 1) Get the IDs of reads in bam file (reads mapping human)
  readsfile=$(basename -s .bam !{bam}).readids.txt
  samtools view !{bam} | cut -f 1 | sort -u > $readsfile

  # 2) Generate fastq 1 without human reads using seqkit, compress with pigz (parallel gzip)
  oname1=$(basename -s .fastq.gz !{fastq[0]}).filt.fastq.gz
  seqkit grep -v -f $readsfile !{fastq[0]} | pigz -p !{params.resources.removeHumanReads.cpus} > $oname1
  
  # 3) Generate fastq 2 without human reads using seqkit, compress with pigz (parallel gzip)
  oname2=$(basename -s .fastq.gz !{fastq[1]}).filt.fastq.gz
  seqkit grep -v -f $readsfile !{fastq[1]} | pigz -p !{params.resources.removeHumanReads.cpus} > $oname2
  '''
}


workflow {

  ch_rawfastq = Channel
  .fromFilePairs(params.raw_fastq)
  //.view{ "input fastq files: $it" }

 //Trim reads
  ch_fastq_processed = ch_rawfastq
  if(params.doTrimmomatic.do_trim){
    doTrimmomatic(ch_rawfastq, 
                params.doTrimmomatic.illuminaclip, 
                params.doTrimmomatic.slidingwindow,
                params.doTrimmomatic.minlen
                )
    ch_fastq_processed = doTrimmomatic.out
     //.view{ "Illumina trimmed reads: $it" }
  }

  //Map paired reads
  ch_fastq_processed_paired = ch_fastq_processed.map{it -> tuple(it[0], it[1])}
  //.view{ "trimmed fastq paired only: $it" }
  ch_alignment_output = Channel.from([])
  ch_bam = Channel.from([])
  if(params.mapping_tool == 'bowtie2'){
        alignBowtie2(params.alignBowtie2.index, 
                    params.alignBowtie2.options,
                    ch_fastq_processed_paired)
        ch_alignment_output = alignBowtie2.out
        ch_bam = ch_alignment_output.map{it -> tuple(it[0], it[1])}
  }

  //Sort and index bam with human reads
  samtoolsSort(ch_bam)
  ch_bam_sorted = samtoolsSort.out
    //.view{ "sorted human bam: $it" }

  //Remove Human reads
  ch_fastq_filtered = Channel.from([])
  if(params.removeHumanReads.remove_human_reads){
        ch_bam_and_fastq = ch_bam_sorted.join(ch_fastq_processed_paired)
          //.view{"join bam and trimmed fastq paired: $it" }
        removeHumanReads(ch_bam_and_fastq)
        ch_fastq_filtered_all = removeHumanReads.out
            //.view{"filtered fastq and list: $it"}
        ch_fastq_filtered = ch_fastq_filtered_all.map{it -> tuple(it[0], it[2])}
            //.view{"filtered fastq only: $it"}
  }else{
        ch_fastq_filtered = ch_fastq_processed_paired
  }


  // FastQC STEPS
  //Initialize channels 
  ch_fastqc = Channel.from([])
  ch_flatfastq = Channel.from([])

  //Get a flattened list of raw fastq
  if(params.getFastQCIllumina.do_fastqc_raw){
    ch_flatfastq = ch_rawfastq.map{it -> it[1]}.flatten()
  }

  //Add trimmed fastq to the flat fastq channel to perform FastQC
    if(params.getFastQCIllumina.do_fastqc_trim){
      ch_flatfastq = ch_fastq_processed.map{it -> it[1]}
      .flatten().concat(ch_flatfastq)
    }
    if(params.getFastQCIllumina.do_fastqc_trim_single){
      ch_flatfastq = ch_fastq_processed.map{it -> it[2]}
      .flatten().concat(ch_flatfastq)
    }

  // Do FastQC for all fastq files
  if(params.getFastQCIllumina.do_fastqc_raw | params.getFastQCIllumina.do_fastqc_trim){
    getFastQCIllumina(ch_flatfastq)
    //.view{ "getFastQCIllumina - FastQC reports: $it" }
    ch_fastqc = getFastQCIllumina.out
  }

}