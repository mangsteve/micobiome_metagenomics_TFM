


include { doTrimmomatic } from '../modules/trimmomatic'
include { getFastQCIllumina        } from '../modules/fastqc'
include { alignBowtie2 } from '../modules/bowtie2'
include { removeHumanReads } from '../modules/removehumanreads'
include { samtoolsSort } from '../modules/samtoolssort'

workflow CLEANFASTQ {
  take: ch_rawfastq
  main:
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

  emit:
  ch_fastq_processed
  ch_fastq_processed_paired
  ch_fastq_filtered_all
  ch_fastq_filtered
  ch_fastqc
  ch_alignment_output
  ch_bam_sorted

}