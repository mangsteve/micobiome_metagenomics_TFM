

include { CLEANFASTQ } from './workflows/cleanfastqwf.nf'
include { KRAKEN2BRACKEN } from './workflows/kraken2brackenwf.nf'
include { MULTIQC } from './workflows/multiqcwf.nf'
include { HUMANN3 } from './workflows/humann3wf.nf'


workflow {

  
   ch_rawfastq = Channel.fromFilePairs(params.raw_fastq)

   //Call clean fastq workflow
   CLEANFASTQ(ch_rawfastq)

   //Get outputs (not all of them are used)
   ch_fastq_processed  = CLEANFASTQ.out.ch_fastq_processed
   ch_fastq_processed_paired = CLEANFASTQ.out.ch_fastq_processed_paired
   ch_fastq_filtered_all = CLEANFASTQ.out.ch_fastq_filtered_all
   ch_fastq_filtered = CLEANFASTQ.out.ch_fastq_filtered
   ch_fastqc = CLEANFASTQ.out.ch_fastqc
   ch_alignment_output = CLEANFASTQ.out.ch_alignment_output
   ch_bam_sorted = CLEANFASTQ.out.ch_bam_sorted

   //Call kraken workflow
   KRAKEN2BRACKEN(ch_fastq_filtered)
   //Get outputs (not all of them are used)
   ch_kraken2_output = KRAKEN2BRACKEN.out.ch_kraken2_output
   ch_bracken_output = KRAKEN2BRACKEN.out.ch_bracken_output
   ch_transform2mpa_output = KRAKEN2BRACKEN.out.ch_transform2mpa_output
   ch_combineMpa_output = KRAKEN2BRACKEN.out.ch_combineMpa_output
   ch_krona_output = KRAKEN2BRACKEN.out.ch_krona_output

  //Call Humann3 workflow
   if(params.resources.doHumann3.do_humann){
    HUMANN3(ch_fastq_filtered)
    ch_humann3 = HUMANN3.out.ch_humann3

   }

   //Call MultiQC workflow
   MULTIQC(
     ch_fastqc,
     ch_fastq_processed,
     ch_alignment_output,
     ch_kraken2_output,
     ch_bracken_output
   )
   ch_multiqc_out = MULTIQC.out.ch_multiqc_out
}