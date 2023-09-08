include { multiQC } from '../modules/multiqc'


workflow MULTIQC{
  take:
  ch_fastqc
  ch_fastq_processed
  ch_alignment_output
  ch_kraken2_output
  ch_bracken_output

  main:

  //prepare MultiQC process input
  fastqc_coll = ch_fastqc.collect().ifEmpty([])
  trim_qc = ch_fastq_processed.map{it -> it[4]}.collect().ifEmpty([])
  bowtie2_err = ch_alignment_output.map{it -> it[2]}.collect().ifEmpty([])
  kraken_err = ch_kraken2_output.map{it -> it[2]}.collect().ifEmpty([])
  bracken_err = ch_bracken_output.map{it -> it[4]}.collect().ifEmpty([])

  //Call multiQC process
  multiQC(params.multiQC.configyaml,
            fastqc_coll, 
            trim_qc, 
            bowtie2_err, 
            kraken_err, 
            bracken_err
            )
  ch_multiqc_out = multiQC.out
emit:
ch_multiqc_out
}