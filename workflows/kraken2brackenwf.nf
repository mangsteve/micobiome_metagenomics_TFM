

include { callKraken2 } from '../modules/kraken2'
include { callBracken } from '../modules/bracken'
include { braken2mpa } from '../modules/bracken2mpa'
include { combineMpa } from '../modules/combinempa'
include { callKronaFromKraken2 } from '../modules/krona'


workflow KRAKEN2BRACKEN {
  take:
  ch_fastq_filtered

  main:
  //Call Kraken2
  callKraken2(params.callKraken2.k2database,
            params.callKraken2.confidence,
            ch_fastq_filtered
  )
  ch_kraken2_output = callKraken2.out
    //.view{"Kraken2 output: $it"}

  //Call Bracken
  ch_bracken_input = ch_kraken2_output
        .map{it -> tuple(it[0], it[2])}
        .combine(Channel.of(['species', 'S'],['genus', 'G'],['phylum', 'P']))
        //.view{"Bracken input: $it"}
  callBracken(params.callKraken2.k2database,
        params.callBracken.threshold,
        params.callBracken.readlen,
        ch_bracken_input
  )
  ch_bracken_output = callBracken.out
    //.view{"Bracken output: $it"}

 //Transform to mpa and merge
  ch_transform2mpa_input = ch_bracken_output
     .map{it -> tuple(it[0], it[1], it[3])}
  braken2mpa(ch_transform2mpa_input)
  ch_transform2mpa_output = braken2mpa.out
    //.view{"Transform to MPA output: $it"}
  ch_combineMpa_input = ch_transform2mpa_output
     .map{it -> tuple(it[1], it[2])}
     .groupTuple()
     .map{it -> tuple(it[0], it[1], it[0][0])} //[1].join(' ') -> it is not necessary to concat files as a string
     //.view{"Combine MPA input: $it"}
  combineMpa(ch_combineMpa_input)
  ch_combineMpa_output = combineMpa.out
     //.view{"Combine MPA output: $it"}
 

  //callKronaFromKraken2: Krona plot from Kraken report
  if(params.resources.callKronaFromKraken2.do_krona){
    ch_krona_input = ch_kraken2_output
        .map{it -> tuple(it[0], it[2])}
    callKronaFromKraken2(ch_krona_input)
    ch_krona_output = callKronaFromKraken2.out
     //.view{ "Krona output: $it" }
  }else{
    ch_krona_output = Channel.from([])
  }
  emit:
  ch_kraken2_output
  ch_bracken_output
  ch_transform2mpa_output
  ch_combineMpa_output
  ch_krona_output
}