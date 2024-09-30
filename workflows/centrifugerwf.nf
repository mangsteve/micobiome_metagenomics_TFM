

include { CentrifugerDownload } from '../modules/centrifuger_download'
include { CentrifugerMakeFileList } from '../modules/centrifuger_filelist'
include { CentrifugerBuildDB } from '../modules/centrifuger_build'

workflow CENTRIFUGER {
  take:
  ch_fastq_filtered

  main:

  //Call CentrifugerDownload
  if(params.CentrifugerDownload.do){
    cfgr_download_in = Channel.fromList(params.CentrifugerDownload.argument_list.tokenize(';'))
                                .map{it -> it.tokenize(":")}
                                .map{it -> [false, it[0], it[1], it[2]]}
                                .concat(Channel.from([true, 'taxonomy', 'taxonomy','taxonomy'])
                                                .collect()
                                                )
                                .view{"CentrifugerDownload input: $it"}
    CentrifugerDownload(cfgr_download_in)
    ch_centrifuger_downloads = CentrifugerDownload.out
      .view{"CentrifugerDownload output: $it"}
      .branch{
        taxonomy: it[0] == true
        sequences: it[0] == false
      }
      .set{ ch_centrifuger_downloads_branched }

      ch_centrifuger_downloads_branched.taxonomy.view{"CentrifugerDownload output taxonomy: $it"}
      ch_centrifuger_downloads_branched.sequences.view{"CentrifugerDownload output sequences: $it"}

      // Create either one or several databases
      if(params.CentrifugerMakeFileList.merge_dbs){
      ch_filelist_input = ch_centrifuger_downloads_branched.sequences
        .map{it -> it[2]}
        .collect()
        .map{it -> [params.CentrifugerMakeFileList.merges_dbs_name, it]}
      ch_seq2tax = 
        ch_centrifuger_downloads_branched.sequences
        .map{it -> it[3]}
        .collect()
        .map{it -> [params.CentrifugerMakeFileList.merges_dbs_name, it]}
      ch_filelist_input = ch_filelist_input.join(ch_seq2tax) 
        .view{"CentrifugerMakeFileList input merged: $it"}
      }else{
        ch_filelist_input = ch_centrifuger_downloads_branched.sequences
        .map{it -> [it[1], it[2], it[3]]}
        .view{"CentrifugerMakeFileList input separated: $it"}
      }
      CentrifugerMakeFileList(ch_filelist_input)
      ch_filelist_output = CentrifugerMakeFileList.out
        .view{"CentrifugerMakeFileList output: $it"}

  }else{
    ch_centrifuger_downloads = Channel.from([])
    ch_filelist_output = Channel.from([])
  }

  if(params.CentrifugerDownload.do){
    ch_centrifuger_index = CentrifugerBuildDB(
                      ch_filelist_output,
                      ch_centrifuger_downloads_branched.taxonomy.map{it -> it[2]}
    )
      .view{"CentrifugerBuildDB output: $it"}
  }else{
    ch_centrifuger_index =  Channel.from([])
  }

  emit:
  ch_centrifuger_downloads
  ch_filelist_output
  ch_centrifuger_index
}