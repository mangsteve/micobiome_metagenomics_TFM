

include { CentrifugerDownload } from '../modules/centrifuger_download'
include { CentrifugerMakeFileList } from '../modules/centrifuger_filelist'

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
      if(params.CentrifugerDownload.merge_fasta){
      ch_filelist_input = ch_centrifuger_downloads_branched.sequences
        .map{it -> it[2]}
        .collect()
        .map{it -> ["mergeddb", it]}
        .view{"CentrifugerMakeFileList input merged: $it"}
      }else{
        ch_filelist_input = ch_centrifuger_downloads_branched.sequences
        .map{it -> [it[0], it[2]]}
        .view{"CentrifugerMakeFileList input separated: $it"}
      }
      CentrifugerMakeFileList(ch_filelist_input)
      ch_filelist_output = CentrifugerMakeFileList.out

  }else{
    ch_centrifuger_downloads = Channel.from([])
    ch_filelist_output = Channel.from([])
  }


  ch_centrifuger_index = Channel.from([])

  emit:
  ch_centrifuger_downloads
  ch_centrifuger_index
}