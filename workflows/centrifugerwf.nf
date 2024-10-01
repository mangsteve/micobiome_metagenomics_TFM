

include { CentrifugerDownload } from '../modules/centrifuger_download'
include { CentrifugerMakeFileList } from '../modules/centrifuger_filelist'
include { CentrifugerBuildDB } from '../modules/centrifuger_build'
include { CentrifugerCall } from '../modules/centrifuger_call'

workflow CENTRIFUGER {
  take:
  ch_fastq_filtered

  main:

  //Call CentrifugerDownload
  if(params.CentrifugerDownload.do && params.CentrifugerBuildDB.do){
    cfgr_download_in = Channel.fromList(params.CentrifugerDownload.argument_list.tokenize(';'))
                                .map{it -> it.tokenize(":")}
                                .map{it -> [false, it[0], it[1], it[2]]}
                                .concat(Channel.from([true, 'taxonomy', 'taxonomy','taxonomy'])
                                                .collect()
                                                )
                                //.view{"CentrifugerDownload input: $it"}
    CentrifugerDownload(cfgr_download_in)
    ch_centrifuger_downloads = CentrifugerDownload.out
      //.view{"CentrifugerDownload output: $it"}
      .branch{
        taxonomy: it[0] == true
        sequences: it[0] == false
      }
      .set{ ch_centrifuger_downloads_branched }

      ch_centrifuger_downloads_branched.taxonomy
        .view{"CentrifugerDownload output taxonomy: $it"}
      ch_centrifuger_downloads_branched.sequences
        .view{"CentrifugerDownload output sequences: $it"}

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
          //.view{"CentrifugerMakeFileList input merged: $it"}
      }else{
        ch_filelist_input = ch_centrifuger_downloads_branched.sequences
          .map{it -> [it[1], it[2], it[3]]}
          //.view{"CentrifugerMakeFileList input separated: $it"}
      }

      // Create file(s) with list of fasta files, and merged seq 2 taxid translating map if needed
      CentrifugerMakeFileList(ch_filelist_input)
      ch_filelist_output = CentrifugerMakeFileList.out
        //.view{"CentrifugerMakeFileList output: $it"}
      
      ch_filelist_output = ch_filelist_output.combine(ch_centrifuger_downloads_branched.taxonomy.map{it -> it[2]})
        .view{"CentrifugerBuildDB input (from download): $it"}

  }else{
    // Otherwise, get fasta file list and seq2taxid map from config file
    ch_centrifuger_downloads = Channel.from([])
    if(params.CentrifugerBuildDB.do){
      ch_filelist_output = Channel.from([
        params.CentrifugerBuildDB.index_name])
        .concat(Channel.fromPath(params.CentrifugerBuildDB.filelist))
        .concat(Channel.fromPath(params.CentrifugerBuildDB.seqid2taxid))
        .concat(Channel.fromPath(params.CentrifugerBuildDB.taxonomy))        
      .collect()
      //.map(it -> tuple(it[0], path(it[1]), path(it[3]), path(it[3])))
       .view{"CentrifugerBuildDB input (from files): $it"}
    }else{
      ch_filelist_output = Channel.from([])
    }
  }

  // Build Centrifuger Index (or indices), or get it from config file (only one in that case)
  if(params.CentrifugerBuildDB.do){
    ch_centrifuger_index = CentrifugerBuildDB(ch_filelist_output)
      .view{"CentrifugerBuildDB output: $it"}
  }else{
    ch_centrifuger_index =  Channel.from([
      params.CentrifugerCall.index_name,
      params.CentrifugerCall.index_path
    ]).collect()
    .view{"Centrifuger input (from files): $it"}
  }

  // Call Centrifuger
  if(params.CentrifugerCall.do){
    // Use the 'combine' operator to calculate the product of all samples x all indices
    centrifuge_call_input = ch_centrifuger_index
      .combine(ch_fastq_filtered)
      .view{"Centrifuger input merged: $it"}
    
    CentrifugerCall(centrifuge_call_input)
    centrifuge_call_output = CentrifugerCall.out
      .view{"Centrifuger output: $it"}

  }else{
    centrifuge_call_output = Channel.from([])
  }
  emit:
  ch_centrifuger_downloads
  ch_filelist_output
  ch_centrifuger_index
  centrifuge_call_output
}