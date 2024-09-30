process CentrifugerBuildDB{
  label 'mg21_centrifuger_build'
  conda params.CentrifugerBuildDB.conda
  cpus params.resources.CentrifugerBuildDB.cpus
  memory params.resources.CentrifugerBuildDB.mem
  queue params.resources.CentrifugerBuildDB.queue
  //array params.resources.array_size
  clusterOptions params.resources.CentrifugerBuildDB.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg21_centrifuger_build", mode: 'symlink'
  input:
    tuple(val(ref_name), path(file_list), path(seqid2taxid))
    path taxonomy 

  output:
    tuple(val(ref_name), path(ref_name))
  
  shell:
  '''
  mkdir !{ref_name}

  mem=$(echo "!{params.resources.CentrifugerBuildDB.mem}" | cut -f1 -d' ')'G'
  centrifuger-build -t !{params.resources.CentrifugerBuildDB.cpus} \
    --conversion-table !{seqid2taxid} \
	--taxonomy-tree !{taxonomy}/nodes.dmp --name-table !{taxonomy}/names.dmp \
	-l !{file_list} -o !{ref_name}/!{ref_name} --build-mem $mem

  
  '''
  stub:
  '''
  mkdir !{ref_name}
  touch !{ref_name}/!{ref_name}1.cfr
  touch !{ref_name}/!{ref_name}2.cfr
  touch !{ref_name}/!{ref_name}3.cfr
  touch !{ref_name}/!{ref_name}4.cfr
  '''
  }