process translateHumann{
  label 'mg15_humanntranslated'
  conda params.doHumann3.conda
  cpus params.resources.translateHumann3.cpus
  memory params.resources.translateHumann3.mem
  queue params.resources.translateHumann3.queue
  clusterOptions params.resources.translateHumann3.clusterOptions
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg15_translateHumann", mode: 'copy'
  input:
    path humann_results_merged
    path cazy_db
  output:
    path("*.tsv")

  shell:
  '''
  # Regroup terms
  humann_regroup_table -i !{humann_results_merged} -o humann3_merged_genetables_GO.tsv -g uniref90_go
  humann_regroup_table -i !{humann_results_merged} -o humann3_merged_genetables_KO.tsv -g uniref90_ko
  humann_regroup_table -i !{humann_results_merged} -o humann3_merged_genetables_RXN.tsv -g uniref90_rxn
  humann_regroup_table -i !{humann_results_merged} -o humann3_merged_genetables_CAZY.tsv -c !{cazy_db}
  
  # Renorm
  humann_renorm_table --input !{humann_results_merged} --output humann3_merged_abundances_CPM.tsv --units cpm --update-snames
  humann_renorm_table --input humann3_merged_genetables_GO.tsv --output humann3_merged_genetables_GO_CPM.tsv --units cpm --update-snames
  humann_renorm_table --input humann3_merged_genetables_KO.tsv --output humann3_merged_genetables_KO_CPM.tsv --units cpm --update-snames
  humann_renorm_table --input humann3_merged_genetables_RXN.tsv --output humann3_merged_genetables_RXN_CPM.tsv --units cpm --update-snames
  humann_renorm_table --input humann3_merged_genetables_CAZY.tsv --output humann3_merged_genetables_CAZY_CPM.tsv --units cpm --update-snames
  
  # Add names
  humann_rename_table --input humann3_merged_abundances_CPM.tsv --output humann3_merged_abundances_CPM_named.tsv --names metacyc-pwy
  humann_rename_table --input humann3_merged_genetables_RXN_CPM.tsv --output humann3_merged_genetables_RXN_CPM_named.tsv --names metacyc-rxn
  humann_rename_table --input humann3_merged_genetables_KO_CPM.tsv --output humann3_merged_genetables_KO_CPM_named.tsv --names kegg-pathway
  humann_rename_table --input humann3_merged_genetables_GO_CPM.tsv --output humann3_merged_genetables_GO_CPM_named.tsv --names go
  
  '''

  stub:
  """
  # Regroup terms
  touch humann3_merged_genetables_GO.tsv
  touch humann3_merged_genetables_KO.tsv
  touch humann3_merged_genetables_RXN.tsv
  touch humann3_merged_genetables_CAZY.tsv
  
  # Renorm
  touch humann3_merged_abundances_CPM.tsv
  touch humann3_merged_genetables_GO_CPM.tsv
  touch humann3_merged_genetables_KO_CPM.tsv
  touch humann3_merged_genetables_RXN_CPM.tsv
  touch humann3_merged_genetables_CAZY_CPM.tsv
  
  # Add names
  touch humann3_merged_genetables_RXN_CPM_named.tsv
  touch humann3_merged_genetables_KO_CPM_named.tsv
  touch humann3_merged_genetables_GO_CPM_named.tsv
  touch humann3_merged_abundances_CPM_named.tsv
  """
}