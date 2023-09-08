process doTrimmomatic{
  label 'mg01_trimmomatic'
  conda params.doTrimmomatic.conda
  cpus params.resources.doTrimmomatic.cpus
  memory params.resources.doTrimmomatic.mem
  queue params.resources.doTrimmomatic.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg01_trimmomatic", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(fastq))
  val illuminaclip
  val slidingwindow
  val minlen
  
  output:
  tuple(val(illumina_id), path('*.trimmed.fastq.gz'), path('*.trimmed.single.fastq.gz'), path('*.trimlog'), path('*.trimlog.out'))

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
        MINLEN:!{minlen} 2> !{illumina_id}.trimlog.out
  '''
}

process getFastQCIllumina{
  label 'mg02_fastqc'
  conda params.getFastQCIllumina.conda
  cpus params.resources.getFastQCIllumina.cpus
  memory params.resources.getFastQCIllumina.mem
  queue params.resources.getFastQCIllumina.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg02_fastqc", mode: 'symlink'
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
  label 'mg03_bowtie'
  conda params.alignBowtie2.conda
  cpus params.resources.alignBowtie2.cpus
  memory params.resources.alignBowtie2.mem
  queue params.resources.alignBowtie2.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg03_bowtie", mode: 'symlink'
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
  label 'mg04_sortbam'
  conda params.samtoolsSort.conda
  cpus params.resources.samtoolsSort.cpus
  memory params.resources.samtoolsSort.mem
  queue params.resources.samtoolsSort.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg04_sortbam", mode: 'symlink'
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
  label 'mg05_filtreads'
  conda params.removeHumanReads.conda
  cpus params.resources.removeHumanReads.cpus
  memory params.resources.removeHumanReads.mem
  queue params.resources.removeHumanReads.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg05_filtreads", mode: 'symlink'
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


process callKraken2{
  label 'mg06_kraken2'
  conda params.callKraken2.conda
  cpus params.resources.callKraken2.cpus
  memory params.resources.callKraken2.mem
  queue params.resources.callKraken2.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg06_kraken2", mode: 'symlink'
  input:
  path k2database
  val confidence
  tuple(val(illumina_id), val(fastq))

  output:
  tuple(val(illumina_id), path("*standard.kraken2"), path("*.standard.kraken2.report") ,path("*tx.fastq.gz"), path("*kraken2.err"))
  

  shell:
  '''
  outfile=!{illumina_id}.standard.kraken2
  report=!{illumina_id}.standard.kraken2.report
  unclassified=$(basename -s .fastq.gz !{fastq[0]})
  summary=!{illumina_id}.standard.kraken2.err

  kraken2 --db !{k2database} \
        --confidence !{confidence} \
        --threads !{params.resources.callKraken2.cpus} \
        --unclassified-out $unclassified#.tx.fastq \
        --paired !{fastq[0]} !{fastq[1]} \
        --output $outfile \
        --report $report 2> $summary
  pigz -p 4 $unclassified'_1.tx.fastq' $unclassified'_2.tx.fastq'

  '''
}


process callBracken{
  label 'mg07_Bracken'
  conda params.callBracken.conda
  cpus params.resources.callBracken.cpus
  memory params.resources.callBracken.mem
  queue params.resources.callBracken.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg07_Bracken", mode: 'symlink'
  input:
  path k2database
  val threshold
  val readlen
  tuple(val(illumina_id), path(k2report), val(taxonomy_level_name), val(taxonomy_level_param))

  output:
  tuple(val(illumina_id), val(taxonomy_level_name), path("*.bracken.txt"), path("*.bracken.report.txt"), path("*bracken.err"))
  

  shell:
  '''
  outfile=!{illumina_id}.!{taxonomy_level_name}.bracken.txt
  report=!{illumina_id}.!{taxonomy_level_name}.bracken.report.txt
  summary=!{illumina_id}.!{taxonomy_level_name}.bracken.err

  bracken -d !{k2database} \
        -i !{k2report} \
        -o $outfile \
        -w $report \
        -r !{readlen}  \
        -t !{threshold} \
        -l !{taxonomy_level_param} > $summary

  '''
}

process braken2mpa{
  label 'mg08_mpa'
  conda params.braken2mpa.conda
  cpus params.resources.braken2mpa.cpus
  memory params.resources.braken2mpa.mem
  queue params.resources.braken2mpa.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg08_mpa", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(taxonomy_level_name), path(bracken_report))

  output:
  tuple(val(illumina_id), val(taxonomy_level_name), path("*.mpa.txt"))
  

  shell:
  '''
  outfile=$(basename -s .txt !{bracken_report}).mpa.txt
  kreport2mpa.py -r !{bracken_report} -o $outfile --display-header
  '''
}


process combineMpa{
  label 'mg09_combinempa'
  conda params.combineMpa.conda
  cpus params.resources.combineMpa.cpus
  memory params.resources.combineMpa.mem
  queue params.resources.combineMpa.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg09_combinempa", mode: 'copy'
  input:
  tuple(val(taxonomy_level_name), val(mpa), val(letter))

  output:
  tuple(val(taxonomy_level_name), path("*.mpa.combined.txt"), path("*.mpa.combined.clean1.txt"), path("*.mpa.combined.clean2.txt"))
  

  shell:
  '''
  outfile=!{taxonomy_level_name}.mpa.combined.txt
  outfile2=!{taxonomy_level_name}.mpa.combined.clean1.txt
  outfile3=!{taxonomy_level_name}.mpa.combined.clean2.txt

  combine_mpa.py -i !{mpa} -o $outfile
  grep -E "(!{letter}__)|(#Classification)" $outfile > $outfile2
  sed -e 's/.species.bracken.report.txt//g' $outfile2 > $outfile3

  '''
}

process callKronaFromKraken2{
  label 'mg10_krona'
  conda params.callKronaFromKraken2.conda
  cpus params.resources.callKronaFromKraken2.cpus
  memory params.resources.callKronaFromKraken2.mem
  queue params.resources.callKronaFromKraken2.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg10_krona", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(kraken_report))

  output:
  tuple(val(illumina_id), path("*.krona.html"))
  

  shell:
  '''
  outfile=!{illumina_id}.krona.html

  ktImportTaxonomy -m 3 -t 5 !{kraken_report} -o $outfile

  '''
}

process multiQC{
  label 'mg11_multiqc'
  conda params.multiQC.conda
  cpus params.resources.multiQC.cpus
  memory params.resources.multiQC.mem
  queue params.resources.multiQC.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg10_multiqc", mode: 'copy'
  input:
    path yaml
    path ch_fastqc
    path trim_qc
    path bowtie2_err
    path kraken_err
    path bracken_err

  output:
  path("multiqc_report.html")
  path("multiqc_data")
  

  shell:
  '''
  mv !{yaml} multiqc_config.yaml
  multiqc .  
  '''
}

process concatFastq{
  label 'mg12_concat_fastq'
  conda params.concatFastq.conda
  cpus params.resources.concatFastq.cpus
  memory params.resources.concatFastq.mem
  queue params.resources.concatFastq.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg12_concat_fastq", mode: 'symlink'
  input:
  tuple(val(illumina_id), val(fastq))

  output:
  tuple(val(illumina_id), path("*_merged.fastq.gz"))
  
  shell:
  '''
  merged_fastq=!{illumina_id}_merged.fastq.gz
  zcat !{fastq[0]} !{fastq[1]} | pigz -p !{params.resources.doHumann3.cpus} > $merged_fastq
  '''
}

process doHumann3{
  label 'mg13_humann3'
  conda params.doHumann3.conda
  cpus params.resources.doHumann3.cpus
  memory params.resources.doHumann3.mem
  queue params.resources.doHumann3.queue
  errorStrategy { task.exitStatus in 1..2 ? 'retry' : 'ignore' }
  maxRetries 10
  publishDir "$results_dir/mg13_humann3", mode: 'symlink'
  input:
    path bowtie2db
    path metaphlan_index
    tuple(val(illumina_id), path(fastq_merged))

  output:
  path("*_humann3results")
  
  shell:
  '''
  outdir=!{illumina_id}_humann3results
  humann --input !{fastq_merged} \
    --output  $outdir \
    --threads !{params.resources.doHumann3.cpus}  \
    --input-format "fastq"  \
    --metaphlan-options "--input_type fastq --nproc !{params.resources.doHumann3.cpus} --index !{metaphlan_index}  --bowtie2db !{bowtie2db}" \
    --resume
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

  //Call Kraken2
  callKraken2(params.callKraken2.k2database,
            params.callKraken2.confidence,
            ch_fastq_filtered
  )
  ch_kraken2_output = callKraken2.out
   // .view{"Kraken2 output: $it"}

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
     .map{it -> tuple(it[0], it[1].join(' '), it[0][0])}
     //.view{"Combine MPA input: $it"}
  combineMpa(ch_combineMpa_input)
  ch_combineMpa_output = combineMpa.out
     //.view{"Combine MPA output: $it"}

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

  //callKronaFromKraken2: Krona plot from Kraken report
  if(params.resources.callKronaFromKraken2.do_krona){
    ch_krona_input = ch_kraken2_output
        .map{it -> tuple(it[0], it[2])}
    callKronaFromKraken2(ch_krona_input)
    ch_krona_output = callKronaFromKraken2.out
     //.view{ "Krona output: $it" }
  }

  //MultiQC
  fastqc_coll = ch_fastqc.collect().ifEmpty([])
  trim_qc = ch_fastq_processed.map{it -> it[4]}.collect().ifEmpty([])
  bowtie2_err = ch_alignment_output.map{it -> it[2]}.collect().ifEmpty([])
  kraken_err = ch_kraken2_output.map{it -> it[2]}.collect().ifEmpty([])
  bracken_err = ch_bracken_output.map{it -> it[4]}.collect().ifEmpty([])
  multiQC(params.multiQC.configyaml,
            fastqc_coll, 
            trim_qc, 
            bowtie2_err, 
            kraken_err, 
            bracken_err
            )

//Humann3. First concatenate fastq and then call the humann3 pipeline
//ch_fastq_filtered.view{ "Humann3 input: $it" }
if(params.resources.doHumann3.do_humann){

    concatFastq(ch_fastq_filtered)
    ch_concat_fastq = concatFastq.out
        //.view{ "concat fastq output: $it" }

    doHumann3(
        params.doHumann3.bowtie2db,
        params.doHumann3.metaphlan_index, 
        ch_concat_fastq
    )
    ch_humann3 = doHumann3.out
        .view{ "Humann3 output: $it" }
}
}