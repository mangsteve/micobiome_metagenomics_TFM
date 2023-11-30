# Metagenomics pipeline to classify and annotate Illumina reads

Trims the Illumina reads using Trimmomatic and removes human sequences by aligning with bowtie2. Then, Kraken2 and Bracken are used for species identification, and the Humann3 pipeline for functional annotation. 

# Running the pipeline:

The pipeline is prepared to run locally or in a Slurm cluster. You can switch between configurations by using different config files (in the **/config** directory)

```
nextflow run all.nf -c config/run_samples_cluster.config -profile conda -resume -with-report report.html -with-dag pipeline_dag.html
```

# Dependencies

The pipeline uses conda environments. All the instalation steps can be found in **scripts/create_environments.sh**

# Structure of the repository

Individual processes (e.g., call FastQC, call Kraken2, etc) are in individual files in the **modules/** directory.
Workflows use a set of related processes (e.g., Trimmomatic -> Bowtie2). Each workflow is in an individual file in the **workflows/** directory.

Configuration files are in **config/**, and **sbatch** files to launch the pipeline in a server are in **sbatch/**. The **scripts/** directory contains the original **sbatch** scripts used to create the pipeline. **scripts2/** contains some useful extra scripts. 
