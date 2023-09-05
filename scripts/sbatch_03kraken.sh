#!/bin/bash

#################################################################################################################
## script para el quality control de los ficheros                                                              ##
##                                                                                                             ##
## SBATCH --job-name=name: Nombre del trabajo                                                                  ##
## SBATCH --output=dir/name_%A_%a.out #archivo_%j.out: Directorio y nombre del fichero de salida               ##
## SBATCH --error=dir/name_%A_%a.err #archivo_%j.err: Directorio y nombre del ficherro de errores              ##
## SBATCH --array=1-n: paralelizacion del trabajo, lanza simultaneamente 16 trabajos                           ##
## SBATCH --partition=long: la particion que queremos usar                                                     ##
##          "slimits" para ver las opciones                                                                    ##
## SBATCH --mem 10G                                                                                            ##
##                                                                                                             ##
##                          Parámetros del sistema de colas que utilizó Teresa                                 ##
##                                                                                                             ##
## #SBATCH --array=222-442 #1-221                                                                              ##
## SBATCH --partition=bigmem                                                                                   ##
## SBATCH --cpus-per-task 10                                                                                   ##
## SBATCH --mem=100G                                                                                           ##
##                                                                                                             ##
#################################################################################################################

#SBATCH --job-name=CovidPersistente_kraken
#SBATCH --output=/home/ysanz/CovidPersistente/03_kraken/CovidPersistente_kraken_%A_%a.out
#SBATCH --error=/home/ysanz/CovidPersistente/03_kraken/CovidPersistente_kraken_%A_%a.err
#SBATCH --array=1-39
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task 10
#SBATCH --mem=100G

# Submitted batch job 2219471

#Load modules
module load Kraken2/2.1.2-gompi-2021b
module load Bracken/2.6.2-GCCcore-11.2.0

# take one line of the sample list each time
SEEDFILE=/home/ysanz/CovidPersistente/data/filenames_CovidPersistente.txt # coge la lista de muestras
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//") # me quedo con el nombre de la muestra
WORKDIR=/home/ysanz/CovidPersistente
THREADS=10

# Database Downloads (UPDATE: 6/7/2022)
##cd /home/ysanz/ddbb
##wget --no-check-certificate https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20220607.tar.gz
##mkdir k2_standard_20220607 && tar -xvzf k2_standard_20220607.tar.gz -C /home/ubuntu/myarchive tar -xvzf 

# Taxonomic assignment of metagenomic reads
#mkdir ${WORKDIR}/03_kraken
kraken2 --db /home/ysanz/ddbb/k2_standard_20220607/ \
        --confidence 0.1 \
        --threads ${THREADS} \
        --paired ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_1.fastq.gz ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_2.fastq.gz \
        --output ${WORKDIR}/03_kraken/${SAMPLE}.standard.kraken2 \
        --report ${WORKDIR}/03_kraken/${SAMPLE}.standard.kraken2.report
#Reading a Kraken report
# 1.Percentage of reads covered by the clade rooted at this taxon
# 2.Number of reads covered by the clade rooted at this taxon
# 3.Number of reads assigned directly to this taxon
# 4.A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply ‘-‘.
# 5.NCBI taxonomy ID
# 6.Indented scientific name


# Abundance estimation 
# https://www.nicholas-ollberding.com/post/taxonomic-and-functional-profiling-using-biobakery-workflows/
## Specie
bracken -d /home/ysanz/ddbb/k2_standard_20220607/ \
        -i ${WORKDIR}/03_kraken/${SAMPLE}.standard.kraken2.report \
        -o ${WORKDIR}/03_kraken/${SAMPLE}.standard.bracken_species.txt \
        -w ${WORKDIR}/03_kraken/${SAMPLE}.standard.bracken.report_species.txt \
        -r 75  \
        -t 10 \
        -l S
## Genus
bracken -d /home/ysanz/ddbb/k2_standard_20220607/ \
        -i ${WORKDIR}/03_kraken/${SAMPLE}.standard.kraken2.report \
        -o ${WORKDIR}/03_kraken/${SAMPLE}.standard.bracken_genus.txt \
        -w ${WORKDIR}/03_kraken/${SAMPLE}.standard.bracken.report_genus.txt \
        -r 75  \
        -t 10 \
        -l G

## Phylum
bracken -d /home/ysanz/ddbb/k2_standard_20220607/ \
        -i ${WORKDIR}/03_kraken/${SAMPLE}.standard.kraken2.report \
        -o ${WORKDIR}/03_kraken/${SAMPLE}.standard.bracken_phylum.txt \
        -w ${WORKDIR}/03_kraken/${SAMPLE}.standard.bracken.report_phylum.txt \
        -r 75  \
        -t 10 \
        -l P
