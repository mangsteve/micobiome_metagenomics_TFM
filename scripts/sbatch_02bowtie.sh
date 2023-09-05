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
## SBATCH --cpus-per-task 4                                                                                    ##
## SBATCH --mem 10G                                                                                            ##
#################################################################################################################

#SBATCH --job-name=CovidPersistente_bowtie
#SBATCH --output=/home/ysanz/CovidPersistente/02_bowtie2/CovidPersistente_bowtie_%A_%a.out
#SBATCH --error=/home/ysanz/CovidPersistente/02_bowtie2/CovidPersistente_bowtie_%A_%a.err
#SBATCH --array=1-39
#SBATCH --partition=long
#SBATCH --cpus-per-task 4
#SBATCH --mem=20G

module load OpenBLAS/0.3.18-GCC-11.2.0
module load Bowtie2/2.4.5-foss-2021b
module load SAMtools/1.14-GCC-11.2.0

# Submitted batch job 2219392

# take one line of the sample list each time
SEEDFILE=/home/ysanz/CovidPersistente/data/filenames_CovidPersistente.txt # coge la lista de muestras
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//") # me quedo con el nombre de la muestra
WORKDIR=/home/ysanz/CovidPersistente
DATADIR=/home/ysanz/CovidPersistente/01_trimmomatic
DDBBDIR=/home/ysanz/ddbb/bowtie-index
THREADS=10

# a) bowtie2 mapping against host sequence
## bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
bowtie2 -x ${DDBBDIR}/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index --threads ${THREADS} --fast-local --phred33 \
-1 ${WORKDIR}/01_trimmomatic/${SAMPLE}_1.trimmed.fastq.gz -2 ${WORKDIR}/01_trimmomatic/${SAMPLE}_2.trimmed.fastq.gz \
-S ${WORKDIR}/02_bowtie2/${SAMPLE}_mapped_and_unmapped.sam 

## convert file .sam to .bam
samtools view -bS ${WORKDIR}/02_bowtie2/${SAMPLE}_mapped_and_unmapped.sam > ${WORKDIR}/02_bowtie2/${SAMPLE}_mapped_and_unmapped.bam
rm ${WORKDIR}/02_bowtie2/${SAMPLE}_mapped_and_unmapped.sam

# b) filter required unmapped reads
## SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
##  -f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
##  -F 256   Do not(-F) extract alignments which are: <not primary alignment>
samtools view -b -f 12 -F 256 ${WORKDIR}/02_bowtie2/${SAMPLE}_mapped_and_unmapped.bam > ${WORKDIR}/02_bowtie2/${SAMPLE}_bothReadsUnmapped.bam
rm ${WORKDIR}/02_bowtie2/${SAMPLE}_mapped_and_unmapped.bam

# c)  split paired-end reads into separated fastq files .._R1 .._R2
## sort bam file by read name (-n) to have paired reads next to each other (3 parallel threads, each using up to 5G memory)
samtools sort -n -m 5G -@ 3 ${WORKDIR}/02_bowtie2/${SAMPLE}_bothReadsUnmapped.bam -o ${WORKDIR}/02_bowtie2/${SAMPLE}_bothReadsUnmapped_sorted.bam
rm ${WORKDIR}/02_bowtie2/${SAMPLE}_bothReadsUnmapped.bam

#c para la compresion usar 7
samtools fastq -c 7 -@ ${THREADS} ${WORKDIR}/02_bowtie2/${SAMPLE}_bothReadsUnmapped_sorted.bam \
    -1 ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_1.fastq.gz \
    -2 ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_2.fastq.gz \
    -0 /dev/null -s /dev/null -n

# FastQC
module load FastQC/0.11.8-Java-1.8.0_92
mkdir ${WORKDIR}/02_bowtie2/00_fastqc
fastqc ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_1.fastq.gz ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_2.fastq.gz \
--threads ${THREADS} --outdir ${WORKDIR}/02_bowtie2/00_fastqc

# MultiQC
#module load MultiQC/1.7-foss-2016b-Python-2.7.12
#multiqc -p --interactive -o multiqc_report_Bowtie2 -n multiqc_report_after-removehost.html *fastqc.zip
