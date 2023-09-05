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

# Importante! los directorios de trabajo donde se crean los output y error tienen que estar creados de antemano.
# Antes de lanzar el script con sbatch sbatch_01trim.sh tengo que hacer: mkdir 01_trimmomatic 
#SBATCH --job-name=CovidPersistente_trim
#SBATCH --output=/home/ysanz/CovidPersistente/01_trimmomatic/CovidPersistente_trim_%A_%a.out
#SBATCH --error=/home/ysanz/CovidPersistente/01_trimmomatic/CovidPersistente_trim_%A_%a.err
#SBATCH --array=1-39
#SBATCH --partition=long
#SBATCH --cpus-per-task 4
#SBATCH --mem=20G

# Submitted batch job 2219606

# load all modules or virtual enviroment you will need.
module load Trimmomatic/0.36-Java-1.8.0_92

# take one line of the sample list each time
SEEDFILE=/home/ysanz/CovidPersistente/data/filenames_CovidPersistente.txt # coge la lista de muestras
FILENAME=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_1.fq.gz//")
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//")

# variables
WORKDIR=/home/ysanz/CovidPersistente
DATADIR=/home/ysanz/CovidPersistente/data
THREADS=20

# Forth step: Launch commands that use the previous variable as ${SAMPLE}
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${THREADS} -phred33 \
-trimlog ${WORKDIR}/01_trimmomatic/${SAMPLE}.trimlog \
${DATADIR}/${FILENAME}_1.fq.gz ${DATADIR}/${FILENAME}_2.fq.gz \
${WORKDIR}/01_trimmomatic/${SAMPLE}_1.trimmed.fastq.gz \
${WORKDIR}/01_trimmomatic/${SAMPLE}_1.trimmed.single.fastq.gz \
${WORKDIR}/01_trimmomatic/${SAMPLE}_2.trimmed.fastq.gz \
${WORKDIR}/01_trimmomatic/${SAMPLE}_2.trimmed.single.fastq.gz \
ILLUMINACLIP:${DATADIR}/Adapters-for-MGI_simpletrimming.txt:2:30:10 \
SLIDINGWINDOW:20:30 \
MINLEN:75

# FastQC
module load FastQC/0.11.8-Java-1.8.0_92
mkdir ${WORKDIR}/01_trimmomatic/00_fastqc_trim
fastqc ${WORKDIR}/01_trimmomatic/${SAMPLE}_1.trimmed.fastq.gz ${WORKDIR}/01_trimmomatic/${SAMPLE}_2.trimmed.fastq.gz --threads ${THREADS} --outdir ${WORKDIR}/01_trimmomatic/00_fastqc


# MultiQC: SIEMPRE A MANO! 
#module load MultiQC/1.7-foss-2016b-Python-2.7.12
#multiqc -p --interactive -o multiqc_report -n multiqc_report_after-trim.html *fastqc.zip

# Para ver como evolucia el trabajo:
#sacct -j JobNumber