#!/bin/bash

#################################################################################################################
## script para el quality control de los ficheros                                                              ##
##                                                                                                             ##
## SBATCH --job-name=name: Nombre del trabajo                                                                  ##
## SBATCH --output=dir/name_%A_%a.out #archivo_%j.out: Directorio y nombre del fichero de salida               ##
## SBATCH --error=dir/name_%A_%a.err #archivo_%j.err: Directorio y nombre del ficherro de errores              ##
## SBATCH --array=1-n: paralelizacion del trabajo, lanza simultaneamente 16 trabajos                           ##
##                                      OJO!!! son el numero de muestras no de FICHEROS                        ##
## SBATCH --partition=long: la particion que queremos usar                                                     ##
##          "slimits" para ver las opciones                                                                    ##
## SBATCH --cpus-per-task 4                                                                                    ##
## SBATCH --mem 10G                                                                                            ##
#################################################################################################################

#SBATCH --job-name=CovidPersistente_QC
#SBATCH --output=/home/ysanz/CovidPersistente/00_fastqc/CovidPersistente_QC%A_%a.out #archivo_%j.out
#SBATCH --error=/home/ysanz/CovidPersistente/00_fastqc/CovidPersistente_QC%A_%a.err #archivo_%j.err
#SBATCH --array=1-39 
#SBATCH --partition=long
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G

# Submitted batch job 2219113

module load FastQC/0.11.8-Java-1.8.0_92

# Lista de ficheros
cd /home/ysanz/CovidPersistente/data/
ls *1.fq.gz > filenames_CovidPersistente.txt

# take one line of the sample list each time
SEEDFILE=/home/ysanz/CovidPersistente/data/filenames_CovidPersistente.txt
FILENAME=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_1.fq.gz//")
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//")


WORKDIR=/home/ysanz/CovidPersistente
DATADIR=/home/ysanz/CovidPersistente/data
THREADS=10

#mkdir -p ${WORKDIR}/00_fastqc
fastqc ${DATADIR}/${FILENAME}_1.fq.gz ${DATADIR}/${FILENAME}_2.fq.gz --threads ${THREADS} --outdir ${WORKDIR}/00_fastqc 

# MultiQC
#module load MultiQC/1.7-foss-2016b-Python-2.7.12
#multiqc -p --interactive -o multiqc_report_CovidPersistente -n multiqc_report_CovidPersistente.html *fastqc.zip