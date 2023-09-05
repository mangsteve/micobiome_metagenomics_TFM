###############################################################################
##  Crear todos los directorios de trabajo para un analisis de metagenomica  ##
##           El nombre del proyecto (project) sera el del WORKDIR            ##
###############################################################################

## Completar el nombre del proyecto
project=CovidPersistente

WORKDIR=/home/ysanz/${project}

## Crear todos los directorios
# fastqc rawdata
mkdir ${WORKDIR}/00_fastqc
# trimmomatic
mkdir ${WORKDIR}/01_trimmomatic
mkdir ${WORKDIR}/01_trimmomatic/00_fastqc
# bowtie2
mkdir ${WORKDIR}/02_bowtie2
mkdir ${WORKDIR}/02_bowtie2/00_fastqc
# kraken & bracken
mkdir ${WORKDIR}/03_kraken
mkdir ${WORKDIR}/04_bracken_abundance_file