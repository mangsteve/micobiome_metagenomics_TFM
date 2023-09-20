#!/bin/bash

#SBATCH --job-name=arrayJob
#SBATCH --output=examples/arrayJob_%A_%a.out
#SBATCH --error=examples/arrayJob_%A_%a.err
#SBATCH --array=1
#SBATCH --partition=long
#SBATCH --cpus-per-task 20
#SBATCH --mem=100G

# 20230814 - Submitted batch job 2220110 2220111



# modules

module load Anaconda3/5.3.0

# Activating conda environment
source activate metagenomics


WORKDIR=/home/ysanz/CELIACOS/examples
#RESULTDIR=${WORKDIR}/09_humann/${SAMPLE}
THREADS=20


# Performing the Humann3 pipeline
echo "Initializing humann"
humann --input ${WORKDIR}/demo.fastq --output  ${WORKDIR}/results/fastq --threads ${THREADS} --input-format "fastq" --metaphlan-options "--input_type fastq --nproc ${THREADS} --index mpa_vOct22_CHOCOPhlAnSGB_202212  --bowtie2db /home/ysanz/.local/lib/python3.9/site-packages/metaphlan/metaphlan_databases/" --resume
