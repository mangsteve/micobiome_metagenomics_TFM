#!/bin/bash
#SBATCH --job-name=arrayJob
#SBATCH --output=06_prokka/arrayJob_%A_%a.out
#SBATCH --error=06_prokka/arrayJob_%A_%a.err
#SBATCH --array=1-129
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task 20
#SBATCH --mem=100G

# Load all the modules or virutal enviroments that you will need
module load Anaconda3/5.3.0

# Initializing Anaconda Virtual Environment (metagenomics)
source activate metagenomics

# Take one line of the sample list each time and establish 
SEEDFILE=/home/ysanz/CELIACOS/data/samplelist.txt
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//")
WORKDIR=/home/ysanz/CELIACOS
THREADS=20

# mkdir ${WORKDIR}/06_prokka/${SAMPLE}
# We get into the directory with the scaffolds.fasta file
cd ${WORKDIR}/05_spades/${SAMPLE}

metaprokka --cpus ${THREADS} --metagenome --force --addgenes --addmrna --gcode 11 --centre IATA --compliant \
                --locustag ${SAMPLE} --prefix ${SAMPLE} --outdir ${WORKDIR}/06_prokka/${SAMPLE} \
                ${SAMPLE}_scaffolds.fasta

