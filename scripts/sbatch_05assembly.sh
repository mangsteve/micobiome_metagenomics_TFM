#!/bin/bash
#SBATCH --job-name=arrayJob
#SBATCH --output=05_spades/arrayJob_%A_%a.out
#SBATCH --error=05_spades/arrayJob_%A_%a.err
#SBATCH --array=6-129 # We continue with the remaining samples, since we have completed the first 5 in the testing procedure
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task 10
#SBATCH --mem=100G

# Load all the modules or virutal enviroments that you will need
module load SPAdes/3.15.3-GCC-11.2.0

# Take one line of the sample list each time and establish the sample name in order to relate it with the fastqc
SEEDFILE=/home/ysanz/CELIACOS/data/samplelist.txt
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//")
WORKDIR=/home/ysanz/CELIACOS
THREADS=10


mkdir ${WORKDIR}/02_bowtie2/${SAMPLE}
mkdir ${WORKDIR}/05_spades/${SAMPLE}

metaspades.py -t ${THREADS} -m 300 --only-assembler \
	-1 ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_1.fastq.gz \
	-2 ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_2.fastq.gz \
	-o ${WORKDIR}/02_bowtie2/${SAMPLE}

mv ${WORKDIR}/02_bowtie2/${SAMPLE}/scaffolds.fasta ${WORKDIR}/05_spades/${SAMPLE}/${SAMPLE}_scaffolds.fasta
rm -r ${WORKDIR}/02_bowtie2/${SAMPLE}
