#!/bin/bash

#SBATCH --job-name=arrayJob
#SBATCH --output=09_humann/arrayJob_%A_%a.out
#SBATCH --error=09_humann/arrayJob_%A_%a.err
#SBATCH --array=3-129
#SBATCH --partition=long
#SBATCH --cpus-per-task 20
#SBATCH --mem=100G

# Submitted batch job 2220106



# Loading the necessary modules
# module load Python/3.9.6-GCCcore-11.2.0
# module load Python/3.6.1-foss-2016b
# module load binutils/2.26-GCCcore-5.4.0
# module load MetaPhlAn3/3.0.13-foss-2016b-Python-3.6.1
# module load DIAMOND/2.0.13-GCC-11.2.0
module load Anaconda3/5.3.0

# Activating conda environment
source activate metagenomics

# take one line of the sample list each time
SEEDFILE=/home/ysanz/CELIACOS/data/samplelist.txt
SAMPLE=$(cat $SEEDFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed "s/_.*//")
WORKDIR=/home/ysanz/CELIACOS
RESULTDIR=${WORKDIR}/09_humann/${SAMPLE}
METAPHLAN_INDEX="mpa_vOct22_CHOCOPhlAnSGB_202212"
BOWTIE2_DB="/home/ysanz/.local/lib/python3.9/site-packages/metaphlan/metaphlan_databases/"
THREADS=20

# Creating result directory
if [[ ! -d "$RESULTDIR" ]]; then
    echo "Creating directories for ${SAMPLE}..."
    mkdir -p "$RESULTDIR"
    mkdir -p "$RESULTDIR/data"
    mkdir -p "$RESULTDIR/results"
    echo "Directories created for ${SAMPLE}."
else
    echo "Directories previously created for ${SAMPLE}."
fi

# Merging the FASTQ.GZ files
if [[ ! -f "$RESULTDIR/data/$SAMPLE_merged.fastq.gz" ]]; then
    echo "Concatenating ${SAMPLE} fastaq.gz files"
    cat ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_1.fastq.gz \
        ${WORKDIR}/02_bowtie2/${SAMPLE}_host_removed_2.fastq.gz > ${RESULTDIR}/data/${SAMPLE}_merged.fastq.gz
    echo "Files concatenated for ${SAMPLE}."
else
    echo "Files already concatenated for ${SAMPLE}."
fi

# Performing the Humann3 pipeline

humann --input ${RESULTDIR}/data/${SAMPLE}_merged.fastq.gz \
    --output  ${RESULTDIR}/results --threads ${THREADS}  \
    --input-format "fastq.gz"  \
    --metaphlan-options "--input_type fastq --nproc $THREADS --index $METAPHLAN_INDEX  --bowtie2db $BOWTIE2_DB" \
    --resume

# TODO: Explain the commands for the pipeline
# metaphlan demo.fastq --input_type fastq -o profiled_metagenome.txt --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --index mpa_vOct22_CHOCOPhlAnSGB_202212  --bowtie2db /home/ysanz/.local/lib/python3.9/site-packages/metaphlan/metaphlan_databases/

# humann --input /home/ysanz/CELIACOS/09_humann/CM01/data/CM01_merged.fastq.gz --output  /home/ysanz/CELIACOS/09_humann/CM01/results --threads 10 --input-format "fastq.gz" --metaphlan-options "--input_type fastq