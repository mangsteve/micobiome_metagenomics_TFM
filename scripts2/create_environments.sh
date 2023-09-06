conda create --name fastqc-env python=3.10
conda activate fastqc-env
conda install -c bioconda fastqc

conda create --name trimmomatic-env python=3.10
conda activate trimmomatic-env
conda install -c bioconda trimmomatic

conda create --name bowtie2-env python=3.10
conda activate bowtie2-env
conda install -c bioconda bowtie2
conda install -c bioconda samtools

conda create --name bwa-env python=3.10
conda activate bwa-env
conda install -c bioconda bwa bwa-mem2 samtools

conda create --name kraken-env python=3.10
conda activate kraken-env
conda install -c bioconda kraken 
conda install -c bioconda krakentools 

conda create --name bracken-env python=3.10
conda activate bracken-env
conda install -c bioconda bracken

conda create --name krona-env python=3.10
conda activate krona-env
conda install -c bioconda krona krakentools

conda create --name multiqc-env python=3.10
conda activate multiqc-env
conda install -c bioconda multiqc

conda create --name seqtk-env python=3.10
conda activate seqtk-env
conda install -c bioconda seqtk pigz

conda create --name seqkit-env python=3.10
conda activate seqkit-env
conda install -c bioconda seqkit pigz