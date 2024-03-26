FROM staphb/seqkit
WORKDIR /usr/

RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y gcc make libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev && \
        apt-get install -y bzip2

RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    bzip2 -d samtools-1.19.2.tar.bz2 && \
    tar -xvf samtools-1.19.2.tar
    
WORKDIR /usr/samtools-1.19.2   

RUN ./configure --prefix=/usr/local/bin/samtools && \
    make && \
    make install 

ENV PATH="$PATH:/usr/local/bin/samtools/bin"
    
#ln -s /usr/local/bin/samtools/samtools /usr/local/bin/samtools

WORKDIR /usr/

RUN rm -rf samtools-1.19.2