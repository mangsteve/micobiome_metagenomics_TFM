FROM staphb/kraken2

#Set the working directory to be used when the docker gets run
WORKDIR /usr

# Do a few updates of the base system and install R (via the r-base package)
RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y pigz && \
        apt-get install -y git


## Add Krakentools to kraken2-env
#get scripts from github
RUN git clone https://github.com/jenniferlu717/KrakenTools && \
    krakendir=$(which kraken2) && \
    condadir=$(dirname $krakendir) && \
    echo $condadir && \
    cp -v KrakenTools/*.py $condadir  && \
    chmod 777 $condadir/*.py && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    rm -rf KrakenTools
