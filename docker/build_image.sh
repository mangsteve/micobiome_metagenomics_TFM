## Set up local registry
#sudo docker run -d -p 5000:5000 --restart always --name registry registry:2
# sudo vim /etc/docker/daemon.json
# { "insecure-registries": ["ccarlos.local:5000"] }

sudo docker build -t ccarlos/registry:kraken_with_pigz -f kraken_with_pigz.dockerfile .
sudo docker build -t ccarlos/registry:seqkit_with_samtools -f seqkit_and_samtools.dockerfile .


# Enter to test
sudo docker run -ti --entrypoint /bin/bash ccarlos/registry:kraken_with_pigz

# Enter and mount your home to the container
sudo docker run -ti --mount type=bind,source=/home,target=/home  --entrypoint /bin/bash nanozoo/krona

