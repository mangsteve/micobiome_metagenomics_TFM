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

# Save images in order to transfer to server
sudo docker save -o dockerimage_kraken_with_pigz.tar ccarlos/registry:kraken_with_pigz 
sudo docker save -o dockerimage_seqkit_with_samtools.tar ccarlos/registry:seqkit_with_samtools

# Transfer .tar files via scp and go to the server
cd /DATA12/COMUN/tmp
sudo docker load -i dockerimage_kraken_with_pigz.tar
sudo docker load -i dockerimage_seqkit_with_samtools.tar
