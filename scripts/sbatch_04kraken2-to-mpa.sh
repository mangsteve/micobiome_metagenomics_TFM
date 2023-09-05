#!/bin/bash

# Submitted batch job  (species)
# Submitted batch job  (genus & phylum)
# Submitted batch job  (04_bracken_abundace)

# Generating combined abundance tables in mpa format
mkdir /home/ysanz/CovidPersistente/03_kraken/species/mpa
mkdir /home/ysanz/CovidPersistente/03_kraken/genus/mpa
mkdir /home/ysanz/CovidPersistente/03_kraken/phylum/mpa

# Variables
species_mpa=/home/ysanz/CovidPersistente/03_kraken/species/mpa/
genus_mpa=/home/ysanz/CovidPersistente/03_kraken/genus/mpa/
phylum_mpa=/home/ysanz/CovidPersistente/03_kraken/phylum/mpa/

## SPECIES
for i in /home/ysanz/CovidPersistente/03_kraken/report_species/*.standard.bracken.report_species.txt
do
  filename=$(basename "$i")
  fname="${filename%.standard.bracken.report_species.txt}"
  python /home/ysanz/scripts/KrakenTools/kreport2mpa.py -r $i -o ${species_mpa}${fname}_mpa.txt --display-header
done

mkdir ${species_mpa}combined
python /home/ysanz/scripts/KrakenTools/combine_mpa.py -i ${species_mpa}*_mpa.txt -o ${species_mpa}combined/combined_species_mpa.txt
grep -E "(s__)|(#Classification)" ${species_mpa}combined/combined_species_mpa.txt > ${species_mpa}combined/bracken_abundance_species_mpa.txt

## GENUS
for i in /home/ysanz/CovidPersistente/03_kraken/report_genus/*.standard.bracken.report_genus.txt
do
  filename=$(basename "$i")
  fname="${filename%.standard.bracken.report_genus.txt}"
  python /home/ysanz/scripts/KrakenTools/kreport2mpa.py -r $i -o ${genus_mpa}${fname}_mpa.txt --display-header
done

mkdir ${genus_mpa}combined
python /home/ysanz/scripts/KrakenTools/combine_mpa.py -i ${genus_mpa}*_mpa.txt -o ${genus_mpa}combined/combined_genus_mpa.txt
grep -E "(g__)|(#Classification)" ${genus_mpa}combined/combined_genus_mpa.txt > ${genus_mpa}combined/bracken_abundance_genus_mpa.txt

## PHYLUM
for i in /home/ysanz/CovidPersistente/03_kraken/report_phylum/*.standard.bracken.report_phylum.txt
do
  filename=$(basename "$i")
  fname="${filename%.standard.bracken.report_phylum.txt}"
  python /home/ysanz/scripts/KrakenTools/kreport2mpa.py -r $i -o ${phylum_mpa}${fname}_mpa.txt --display-header
done

mkdir ${phylum_mpa}combined
python /home/ysanz/scripts/KrakenTools/combine_mpa.py -i ${phylum_mpa}*_mpa.txt -o ${phylum_mpa}combined/combined_phylum_mpa.txt
grep -E "(p__)|(#Classification)" ${phylum_mpa}combined/combined_phylum_mpa.txt > ${phylum_mpa}combined/bracken_abundance_phylum_mpa.txt

#Cleaning up sample names
sed -i -e 's/.standard.bracken.report_species.txt//g' ${species_mpa}combined/bracken_abundance_species_mpa.txt
sed -i -e 's/.standard.bracken.report_genus.txt//g' ${genus_mpa}combined/bracken_abundance_genus_mpa.txt
sed -i -e 's/.standard.bracken.report_phylum.txt//g' ${phylum_mpa}combined/bracken_abundance_phylum_mpa.txt


#Cleaning up top-level folders
mkdir /home/ysanz/CovidPersistente/04_bracken_abundance_files
cp ${species_mpa}combined/bracken_abundance_species_mpa.txt /home/ysanz/CovidPersistente/04_bracken_abundance_files/.
cp ${genus_mpa}combined/bracken_abundance_genus_mpa.txt /home/ysanz/CovidPersistente/04_bracken_abundance_files/.
cp ${phylum_mpa}combined/bracken_abundance_phylum_mpa.txt /home/ysanz/CovidPersistente/04_bracken_abundance_files/.