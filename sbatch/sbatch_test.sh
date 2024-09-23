#!/bin/bash

#################################################################################################################
## script para el quality control de los ficheros                                                              ##
##                                                                                                             ##
## SBATCH --job-name=name: Nombre del trabajo                                                                  ##
## SBATCH --output=dir/name_%A_%a.out #archivo_%j.out: Directorio y nombre del fichero de salida               ##
## SBATCH --error=dir/name_%A_%a.err #archivo_%j.err: Directorio y nombre del ficherro de errores              ##
## SBATCH --array=1-n: paralelizacion del trabajo, lanza simultaneamente 16 trabajos                           ##
## SBATCH --partition=long: la particion que queremos usar                                                     ##
##          "slimits" para ver las opciones                                                                    ##
## SBATCH --cpus-per-task 4                                                                                    ##
## SBATCH --mem 10G                                                                                            ##
#################################################################################################################

#SBATCH --job-name=test_nf
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --partition=long
#SBATCH --cpus-per-task 4
#SBATCH --mem=2G

module load  Nextflow/23.04.2
module load Anaconda3/5.3.0

nextflow run all.nf -c config/run_samples_testMetaphlan_garnatxa.config -profile conda -resume -with-timeline timeline.html -with-report report.html -with-dag pipeline_dag.html