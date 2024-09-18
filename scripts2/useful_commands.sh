
squeue # procesos de todos
squeue -u <usuario>

sbatch <fichero_sbatch>
scancel <num proceso> # cancelar un proceso
scancel -u <usuario> #cancelar todos los procesos de un usuario
sacct -j <num proceso> # informacion detallada de proceso


git clone https://github.com/mangsteve/micobiome_metagenomics_TFM

# See branches
git branch -a

# In the server, swhitch to remote branch to run current development version
git switch -c incluir_metaphlan  origin/feature/incluir_metaphlan

git add <ficheros> 
git commit -m "mensaje del commit"
git push origin master # hacerlo desde gitkraken
aaaaaaaaaaaaaaaaaaa