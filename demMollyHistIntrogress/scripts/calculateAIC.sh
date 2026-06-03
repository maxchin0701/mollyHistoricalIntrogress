#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o calcAIC

#source module command
source /etc/profile.d/modules.sh
source ~/.bashrc

#activate conda environment
conda activate R 

#run Rscript
Rscript scripts/calculateAIC.R ${snakemake_wildcards[model]}
