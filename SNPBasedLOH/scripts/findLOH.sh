#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o 

#load conda
module load anaconda3/2022.10

#activate R env
conda activate R 

#run scripts
Rscript scripts/findLOH.R ${snakemake_wildcards[chr]} ${snakemake_wildcards[pars]} ${snakemake_wildcards[gqcut]} ${snakemake_wildcards[dpcut]} ${snakemake_input[0]} 