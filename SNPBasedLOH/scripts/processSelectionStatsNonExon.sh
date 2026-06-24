#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=20
#SBATCH -o scriptOuts/runSinger

module load anaconda3/2022.10

#activate R env
conda activate R 

#run scripts
Rscript scripts/processSelectionStatsNonExon.R ${snakemake_input[0]} ${snakemake_output[0]} ${snakemake_output[1]} ${snakemake_output[2]}