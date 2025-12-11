#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o indexCombVCF

#activate conda env
module load GATK
module load anaconda3/2022.10
conda activate genomics

#index with gatk
gatk IndexFeatureFile \
	-I "${snakemake_input[sampVCF]}"