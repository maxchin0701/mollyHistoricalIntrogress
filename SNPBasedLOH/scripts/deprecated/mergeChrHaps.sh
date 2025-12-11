#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o mergeChrHaps

module load bcftools
module load anaconda3/2022.10
conda activate genomics

#merge
bcftools concat ${snakemake_input[chrVCF]} \
	-o "${snakemake_output[0]}" \
	-O z
