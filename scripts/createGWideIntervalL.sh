#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o createGWideIntervalL

module load GATK
module load anaconda3/2022.10
conda activate genomics

#run python script to generate bed for chr
cat ${snakemake_input[bedF]} > "${snakemake_output[0]}"

#generate interval_list from bed
gatk BedToIntervalList I="${snakemake_output[0]}" \
	O="${snakemake_output[1]}" \
	SD=data/refGenome/PForSch01.dict



