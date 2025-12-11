#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o genIntervalList

module load GATK
module load anaconda3/2022.10
conda activate genomics

#run python script to generate bed for chr
python scripts/createBed.py --inFasta "${snakemake_input[0]}" --scaff "${snakemake_wildcards[chr]}" --intSize 1000 --outBed data/refGenome/intervalLists/PFor_"${snakemake_wildcards[chr]}"\.bed

#generate interval_list from bed
gatk BedToIntervalList I=data/refGenome/intervalLists/PFor_"${snakemake_wildcards[chr]}"\.bed \
	O=data/refGenome/intervalLists/PFor_"${snakemake_wildcards[chr]}"\.interval_list \
	SD="${snakemake_input[1]}"



