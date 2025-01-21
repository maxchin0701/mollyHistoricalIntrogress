#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=20
#SBATCH -o indexRef

#activate conda env
module load GATK
module load samtools
module load anaconda3/2022.10
conda activate genomics

#index genome file (bwa)
bwa-mem2 index "${snakemake_input[0]}"

#index genome (samtools)
samtools faidx "${snakemake_input[0]}"

#generate dictionary 
gatk CreateSequenceDictionary \
		R="${snakemake_input[0]}" \
		O="${snakemake_output[1]}"
