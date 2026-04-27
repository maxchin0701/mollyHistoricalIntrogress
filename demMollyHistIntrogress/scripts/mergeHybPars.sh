#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node=8
#SBATCH -o scriptOuts/mergeHybParsSeparate

#load conda env
module load bcftools/1.19

#set chromosome
CHR="${snakemake_wildcards[chr]}"

#merge hybrids
bcftools merge -Oz -o "${snakemake_output[0]}" data/phasedHapsWhatshap/${CHR}_Pfo_*\.gz data/phasedHapsWhatshap/${CHR}_Pla_*\.gz data/phasedHapsWhatshap/${CHR}_Pme_*\.gz

#index 
bcftools tabix "${snakemake_output[0]}"

#remove split sample files
rm -r data/phasedHapsWhatshap/${CHR}_Pfo* data/phasedHapsWhatshap/${CHR}_Pla* data/phasedHapsWhatshap/${CHR}_Pme*


