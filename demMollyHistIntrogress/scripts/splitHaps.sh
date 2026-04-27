#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#source module command
source ~/.bashrc
source /etc/profile.d/modules.sh

#activate conda env
module load tabix/2013-12-16
conda activate R

#split haplotypes
Rscript scripts/splitHapsVCF.R ${snakemake_input[0]} ${snakemake_output[0]}

#unzip and bgzip 
gunzip ${snakemake_output[0]}
bgzip data/splitHap/splitParHyb.vcf
