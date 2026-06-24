#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=20
#SBATCH -o scriptOuts/runSinger

#initialize conda
module load anaconda3/2022.10
conda activate hts_popgen

#export to path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/hts_popgen

#variables
POPS=${snakemake_input[0]}
VCF=${snakemake_input[1]}
BED=${snakemake_input[2]}
OUT=${snakemake_output[0]}

#run SFS stats
sfs --vcf ${VCF} --pops ${POPS} --bed ${BED} --outgroup Pret > ${OUT}