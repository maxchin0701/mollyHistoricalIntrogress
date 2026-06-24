#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node=20
#SBATCH -o scriptOuts/runSinger

#load bcftools
module load bcftools

#concatenate
bcftools concat -Oz -o data/mergedHapCalls/parHapCalls.vcf.gz data/hapCalls/chr*\CombinedFiltPar*