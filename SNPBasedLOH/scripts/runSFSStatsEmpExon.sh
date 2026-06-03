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
POPS=data/popsHTSPopgen.tsv
VCF=data/mergedHapCalls/parHapCalls.vcf.gz
BED=output/LOHRegionsSharedCombined/LOHRegionsSharedCombinedExon.bed
OUT=output/selectionStats/LOHRegionsEmpStatsExon.tsv

#run SFS stats
sfs --vcf ${VCF} --pops ${POPS} --bed ${BED} --outgroup Pret > ${OUT}