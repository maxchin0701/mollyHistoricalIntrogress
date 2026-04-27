#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node=8
#SBATCH -o scriptOuts/mergeHybParsSeparate

#source module command
source /etc/profile.d/modules.sh

#load conda env
module load bcftools/1.19

#set chromosome
CHR="${snakemake_wildcards[chr]}"

#filter vcf for parent and hybrids
bcftools filter \
	"${snakemake_input[0]}" \
	-i 'QUAL>30 & MIN(FORMAT/GQ[0-])>10 && Type="snp"' | \
	bcftools view -m2 -M2 -v snps -Oz -o "${snakemake_output[0]}"

#index
bcftools tabix "${snakemake_output[0]}"

