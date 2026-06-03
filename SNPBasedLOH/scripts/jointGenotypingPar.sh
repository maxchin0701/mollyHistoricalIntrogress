#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o jointGenotyping

#activate conda env
module load GATK
module load anaconda3/2022.10
conda activate genomics
module load bcftools

#call genotypes for parent and hybrids
gatk GenotypeGVCFs \
	-R "${snakemake_input[ref]}" \
	-V gendb://"${snakemake_params[genDBParDir]}" \
	-O "${snakemake_output[0]}"

#filter vcf for parent and hybrids
bcftools filter \
	"${snakemake_output[0]}" \
	-i 'QUAL>30 & MIN(FORMAT/DP[0-])>5 & MIN(FORMAT/GQ[0-])>10 && Type="snp"' | \
	bcftools view -m2 -M2 -v snps -o "${snakemake_output[1]}"

