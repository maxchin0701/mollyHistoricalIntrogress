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

#call genotypes for all
gatk GenotypeGVCFs \
	-R "${snakemake_input[ref]}" \
	-V gendb://"${snakemake_params[genDBDir]}" \
	-O "${snakemake_output[0]}"

#call genotypes for parent and hybrids
gatk GenotypeGVCFs \
	-R "${snakemake_input[ref]}" \
	-V gendb://"${snakemake_params[genDBParHybDir]}" \
	-O "${snakemake_output[1]}"

#filter vcf for all
bcftools filter \
	-i 'QUAL>30 & MIN(FORMAT/DP[0-])>5 && Type="snp"' \
	"${snakemake_output[0]}" | \
	bcftools view -v snps -o "${snakemake_output[2]}"

#filter vcf for parent and hybrids
bcftools filter \
	"${snakemake_output[1]}" \
	-i 'QUAL>30 & MIN(FORMAT/DP[0-])>5 && Type="snp"' | \
	bcftools view -m2 -M2 -v snps -o "${snakemake_output[3]}"



