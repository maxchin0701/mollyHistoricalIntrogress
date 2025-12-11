#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 10:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o collectReadCounts

#activate conda env
module load GATK
module load anaconda3/2022.10
conda activate genomics

#take in snakemake vars
bam="${snakemake_input[bam]}"
intervals="${snakemake_input[intervals]}"
outFile="${snakemake_output[0]}"

#create dir for sample if it doesnt already exist
if [ ! -d "data/gCNV/readCounts/"${snakemake_wildcards[sample]}"" ]; then
	mkdir data/gCNV/readCounts/"${snakemake_wildcards[sample]}"
fi

#collect read counts
gatk CollectReadCounts \
		-I $bam \
		-L $intervals \
		--interval-merging-rule OVERLAPPING_ONLY \
		-O $outFile


