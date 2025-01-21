#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -o createSampleMap

#activate conda
module load anaconda3/2022.10
conda activate genomics

#save variables
sampN=${snakemake_params[sampNames]}
out=${snakemake_output[0]}

#loop through and build file
for samp in $sampN; do
	filePath=data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz
	printf "%s\t%s\n" "$samp" "$filePath" >> "$out"
done
