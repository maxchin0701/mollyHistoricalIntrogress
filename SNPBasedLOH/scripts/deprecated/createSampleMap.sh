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
sampPHN=${snakemake_params[sampNamesParHyb]}
sampParsN=${snakemake_params[sampNamesPars]}
out=${snakemake_output[0]}
outPH=${snakemake_output[1]}
outP=${snakemake_output[2]}


#loop through and build file for all samps
for samp in $sampN; do
	if [ "$samp" = "Psul" ] || [ "$samp" = "Pme_sul" ]; then
		mv data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\-hifi.vcf.gz data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz
		mv data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\-hifi.vcf.gz.tbi data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz.tbi
	fi
	filePath=data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz
	printf "%s\t%s\n" "$samp" "$filePath" >> "$out"
done

#loop through and build file for all samps
for samp in $sampPHN; do
	if [ "$samp" = "Psul" ] || [ "$samp" = "Pme_sul" ]; then
		mv data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\-hifi.vcf.gz data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz
		mv data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\-hifi.vcf.gz.tbi data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz.tbi
	fi
	filePath=data/hapCalls/$samp\/$samp"${snakemake_wildcards[chr]}"\.vcf.gz
	printf "%s\t%s\n" "$samp" "$filePath" >> "$outPH"
done

