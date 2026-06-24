#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 6:00:00
#SBATCH --ntasks-per-node=6
#SBATCH -o callHaplotypes

#activate conda env
module load GATK
module load anaconda3/2022.10
conda activate gatk

#create dir for sample if it doesnt already exist
if [ ! -d "data/hapCalls/"${snakemake_wildcards[sample]}"" ]; then
	mkdir data/hapCalls/"${snakemake_wildcards[sample]}"
fi

#create dir for sample if it doesnt already exist
'''
if [ ! -d "data/realignedBAM/"${snakemake_wildcards[sample]}"" ]; then
	mkdir data/realignedBAM/"${snakemake_wildcards[sample]}"
fi
'''

#call variants (save in temp file)
gatk --java-options "-Xms40G -Xmx40G -XX:ParallelGCThreads=2" HaplotypeCaller \
	-R  "${snakemake_input[ref]}"\
	-I "${snakemake_input[bam]}" \
	-L "${snakemake_input[intervals]}" \
	-O data/hapCalls/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}""${snakemake_wildcards[chr]}"_TEMP.vcf.gz \
	-ERC GVCF \
	-bamout data/realignedBAM/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}""${snakemake_wildcards[chr]}"_realigned.bam

#move to final
mv data/hapCalls/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}""${snakemake_wildcards[chr]}"_TEMP.vcf.gz "${snakemake_output[0]}"
mv data/hapCalls/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}""${snakemake_wildcards[chr]}"_TEMP.vcf.gz.tbi "${snakemake_output[0]}".tbi