#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 6:00:00
#SBATCH --ntasks-per-node=6
#SBATCH -o callHaplotypes

#activate conda env
module load anaconda3/2022.10
conda activate gatk

#export gatk path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/gatk-4.1.9.0

#create dir for sample if it doesnt already exist
if [ ! -d "data/hapCalls/"${snakemake_wildcards[sample_hifi]}"" ]; then
	mkdir data/hapCalls/"${snakemake_wildcards[sample_hifi]}"
fi

#call variants (save in temp file)
gatk --java-options "-Xms80g -Xmx80g -XX:ParallelGCThreads=2" HaplotypeCaller \
	-R  "${snakemake_input[ref]}"\
	-I "${snakemake_input[bam]}" \
	-L "${snakemake_input[intervals]}" \
	-O data/hapCalls/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}""${snakemake_wildcards[chr]}"_TEMP.vcf.gz \
	--native-pair-hmm-threads 2 \
	--pair-hmm-implementation LOGLESS_CACHING \
	-ERC GVCF

#move to final
mv data/hapCalls/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}""${snakemake_wildcards[chr]}"_TEMP.vcf.gz "${snakemake_output[0]}"
mv data/hapCalls/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}""${snakemake_wildcards[chr]}"_TEMP.vcf.gz.tbi "${snakemake_output[0]}".tbi

