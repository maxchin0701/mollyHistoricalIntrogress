#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o trimFastq

#activate conda env
module load BLAST
module load anaconda3/2022.10
conda activate genomics

#export hifiadapterfilt scripts
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/HiFiAdapterFilt
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/HiFiAdapterFilt/DB

#define vars
inReads=(${snakemake_input[reads]})
outReads=(${snakemake_output[0]})

#create dir for sample if it doesnt already exist
if [ ! -d "data/trimSeq/"${snakemake_wildcards[sample_hifi]}"" ]; then
	mkdir data/trimSeq/"${snakemake_wildcards[sample_hifi]}"
fi

#change to working directory
cd data/rawSeq/"${snakemake_wildcards[sample_hifi]}"

#run hifiadapterfilt
bash hifiadapterfilt.sh -p "${snakemake_wildcards[sample_hifi]}"_Hifi -m 95 -o ../../trimSeq/"${snakemake_wildcards[sample_hifi]}"

#rename filt fastq
mv ../../trimSeq/"${snakemake_wildcards[sample_hifi]}"\/"${snakemake_wildcards[sample_hifi]}"_Hifi.filt.fastq.gz \
	../../trimSeq/"${snakemake_wildcards[sample_hifi]}"\/"${snakemake_wildcards[sample_hifi]}"\Filt.fastq.gz


