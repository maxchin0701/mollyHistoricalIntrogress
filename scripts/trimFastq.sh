#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o trimFastq

#activate conda env
module load anaconda3/2022.10
conda activate genomics

#define vars
inFor=(${snakemake_input[forRead]})
inRev=(${snakemake_input[revRead]})
outFor=(${snakemake_output[0]})
outRev=(${snakemake_output[1]})

#create dir for sample if it doesnt already exist
if [ ! -d "data/trimSeq/"${snakemake_wildcards[sample]}"" ]; then
	mkdir data/trimSeq/"${snakemake_wildcards[sample]}"
fi

#run trimmomatic
trimmomatic PE -threads 4 -phred33 $inFor $inRev \
	$outFor data/trimSeq/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\_forward_unpaired.fq.gz \
	$outRev data/trimSeq/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\_reverse_unpaired.fq.gz \
	ILLUMINACLIP:data/trimSeq/TruSeq3-PE-2.fa:2:30:7 LEADING:3 TRAILING:3 MINLEN:50





