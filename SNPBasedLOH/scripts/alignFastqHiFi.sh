#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=36
#SBATCH -o align

#activate conda env
module load bedtools
module load GATK
module load anaconda3/2022.10
conda activate genomics

#export minimap path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/minimap2

#create dir for sample if it doesnt already exist
if [ ! -d "data/aligned/"${snakemake_wildcards[sample_hifi]}"" ]; then
	mkdir data/aligned/"${snakemake_wildcards[sample_hifi]}"
fi

#do alignment
minimap2 -ax map-hifi \
	-t 35 \
	"${snakemake_input[ref]}" "${snakemake_input[reads]}" > data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\.sam

#convert sam to bam
gatk SamFormatConverter -I data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\.sam \
	-O data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\.bam

#remove sam
rm -r data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\.sam

#add read groups to bam
gatk AddOrReplaceReadGroups \
    I=data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\.bam \
    O=data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\Sorted.bam \
    RGID=group"${snakemake_wildcards[sample_hifi]}" \
    RGLB= lib"${snakemake_wildcards[sample_hifi]}" \
    RGPL=illumina \
    RGPU=unit"${snakemake_wildcards[sample_hifi]}" \
    RGSM="${snakemake_wildcards[sample_hifi]}" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

#remove unsorted bam
rm -r data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\.bam

#mark duplicates
gatk MarkDuplicatesSpark \
	-I data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\Sorted.bam \
	-O data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\SortedDupMarked_hifi.bam

#remove non dup-marked bam
rm -r data/aligned/"${snakemake_wildcards[sample_hifi]}"/"${snakemake_wildcards[sample_hifi]}"\Sorted.bam



