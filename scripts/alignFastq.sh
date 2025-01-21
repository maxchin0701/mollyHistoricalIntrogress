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

#create dir for sample if it doesnt already exist
if [ ! -d "data/aligned/"${snakemake_wildcards[sample]}"" ]; then
	mkdir data/aligned/"${snakemake_wildcards[sample]}"
fi

#do alignment
bwa-mem2 mem -t 36 \
	"${snakemake_input[ref]}" \
	"${snakemake_input[forRead]}" \
	"${snakemake_input[revRead]}" > data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\.sam

#convert sam to bam
gatk SamFormatConverter -I data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\.sam \
	-O data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\.bam

#remove sam
rm -r data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\.sam

#add read groups to bam
gatk AddOrReplaceReadGroups \
    I=data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\.bam \
    O=data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\Sorted.bam \
    RGID=group"${snakemake_wildcards[sample]}" \
    RGLB= lib"${snakemake_wildcards[sample]}" \
    RGPL=illumina \
    RGPU=unit"${snakemake_wildcards[sample]}" \
    RGSM="${snakemake_wildcards[sample]}" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

#remove unsorted bam
rm -r data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\.bam

#mark duplicates
gatk MarkDuplicatesSpark \
	-I data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\Sorted.bam \
	-O data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\SortedDupMarked.bam

#remove non dup-marked bam
rm -r data/aligned/"${snakemake_wildcards[sample]}"/"${snakemake_wildcards[sample]}"\Sorted.bam



