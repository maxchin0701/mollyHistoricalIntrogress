#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o downloadSRA

#activate conda env
module load anaconda3/2022.10
conda activate braker

#get into scripts directory
cd scripts

#export path to sratoolkit
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/sratoolkit.3.1.1-centos_linux64/bin

#echo
echo "Downloading ${snakemake_wildcards[srx]}"

#prefetch reads
prefetch "${snakemake_wildcards[srx]}" -O ../data/rawSeq/"${snakemake_wildcards[species]}"

#extract reads
cd ../data/rawSeq/"${snakemake_wildcards[species]}"
fasterq-dump "${snakemake_wildcards[srx]}" --include-technical -S --outdir ./"${snakemake_wildcards[srx]}" --outfile "${snakemake_wildcards[sample]}"

#remove unnecessary files
gzip *.fastq
rm -r ./"${snakemake_wildcards[srx]}"/"$snakemake_wildcards[srx].sra" *.fastq