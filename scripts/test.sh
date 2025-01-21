#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 2:00:00
#SBATCH --ntasks-per-node=10
#SBATCH -o test

#activate conda env
module load anaconda3/2022.10
conda activate braker

#export path to sratoolkit
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/sratoolkit.3.1.1-centos_linux64/bin

species=$1
srx=$2
sample=$3

#extract reads
cd ../data/rawSeq/"$species"
fasterq-dump "$2" --include-technical -S --outdir ./"$2" --outfile "$3"

#remove unnecessary files
gzip ./"$srx"/"$sample"*\.fastq
rm -r ./"$srx"/"$srx".sra ./"$srx"/"$sample"*\.fastq