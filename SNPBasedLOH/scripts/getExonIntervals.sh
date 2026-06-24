#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/getExonsBed

#set variables
GFF=${snakemake_input[0]}
BED=${snakemake_output[0]}

#convert gff exons to bed 
awk 'BEGIN {OFS="\t"}; $3=="exon" {print $1, $4, $5}' ${GFF} > ${BED}