#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/getExonsBed

#set variables
GFF=data/funcAnnotate/eMapperOut/funcAnnotatePFor.emapper.decorated.gff
BED=data/funcAnnotate/PForExons.bed

#convert gff exons to bed 
awk 'BEGIN {OFS="\t"}; $3=="exon" {print $1, $4, $5}' data/funcAnnotate/eMapperOut/funcAnnotatePFor.emapper.decorated.gff > data/funcAnnotate/PForExons.bed