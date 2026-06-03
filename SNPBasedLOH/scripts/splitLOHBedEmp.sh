#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/getExonsBed

#load bedtools
module load bedtools

#set variables
EXONS=data/funcAnnotate/PForExons.bed
LOH=output/LOHRegionsSharedCombined/LOHRegionsSharedCombined.bed
LOHNOEXONS=output/LOHRegionsSharedCombined/LOHRegionsSharedCombinedNonExon.bed
LOHEXONS=output/LOHRegionsSharedCombined/LOHRegionsSharedCombinedExon.bed

#get non-exon bed 
bedtools subtract -a $LOH -b $EXONS > $LOHNOEXONS

#get exon bed
bedtools subtract -a $LOH -b $LOHNOEXONS > $LOHEXONS
