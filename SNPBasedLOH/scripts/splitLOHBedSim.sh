#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/getExonsBed

#load bedtools
module load bedtools

for i in {1..1000}; do
	echo "Processing iteration ${i}"
	EXONS=data/funcAnnotate/PForExons.bed
	LOH=output/LOHRegionsSharedSim/simLOHShared_${i}Sorted.bed
	LOHNOEXONS=output/LOHRegionsSharedSim/simLOHSharedNonExon_${i}Sorted.bed
	LOHEXONS=output/LOHRegionsSharedSim/simLOHSharedExon_${i}Sorted.bed

	#get non-exon bed 
	bedtools subtract -a $LOH -b $EXONS > $LOHNOEXONS

	#get exon bed
	bedtools subtract -a $LOH -b $LOHNOEXONS > $LOHEXONS
done