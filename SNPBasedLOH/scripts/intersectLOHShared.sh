#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 30:00:00
#SBATCH --ntasks-per-node=36
#SBATCH -o align

#load modules
module load bedtools

#set output
LOHBed=${snakemake_input[0]}
GFF=${snakemake_input[1]}
LOHBedSorted=${snakemake_output[0]}
LOHIntersect=${snakemake_output[1]}

#sort regions
sort -k1,1 -k2,2n ${LOHBed} > ${LOHBedSorted}

#bedtools intersect
bedtools intersect -a ${LOHBedSorted} \
	-b ${GFF} -wa -wb -sorted \
	-g ../../data/PForHapMex.genome > ${LOHIntersect}
