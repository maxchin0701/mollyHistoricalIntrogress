#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/getExonsBed

#load bedtools
module load bedtools

#set variables
EXONS=${snakemake_input[0]}
LOH=${snakemake_input[1]}
LOHNOEXONS=${snakemake_output[0]}
LOHEXONS=${snakemake_output[1]}

#get non-exon bed 
bedtools subtract -a $LOH -b $EXONS > $LOHNOEXONS

#get exon bed
bedtools subtract -a $LOH -b $LOHNOEXONS > $LOHEXONS
