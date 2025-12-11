#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 2:00:00
#SBATCH --ntasks-per-node=32
#SBATCH -o maskRepeats-%j

#load conda env
module load anaconda3/2022.10
conda activate repeatMasker

#read in hap
SPECIES=$1
REPDAT=$2

#get into wd
cd ../software/RepeatMasker

#run repeat masker for pFor
./RepeatMasker -dir ../../AGRHistIntrogress/data/maskAssembly/${SPECIES}Mask \
	-lib  ../../AGRHistIntrogress/data/repeatDatabases/${REPDAT} \
	-pa 32 -e rmblast --xsmall ../../AGRHistIntrogress/data/ref/${SPECIES}.fasta
