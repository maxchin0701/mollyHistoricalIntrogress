#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=4
#SBATCH -o scriptOuts/getExonsBed

#load bedtools
module load anaconda3/2022.10

#activate R env
conda activate python3

#run script
python sumCDSLOHIntersect.py --inIntersect ${snakemake_input[0]} \
	--inSimIntersect output/LOHRegionsSimIntersectOut \
	--out ${snakemake_output[0]}