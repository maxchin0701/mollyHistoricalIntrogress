#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getPartitions

#load bcftools
source /etc/profile.d/modules.sh
module load bcftools

#get iteration
ITER=(${snakemake_wildcards[iter]})

#split
bcftools +prune -n 1 -w 10kb -N "rand" --random-seed ${ITER} -O z -o ${snakemake_output[0]} ${snakemake_input[0]}

