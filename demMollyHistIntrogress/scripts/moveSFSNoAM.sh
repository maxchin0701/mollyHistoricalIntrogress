#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pmex

#get model
model=(${snakemake_wildcards[modelNoAM]})

#get into right wd
cd data/${model}

#symlink
ln -s ../sfsNoAM/fastsimcoal2/demModelMolly_jointMAFpop1_0.obs ./${model}_jointMAFpop1_0.obs



