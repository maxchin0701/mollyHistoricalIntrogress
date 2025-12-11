#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1:30:00
#SBATCH --ntasks-per-node=4
#SBATCH -o getSFS

#pop0=Plat
#pop1=Pfor
#pop2=Pmex

#activate conda env
module load anaconda3/2022.10
conda activate easySFS

#export path
export PATH=$PATH:/ocean/projects/bio230047p/mchin/software/easySFS

#get population
#pop=(${snakemake_wildcards[popCombo]})
inVCF=(${snakemake_input[inVCF]})
pops=(${snakemake_input[pops]})
nSampPop0=(${snakemake_params[nSampPop0]})
nSampPop1=(${snakemake_params[nSampPop1]})
nSampPop2=(${snakemake_params[nSampPop2]})


#get real sfs
easySFS.py -i ${inVCF} -p ${pops} --proj ${nSampPop0},${nSampPop1},${nSampPop2} -a --prefix demModelMolly -o data/sfs -f --order pop0,pop1,pop2

#rename
#mv data/sfs/${popCombo}/jointMAF${popCombo}.sfs data/sfs/${popCombo}/jointMAF${popCombo}.obs



